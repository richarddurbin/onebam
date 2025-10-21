/*  File: onebamhts.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 22 00:30 2025 (rd109)
 * Created: Wed Jul  2 13:39:53 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"
#include "sam.h"
#include "array.h"
#include "seqio.h"

typedef struct {
  samFile *f ;
  sam_hdr_t *h ;
  bam1_t *b ;
} BamFile ;

static char *cramReference = 0 ;
void setCramReference (char *cramRef) { cramReference = strdup (cramRef) ; }

static void bamFileClose (BamFile *bf)
{ if (bf->b) bam_destroy1 (bf->b) ;
  if (bf->h) sam_hdr_destroy (bf->h) ;
  if (bf->f) sam_close (bf->f) ;
  free (bf) ;
}

static BamFile *bamFileOpenRead (char* filename)
{
  BamFile *bf = new0 (1, BamFile) ;

  bf->f = sam_open (filename, "r") ;
  if (!bf->f) return false ;

  uint32_t rf = SAM_FLAG | SAM_SEQ ;
  //  if (si->isQual) rf |= SAM_QUAL ;
  if (hts_set_opt (bf->f, CRAM_OPT_REQUIRED_FIELDS, rf))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_REQUIRED_FIELDS\n") ;
      bamFileClose (bf) ;
      return 0 ;
    }
  if (hts_set_opt (bf->f, CRAM_OPT_DECODE_MD, 0))
    { fprintf (stderr, "BamFileOpen failed to set CRAM_OPT_DECODE_MD\n") ;
      bamFileClose (bf) ;
      return 0 ;
    }
  if (cramReference && hts_set_fai_filename (bf->f, cramReference) != 0)
    { fprintf (stderr, "BamFileOpen failed to set cram reference from %s\n", cramReference) ;
      bamFileClose (bf) ;
      return 0 ;
    }

  if (!(bf->h = sam_hdr_read (bf->f)))
    { fprintf(stderr, "BamFileOpen failed to read header\n") ;
      bamFileClose (bf) ;
      return 0 ;
    }

  if (!(bf->b = bam_init1 ()))
    { fprintf(stderr, "BamFileOpen failed to create record buffer\n") ;
      bamFileClose (bf) ;
      return 0 ;
    }

  return bf ;
}

/**************************/

typedef struct {
  char  oneCode ;
  char  tag[3] ;
  char  fmt ;
} AuxField ;
static Array auxFields  = 0 ;
static char *auxType[128] ;

void auxAdd (char *oneCode, char *spec)
{
  if (!auxFields)
    { auxFields = arrayCreate (32, AuxField) ;
      int i ; for (i = 0 ; i < 128 ; ++i) auxType[i] = 0 ;
      auxType['A'] = "CHAR" ;
      auxType['C'] = auxType['c'] = auxType['S'] = auxType['s'] = auxType['I'] = auxType['i'] = "INT" ;
      auxType['f'] = "FLOAT" ;
      auxType['H'] = auxType['Z'] = "STRING" ;
      auxType['D'] = auxType['d'] = auxType['T'] = auxType['t'] = auxType['J'] = auxType['j'] = "INT_LIST" ;
      auxType['g'] = "FLOAT_ARRAY" ;
    }
  
  if (strlen(oneCode) != 1 || *oneCode < 'a' || *oneCode > 'z')
    die ("auxField ONEcode %s must be a single lower-case letter", oneCode) ;
  if (strlen(spec) < 3 || spec[2] != ':' ||
      !((strlen(spec) == 4 && auxType[spec[3]]) ||
	(strlen(spec) == 5 && spec[3] == 'B' && auxType[spec[4]])))
    die ("auxField specification %s must be \"<tag>:<fmt>\" "
	 "where <tag> is a 2-letter bam tag and <fmt> is its format from [ACcSsIifHZ] or B[cCsSiIf]", spec) ;

  AuxField *a = arrayp(auxFields, arrayMax(auxFields), AuxField) ;
  a->oneCode = *oneCode ;
  spec[2] = 0 ; strcpy (a->tag, spec) ; spec[2] = ':' ;
  a->fmt = (spec[3] == 'B') ? spec[4] + 1 : spec[3] ;
}

static I32 accParse (char *acc) // truncates acc to its prefix before terminal digits
{
  char *s = acc ; while (*s) ++s ; // go to end of acc
  while (--s >= acc && *s >= '0' && *s <= '9') ; // go to last char before terminal digits (if any)
  if (s >= acc && *s == '.') *s = 0 ; // delete the version
  while (--s >= acc && *s >= '0' && *s <= '9') ; // go to last char before terminal digits (if any)
  char *term = ++s ;
  U64 n = 0 ; while (*s >= '0' && *s <= '9') n = n * 10 + (*s++ - '0') ;
  if (n > I32MAX) die ("n %lld too big in accParse", (long long) n) ;
  *term = 0 ; // truncate
  return (I32) n ;
}

static char **accTaxName ;
static I32   *accTaxN ;
static int accTaxOrder (const void *a, const void *b)
{ int retVal = strcmp(accTaxName[*(int*)a], accTaxName[*(int*)b]) ;
  if (retVal) return retVal ;
  else return accTaxN[*(int*)a] - accTaxN[*(int*)b] ;
}

static const char binaryAmbigComplement[16] =
  { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 } ;
static const char binaryAmbig2text[] = "=ACMGRSVTWYHKDBN" ;

#define ENSURE_BUF_SIZE(buf,size,newSize,Type)  if (newSize>size) { Type* newBuf = new0(newSize,Type) ; \
    memcpy(newBuf,buf,size*sizeof(Type)) ; newFree(buf,size,Type) ; size = newSize ; buf = newBuf ; }

bool bam21bam (char *bamFileName, char *outFileName, char *accTaxFileName, bool isNames)
{
  BamFile *bf = bamFileOpenRead (bamFileName) ;
  if (!bf) return false ;
  int nTargets = bf->h->n_targets ;
  printf ("read %d references: ", nTargets) ; timeUpdate (stdout) ;

  OneFile *of ;
  int   i ;
  char *s ;

  if (auxFields && arrayMax(auxFields)) // first we must extend the schema with the auxiliary linetypes
    { char *newSchemaText = new0 (strlen(schemaText) + arrayMax(auxFields)*16, char) ;
      strcpy (newSchemaText, schemaText) ;
      for (i = 0 ; i < arrayMax(auxFields) ; ++i)
	{ AuxField *a = arrp(auxFields,i,AuxField) ;
	  sprintf (newSchemaText + strlen(newSchemaText),
		   "D %c 1 %d %s\n", a->oneCode, (int)strlen(auxType[a->fmt]), auxType[a->fmt]) ;
	}
      schemaText = newSchemaText ;
    }
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  if (outFileName) of = oneFileOpenWriteNew (outFileName, schema, "bam", true, 1) ;
  else of = oneFileOpenWriteNew (derivedName (bamFileName, "1bam"), schema, "bam", true, 1) ;
  if (!of) die ("failed to open ONEfile %s to write", outFileName) ;
  oneAddProvenance (of, "onebam", VERSION, getCommandLine()) ;
  if (auxFields && arrayMax(auxFields))
    for (i = 0 ; i < arrayMax(auxFields) ; ++i)
      { AuxField *a = arrp(auxFields,i,AuxField) ;
	oneChar(of,0) = a->oneCode ;
	oneChar(of,1) = a->fmt ;
	oneWriteLine (of, 'Z', 2, a->tag) ;
      }
  
  // deal with the targets - first sort them, allowing for tid 0 = * (no target)
  int *tidMap = new(nTargets+1,int), *revMap = new(nTargets,int) ;
  accTaxName = bf->h->target_name ;
  accTaxN = new(nTargets, I32) ;
  I64 totTargetName = 0 ;
  for (i = 0; i < nTargets; ++i)
    { revMap[i] = i ;
      totTargetName += strlen (accTaxName[i]) ;
      accTaxN[i] = accParse (accTaxName[i]) ; // also truncates to stem as side effect
    }
  qsort(revMap, nTargets, sizeof(int), accTaxOrder);
  printf ("sorted %d targets total length %lld\n", nTargets, (long long)totTargetName) ;
  
  int *taxid = new0(nTargets+1,int) ;
  if (accTaxFileName) // merge the header targets to this to identify the taxids
    { OneFile *oa = oneFileOpenRead (accTaxFileName, schema, "acctax", 1) ;
      if (!oa) die ("failed to open .1acctax file %s", accTaxFileName) ;
      if (!oneGoto (oa, 'A', 1)) die ("no A lines in .1acctax file %s", accTaxFileName) ;
      oneReadLine (oa) ; // must be the first A line
      char accBuf[64] ; if (oneLen(oa) > 63) die ("acc root %s too long > 63", oneString(oa)) ;
      strcpy (accBuf, oneString(oa)) ;
      int  nFound = 0, nBlock = 0, nameComp ;
      bool isNewBlock = false ;
      for (i = 0; i < nTargets; ++i)
	{ char *acc = accTaxName[revMap[i]] ;
	  I32   k   = accTaxN[revMap[i]] ;
	  while ((nameComp = strcmp(accBuf, acc)) < 0)
	    { while (oneReadLine (oa) && !oneLen(oa)) { ; } // oneLen(oa) == 0 implies same string
	      isNewBlock = true ;
	      if (!oa->lineType) break ; // end of file
	      strcpy (accBuf + oneInt(oa,3), oneString(oa)) ; // replace suffix
	    }
	  if (!nameComp) // a match
	    { while (oneInt(oa,1) + oneInt(oa,2) <= k)
		{ if (!oneReadLine (oa)) { nameComp = 2 ; break ; } // use nameComp as flag
		  isNewBlock = true ;
		  if (oneLen(oa))
		    { strcpy (accBuf + oneInt(oa,3), oneString(oa)) ; // replace suffix
		      nameComp = 1 ;
		      break ;
		    }
		}
	      if (nameComp == 2) break ; // finish main i loop
	      if (!nameComp && oneInt(oa,1) <= k)
		{ taxid[revMap[i]+1] = oneInt(oa,0) ;
		  ++nFound ;
		  if (isNewBlock) { ++nBlock ; isNewBlock = false ; }
		}
	    }
	}
      oneFileClose (oa) ;
    }

  tidMap[0] = 0 ;
  for (i = 0 ; i < nTargets ; ++i) tidMap[revMap[i]+1] = i ;
  for (i = 0 ; i < nTargets ; ++i)
    { oneInt(of,0) = bf->h->target_len[revMap[i]] ;
      oneInt(of,1) = taxid[revMap[i]+1] ;
      char *targetName = bf->h->target_name[revMap[i]] ;
      oneWriteLine (of, 'R', strlen(targetName), targetName) ;
    }
  
  // now we process the main records
  int lastqNameSize = 128 ;   char  *lastqName = new (lastqNameSize,char) ; *lastqName = 0 ;
  int seqBufSize    = 1024 ;  char  *seq       = new (seqBufSize,char) ;
  int readGroupSize = 128 ;   char  *readGroup = new (readGroupSize,char) ; *readGroup = 0 ;
  int iBufSize      = 1024 ;  I64   *iBuf      = new (iBufSize,I64) ;
  int fBufSize      = 1024 ;  float *fBuf      = new (fBufSize,float) ;

  I64 nRecord = 0 ;
  while (true) // main loop to read sequences
    { ++nRecord ;
      int res = sam_read1 (bf->f, bf->h, bf->b) ;
      if (res < -1) die ("bamProcess failed to read bam record %lld", (long long)nRecord) ;
      if (res == -1) break ; // end of file

      // get the flag now - we will need it for various things
      I64 flag = bf->b->core.flag ; 

      char *qName = bam_get_qname(bf->b) ;
      if (strcmp (lastqName, qName)) // new query sequence !! what about if qName == * ?
	{ // first update the readGroup if necessary, then the sequence = the object
	  U8 *aux ; 
	  if ((aux = bam_aux_get (bf->b, "RG")) && (s = bam_aux2Z (aux)) && strcmp (s, readGroup))
	    { I64 len = strlen(s) ;
	      oneWriteLine (of, 'G', len, s) ; // write the new read group
	      ENSURE_BUF_SIZE (readGroup, readGroupSize, len+1, char) ;
	      strcpy (readGroup, s) ;
	    }
	  // now get the new sequence, reverse-complementing it if necesssary
	  I64 seqLen = bf->b->core.l_qseq ;
	  ENSURE_BUF_SIZE (seq, seqBufSize, seqLen+1, char) ;
	  char *bseq = (char*) bam_get_seq (bf->b) ;
	  if (flag & BAM_FREVERSE)
	    for (i = seqLen, s = seq ; i-- ; )
	      *s++ = binaryAmbig2text[(int)binaryAmbigComplement[(int)bam_seqi(bseq,i)]] ;
	  else
	    for (i = 0, s = seq ; i < seqLen ; ++i)
	      *s++ = binaryAmbig2text[bam_seqi(bseq,i)] ;
	  oneWriteLine (of, 'S', seqLen, seq) ;
	  // now we write the name change, if requested
	  ENSURE_BUF_SIZE (lastqName, lastqNameSize, bf->b->core.l_qname, char) ; // needed whether isNames or not
	  if (isNames) // write the new name suffix in a J line - saves space over I lines with full name
	    oneWriteLine (of, 'I', strlen(qName), qName) ;
	  strcpy (lastqName, qName) ;
	  for (i = 0 ; i < seqLen ; ++i) // write exceptions for non-ACGT characters
	    if (!acgtCheck[(int)seq[i]])
	      { oneInt(of,0) = i ;
		char c = oneChar(of,1) = seq[i] ;
		I64 n = 1 ; while (++i < seqLen && seq[i] == c) n++ ; oneInt(of,2) = n ;
		oneWriteLine (of, 'N', 0, 0) ;
	      }
	  // and next the qualities, also reversing them if necessary
	  char *bq = (char*) bam_get_qual (bf->b) ;	  
	  if (*bq != '\xff') // else there are no qualities
	    { if (flag & BAM_FREVERSE) // reverse qualities in place
		{ char *q = bq, *p = bq + seqLen, t ;
		  while (p-- > q) { t = *q ; *q++ = *p+33 ; *p = t+33 ; }
		}
	      else
		for (i = 0 ; i < seqLen ; ++i) bq[i] += 33 ;
	      oneWriteLine (of, 'Q', seqLen, bq) ;
	    }
	} // end of new query sequence block

      // write B line, which is the core bam line
      oneInt(of,0) = flag ;
      oneInt(of,1) = tidMap[bf->b->core.tid+1] ; // target
      oneInt(of,2) = bf->b->core.pos ;
      oneInt(of,3) = bf->b->core.qual ;
      I64 nCigar = bf->b->core.n_cigar ; // we store the cigar string as BAM U32s, but need to expand to I64
      ENSURE_BUF_SIZE (iBuf, iBufSize, nCigar, I64) ;
      U32 *cigarArray = bam_get_cigar(bf->b) ;
      for (i = 0 ; i < nCigar ; ++i) iBuf[i] = cigarArray[i] ;
      oneWriteLine (of, 'B', nCigar, iBuf) ;
      // if there is a next fragment from this template then write an X line
      if (bf->b->core.mtid >= 0)
	{ oneInt(of,0) = tidMap[bf->b->core.mtid+1] ;
	  oneInt(of,1) = bf->b->core.mpos ;
	  oneInt(of,2) = bf->b->core.isize ;
	  oneWriteLine (of, 'X', 0, 0) ;
	}

      // now do any requested aux fields
      if (auxFields)
	for (i = 0 ; i < arrayMax(auxFields) ; ++i)
	  { AuxField *a = arrp(auxFields,i,AuxField) ;
	    U8 *aux ; 
	    if ((aux = bam_aux_get (bf->b, a->tag)))
	      switch (a->fmt)
		{
		case 'A':
		  oneChar(of,0) = bam_aux2A(aux) ; oneWriteLine (of, a->oneCode, 0, 0) ; break ;
		case 'C': case 'c':
		case 'S': case 's':
		case 'I': case 'i':
		  oneInt(of,0) = bam_aux2i(aux) ; oneWriteLine (of, a->oneCode, 0, 0) ; break ;
		case 'f':
		  oneReal(of,0) = bam_aux2f(aux) ; oneWriteLine (of, a->oneCode, 0, 0) ; break ;
		case 'Z': case 'H':
		  s = bam_aux2Z(aux) ; oneWriteLine (of, a->oneCode, strlen(s), s) ; break ;
		case 'D': case 'd':
		case 'T': case 't':
		case 'J': case 'j':
		  { I64 len = bam_auxB_len(aux) ; ENSURE_BUF_SIZE (iBuf, iBufSize, len, I64) ;
		    for (i = 0 ; i < len ; ++i) iBuf[i] = bam_auxB2i(aux,i) ;
		    oneWriteLine (of, a->oneCode, len, iBuf) ; break ;
		  }
		case 'g':
		  { I64 len = bam_auxB_len(aux) ; ENSURE_BUF_SIZE (fBuf, fBufSize, len, float) ;
		    for (i = 0 ; i < len ; ++i) fBuf[i] = bam_auxB2f(aux,i) ;
		    oneWriteLine (of, a->oneCode, len, fBuf) ; break ;
		  }
		}
	  }
    }
  printf ("processed %lld records\n", (long long) nRecord) ;
  
  oneFileClose(of) ;
  newFree (tidMap, nTargets+1, int) ;
  newFree (taxid, nTargets, int) ;
  bamFileClose (bf) ;
  return true ;
}

bool onebam2bam (char *oneFileName, char *bamFileName) // reverses the previous function
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *of     = oneFileOpenRead (oneFileName, schema, "bam", 1) ;
  if (!of) die ("failed to open ONEfile %s to read", oneFileName) ;
  BamFile   *bf     = new0 (1, BamFile) ;
  if (bamFileName) bf->f = sam_open (bamFileName, "w") ;
  else bf->f = sam_open (derivedName (bamFileName, "bam"), "w") ;
  return true ;
}

/********** threaded code to read .1bam records and output .1read  ************/

typedef struct {
  OneFile   *ofIn, *ofOut ;
  I64        start, end ;
  int       *taxid ;
} B2rThread ;

typedef struct {
  U32 txid ;
  U16 count ;
  I16 bestScore ;
} TaxInfo ; // info per taxid for this thread

int txCompare (const void *a, const void *b)
{ TaxInfo *txa = (TaxInfo*)a, *txb = (TaxInfo*)b ;
  return txa->txid - txb->txid ;
}

typedef struct {
  U32 i ;
  U16 c ;
  U16 n ;
} Exception ; // must be the same size as TaxInfo, so we can abuse aTx for this

#include <ctype.h>	// for tolower(), toupper()

static void *b2rThread (void *arg)
{
  B2rThread *ti  = (B2rThread*) arg ;
  Hash       hTx = hashCreate (8192) ;
  Array      aTx = arrayCreate (2048, TaxInfo) ;

  I64   seqLenMax ; oneStats (ti->ofIn, 'S', 0, &seqLenMax, 0) ;
  
  char *mLine = new (seqLenMax, char), *mMD = new (seqLenMax, char) ;
  I32   mScore, nmCigar, nmMD ;
  I64  *mCigar = new(seqLenMax,I64) ; // seqLenMax is surely sufficient!
  bool  mFlags ;

  char  mMap[128] ;
  { char  *s = ".-acgtnACGTN", *t = ".-tgcanTGCAN" ; while (*s) mMap[*s++] = *t++ ;  }

  int i, n = 0 ;
  oneGoto (ti->ofIn, 'S', ti->start) ;
  oneReadLine (ti->ofIn) ;
  while (ti->start <= ti->end)
    { ++ti->start ;
      I64      seqLen = oneLen(ti->ofIn) ;
      char    *seq = oneDNAchar(ti->ofIn), *qual ;
      oneWriteLine (ti->ofOut, 'S', seqLen, seq) ; // write back out the sequence
      int      flags, txid, score ;
      TaxInfo *tx ;
      I64      nCigar, *i64cigar ;
      hashClear(hTx) ; arrayMax(aTx) = 0 ;
      mScore = -(1<<30) ;
      while (oneReadLine (ti->ofIn) && ti->ofIn->lineType != 'S')
	switch (ti->ofIn->lineType)
	  {
	  case 'Q': // save it for mismatch record and copy it
	    qual = oneString(ti->ofIn) ;
	    oneWriteLine (ti->ofOut, 'Q', oneLen(ti->ofIn), qual) ;
	    break ;
	  case 'I': // copy it
	    oneInt(ti->ofOut,0) = oneInt(ti->ofIn,0) ;
	    oneWriteLine (ti->ofOut,'I',oneLen(ti->ofIn),oneString(ti->ofIn)) ;
	    break ;
	  case 'B':
	    flags  = oneInt(ti->ofIn,0) ;
	    txid   = ti->taxid[oneInt(ti->ofIn,1)] ;
	    if (hashAdd (hTx, hashInt(txid), &i))
	      { tx = arrayp(aTx,i,TaxInfo) ;
		tx->txid = txid ; tx->count = 0 ; tx->bestScore = -(1<<14) ;
	      }
	    else
	      tx = arrp(aTx,i,TaxInfo) ;
	    ++tx->count ;
	    nCigar = oneLen(ti->ofIn) ; i64cigar = oneIntList(ti->ofIn) ;
	    break ;
	  case 's':
	    score = oneInt(ti->ofIn,0) ;
	    if (score > tx->bestScore) tx->bestScore = score ;
	    break ;
	  case 'm':
	    if (score > mScore) // mScore is best score for the sequence
	      { mScore = score ;
		nmCigar = nCigar ;
		memcpy (mCigar, i64cigar, nCigar*sizeof(I64)) ;
		mFlags = flags ;
		nmMD = oneLen(ti->ofIn) ;
		memcpy (mMD, oneString (ti->ofIn), nmMD) ;
	      }
	    break ;
	  default: break ;
	  }

      // code to write out the max score (mScore) and edit line (mLine)
      // we actually want the edits in the read space, not the reference space
      // first reconstruct reference string, then map that to the read with cigar
      // see  https://vincebuffalo.com/notes/2014/01/17/md-tags-in-bam-files.html
      int iS = 0 ; // position in sequence, reference
      int nM = nmMD ; char *md = mMD ; 
      int nR = 0 ; while (nM && (*md >= '0' && *md <= '9')) { nR = nR*10 + (*md++ - '0') ; --nM ; }
      int nC = nmCigar ; I64 *cigar = mCigar ;
      while (nC--)
	{ int nS = bam_cigar_oplen(*cigar), op = bam_cigar_op(*cigar++) ;
	  switch (bam_cigar_type(op))
	    {
	    case 0: break ; // do nothing
	    case 1: // INSERT (or SOFT_CLIP)
	      for (i = 0 ; i < nS ; ++i) mLine[iS++] = '-' ;
	      break ;
	    case 2: // DELETE (or REF_SKIP)
	      if (nR)
		die ("bad md match nR %d object %d", nR, (int)oneObject(ti->ofIn,'S')) ;
	      if (*md != '^')
		die ("bad delete char %c object %d", *md, (int)oneObject(ti->ofIn,'S')) ;
	      ++md ; --nM ; // eat the ^ character
	      if (nM < nS)
		die ("%d < %d delete chars object %d", nM, nS, (int)oneObject(ti->ofIn,'S')) ;
	      md += nS ; nM -= nS ; // eat the deleted chars
	      while (nM && (*md >= '0' && *md <= '9')) { nR = nR*10 + (*md++ - '0') ; --nM ; }
	      break ;
	    case 3: // MATCH (or EQUAL or DIFF)
	      while (nS)
		{ while (nS && nR) { mLine[iS++] = '.' ; --nR ; --nS ; }
		  if (nS && !nR)
		      { // encode the edit
			mLine[iS] = (qual[iS]<25)?tolower(*md++):toupper(*md++) ; --nM ;
			++iS ; --nS ;
			while (nM && (*md >= '0' && *md <= '9')) {nR = nR*10 +(*md++ - '0'); --nM; }
		      }
		}
	      break ;
	    }
	}
      if (iS != seqLen) die ("iS %d != seqLen %d", iS, seqLen) ;
      if (mFlags & BAM_FREVERSE) // need to reverse-complement mLine
	{ int i, j ;
	  for (i = 0, j = seqLen ; i < --j ; ++i)
	    { char t = mLine[i] ; mLine[i] = mMap[mLine[j]] ; mLine[j] = mMap[t] ; }
	  if (i == j) mLine[i] = mMap[mLine[i]] ;
	}
      oneInt(ti->ofOut,0) = mScore ;
      oneWriteLine (ti->ofOut, 'M', seqLen, mLine) ;
      
      if (arrayMax(aTx) > 1) arraySort (aTx, txCompare) ;
      for (i = 0 ; i < arrayMax(aTx) ; ++i)
	{ tx = arrp(aTx,i,TaxInfo) ;
	  oneInt(ti->ofOut,0) = tx->txid ;
	  oneInt(ti->ofOut,1) = tx->bestScore ;
	  oneInt(ti->ofOut,2) = tx->count ;
	  oneWriteLine (ti->ofOut, 'T', 0, 0) ;
	}
    }

  newFree (mLine, seqLenMax, char) ;
  newFree (mMD, seqLenMax, char) ;
  newFree (mCigar, seqLenMax, I64) ;
  hashDestroy (hTx) ; arrayDestroy (aTx) ;
  return 0 ;
}

// The next function makes a 1read file from a 1bam file.  This can be threaded.

bool bamMake1read (char *inFileName, char *outFileName)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *ofIn = oneFileOpenRead (inFileName, 0, "bam", NTHREAD) ;
  if (!ofIn) { warn ("failed to open .1bam file %s", inFileName) ; return false ; }
  if (!oneFileCheckSchemaText (ofIn, "P 3 seq\nD s 1 3 INT\nD m 1 6 STRING\n"))
    die ("input .1bam file must have s (score) and m (MD) record types") ;
  if (!outFileName) outFileName = derivedName (inFileName, "1read") ;
  OneFile *ofOut = oneFileOpenWriteNew (outFileName, schema, "read", true, NTHREAD) ;
  if (!ofOut) { oneFileClose(ofIn) ; warn ("failed to open %s", outFileName) ; return false ; }
  oneAddProvenance (ofOut, "onebam", VERSION, getCommandLine()) ;

  I64 nRef ; oneStats (ofIn, 'R', &nRef, 0, 0) ;
  int *taxid = new (nRef+1, int) ;
  int i = 0 ;
  taxid[i++] = 0 ; // initial taxid is 0
  oneGoto (ofIn, 'R', 1) ;
  while (oneReadLine(ofIn) && ofIn->lineType == 'R') taxid[i++] = oneInt(ofIn,1) ;
  fprintf (stderr, "read %d references\n", (int) nRef) ;
  
  pthread_t *threads = new (NTHREAD, pthread_t) ;
  B2rThread *ti      = new0 (NTHREAD, B2rThread) ;
  I64 nSeq ; oneStats (ofIn, 'S', &nSeq, 0, 0) ;
  for (i = 0 ; i < NTHREAD ; ++i)
    { ti[i].ofIn  = ofIn + i ;
      ti[i].ofOut = ofOut + i ;
      ti[i].start = 1 + (nSeq * i) / NTHREAD ;
      ti[i].end = (nSeq * (i+1)) / NTHREAD ;
      ti[i].taxid = taxid ;
    }
  for (i = 0 ; i < NTHREAD ; ++i) // create threads
    pthread_create (&threads[i], 0, b2rThread, &ti[i]) ;
  for (i = 0 ; i < NTHREAD ; ++i)
    pthread_join (threads[i], 0) ; // wait for threads to complete
  newFree (threads, NTHREAD, pthread_t) ;
  newFree (ti, NTHREAD, B2rThread) ;
  
  oneFileClose (ofIn) ;
  oneFileClose (ofOut) ;
  newFree (taxid, nRef, int) ;
  fprintf (stderr, "built read records for %lld sequences: ", (long long) nSeq) ;
  timeUpdate (stderr) ;
  return true ;
}

/*******************************************************************************************/
/*********** next direct bam to 1read ******************************************************/
/*********** this is based on a merge of bam21bam() and bamMake1read() by Claude ***********/

// this routine calculates the mismatch line in read space from the cigar and MD string
// see  https://vincebuffalo.com/notes/2014/01/17/md-tags-in-bam-files.html for MD interpretation

static void generateMline (char *mLine, int seqLen, I32 *cigar, char *md, bool isReverse)
{
  static char mMap[128] = "" ;
  static bool isFirst = true ;
  if (isFirst)
    { char *s = ".-acgtnACGTN", *t = ".-tgcanTGCAN";
      while (*s) mMap[*s++] = *t++;
      isFirst = false ;
    }

  int iS = 0 ; // position in sequence
  int nCigar = 0 ; while (cigar[nCigar]) ++nCigar ;
  int nM = strlen(md) ;
  int nR = 0 ; while (nM && (*md >= '0' && *md <= '9')) { nR = nR*10 + (*md++ - '0') ; --nM ; }
  while (nCigar--)
    { int nS = bam_cigar_oplen(*cigar), op = bam_cigar_op(*cigar++) ;
      switch (bam_cigar_type(op))
	{
	case 0: // do nothing
	  break ;
	case 1: // INSERT (or SOFT_CLIP)
	  { int i ; for (i = 0 ; i < nS ; ++i) mLine[iS++] = '-' ; }
	  break ;
	case 2: // DELETE (or REF_SKIP)
	  if (nR)         die ("bad md match nR %d", nR) ;
	  if (*md != '^') die ("bad delete char %c", *md) ;
	  ++md ; --nM ; // eat the ^ character
	  if (nM < nS)    die ("%d < %d delete chars", nM, nS) ;
	  md += nS ; nM -= nS ; // eat the deleted chars
	  while (nM && (*md >= '0' && *md <= '9')) { nR = nR*10 + (*md++ - '0') ; --nM ; }
	  break ;
	case 3: // MATCH (or EQUAL or DIFF)
	  while (nS)
	    { while (nS && nR) { mLine[iS++] = '.' ; --nR ; --nS ; }
	      if (nS && !nR)
		{ // encode the edit
		  mLine[iS] = *md++ ; --nM ;
		  ++iS ; --nS ;
		  while (nM && (*md >= '0' && *md <= '9')) {nR = nR*10 +(*md++ - '0'); --nM; }
		}
	    }
	  break ;
	}
    }
  if (iS != seqLen) die ("generateMline: iS %d != seqLen %d", iS, seqLen) ;
  
  if (isReverse) // need to reverse-complement mLine
    { int i, j ;
      for (i = 0, j = seqLen ; i < --j ; ++i)
	{ char t = mLine[i] ; mLine[i] = mMap[mLine[j]] ; mLine[j] = mMap[t] ; }
      if (i == j) mLine[i] = mMap[mLine[i]] ;
    }
}

#include <unistd.h>   // for getpid()

static inline U8* bufPushChar (U8* b, char x) { char *bb = (char*)b ; *bb++ = x ; return (U8*)bb ; }
static inline U8*  bufPushI32 (U8* b, I32 x) { I32 *bb = (I32*)b ; *bb++ = x ; return (U8*)bb ; }
static inline char bufPopChar (U8* *b) { char x = *(char*)*b ; *b += sizeof(char) ; return x ; }
static inline I32   bufPopI32 (U8* *b) { I32 x = *(I32*)*b ; *b += sizeof(I32) ; return x ; }

// write the current contents of buf into a OneFile
static void write1read (Array oneFileNames, Array bLoc,
			  char *firstName, int minNameShare, int maxNameLen)
{
  static OneSchema *schema = 0 ;
  static char nameBase[64], *nameTail ;
  if (!schema) // initialise
    { schema = oneSchemaCreateFromText (schemaText) ;
      sprintf (nameBase, "./onebam-%d-", getpid()) ;
      nameTail = nameBase + strlen(nameBase) ;
    }

  if (!arrayMax(bLoc)) die ("write1read() called empty - no BAM records or buffer too small") ;
  printf ("entering write1read %d for %d items: ",
	  (int)arrayMax(oneFileNames), (int)arrayMax(bLoc)) ;
  timeUpdate (stdout) ;
  { char c = firstName[minNameShare] ; firstName[minNameShare] = 0 ;
    printf ("  name prefix %s maxNameLen %d\n", firstName, maxNameLen) ;
    firstName[minNameShare] = c ;
  }

  // open the OneFile and store its name in names
  sprintf (nameTail, "%d", (int)arrayMax(oneFileNames)) ;
  OneFile *of = oneFileOpenWriteNew (nameBase, schema, "read", true, 1) ;
  if (!of) die ("failed to open ONEfile %s to write", nameBase) ;
  array(oneFileNames,arrayMax(oneFileNames),char*) = strdup (nameBase) ;
  oneWriteLine (of, 'P', minNameShare, firstName) ;

  // now set up an array to sort the query names and sort using Gene's msd_sort()
  int ksize = maxNameLen + 1 - minNameShare ; // add the +1 so 0-terminated
  int rsize = ksize + sizeof(I32) ;
  U8 *a = new0(rsize*arrayMax(bLoc), U8) ; // need new0 to give terminating NULs on names
  I32 i, n = arrayMax(bLoc) ;
  for (i = 0 ; i < n ; ++i)
    { U8 *b = arr(bLoc,i,U8*) ;
      if (*b++ != 'S') die ("no S found at start of sequence %d", i) ;
      U8 nameLen = *b++, nameShare = *b++ ;
      char *s = (char*)(a + rsize*i) ;
      int shared = nameShare - minNameShare ;
      if (shared) { memcpy (s, firstName+minNameShare, shared) ; s += shared ; }
      memcpy (s, (char*)b, nameLen - nameShare) ;
      *(I32*)(a + rsize*i + ksize) = i ;
    }
  msd_sort (a, n, rsize, ksize, 0, 0, 1) ;

  // now write them out, in the new sorted order
  // we may have to merge multiple entries for the same read
  char  lastName[256] = "\n" ;
  I32   maxScore = - (1<<30) ;
  char *mLine = 0 ;
  Array aTx = arrayCreate (1024, TaxInfo) ;
  I32   seqLen ;
  for (i = 0 ; i < n ; ++i)
    { char *nameEnd = (char*)(a + rsize*i) ;
      I32 j = *(I32*)(nameEnd + ksize) ;
      U8 *b = arr(bLoc,j,U8*) ; ++b ; // swallow the 'S' - checked above
      U8  nameLen = *b++, nameShare = *b++ ; b += nameLen - nameShare ; // and the name stuff
      if (strcmp (lastName, nameEnd)) // new read
	{ if (arrayMax (aTx)) // write out max score and taxid lines
	    { oneInt(of,0) = maxScore ;
	      oneWriteLine (of, 'M', seqLen, mLine) ;
	      if (arrayMax(aTx) > 1) // sort and pack aTx
		{ arraySort (aTx, txCompare) ;
		  TaxInfo *iTx = arrp(aTx,-1,TaxInfo), *jTx = arrp(aTx,0,TaxInfo) ;
		  TaxInfo *nTx = jTx + arrayMax(aTx) ; // end of array
		  while (jTx < nTx)
		    { if (++iTx < jTx) *iTx = *jTx ;
		      while (++jTx < nTx && jTx->txid == iTx->txid)
			{ iTx->count += jTx->count ;
			  if (jTx->bestScore > iTx->bestScore) iTx->bestScore = jTx->bestScore ;
			}
		    }
		  arrayMax(aTx) = iTx - arrp(aTx,0,TaxInfo) + 1 ;
		}
	      int  it ;
	      TaxInfo *tx =  arrp(aTx,0,TaxInfo) ;
	      for (it = 0 ; it < arrayMax(aTx) ; ++it, ++tx)
		{ oneInt(of,0) = tx->txid ; oneInt(of,1) = tx->bestScore ; oneInt(of,2) = tx->count ;
		  oneWriteLine (of, 'T', 0, 0) ;
		}
	      arrayMax(aTx) = 0 ; // clear the array
	      maxScore = - (1<<30) ; // and reset maxScore
	    }
	  seqLen = bufPopI32 (&b) ;
	  oneWriteLine (of, 'S', seqLen, b) ; b += seqLen ;
	  oneWriteLine (of, 'I', nameLen-minNameShare, nameEnd) ;
	  strncpy (lastName, nameEnd, ksize) ;
	  I32 nException = bufPopI32 (&b) ;
	  for (j = 0 ; j < nException ; ++j)
	    { Exception e = *(Exception*)b ; b += sizeof(Exception) ;
	      oneInt(of,0) = e.i ; oneChar(of,1) = e.c ; oneInt(of,2) = e.n ;
	      oneWriteLine (of, 'N', 0, 0) ;
	    }
	  if (*b == 'Q') // qualities
	    { oneWriteLine (of, 'Q', seqLen, ++b) ;
	      b += seqLen ;
	    }
	}
      else // have to bypass exceptions and qualities - don't check consistency
	{ seqLen = bufPopI32 (&b) ;
	  b += seqLen ;
	  I32 nException = bufPopI32 (&b) ;
	  b += nException * sizeof(Exception) ;
	  if (*b == 'Q') b += 1 + seqLen ;
	}
      if (*b == 'T') // add any taxon information onto aTx
	{ ++b ;
	  I32 mScore = bufPopI32 (&b) ;
	  if (mScore > maxScore) { maxScore = mScore ; mLine = (char*) b ; }
	  b += seqLen ;
	  int nt =  bufPopI32 (&b) ; // number of taxids
	  int oldMax = arrayMax(aTx) ;
	  arrayp(aTx,oldMax+nt-1,TaxInfo)->txid = 1 ; // set arrayMax and make space
	  memcpy (arrp(aTx,oldMax,TaxInfo), (TaxInfo*)b, nt*sizeof(TaxInfo)) ;
	}
    }

  if (arrayMax (aTx)) // need to write the final block
    { oneInt(of,0) = maxScore ;
      oneWriteLine (of, 'M', seqLen, mLine) ;
      if (arrayMax(aTx) > 1) // sort and pack aTx
	{ arraySort (aTx, txCompare) ;
	  TaxInfo *iTx = arrp(aTx,-1,TaxInfo), *jTx = arrp(aTx,0,TaxInfo) ;
	  TaxInfo *nTx = jTx + arrayMax(aTx) ; // end of array
	  while (jTx < nTx)
	    { if (++iTx < jTx) *iTx = *jTx ;
	      while (++jTx < nTx && jTx->txid == iTx->txid)
		{ iTx->count += jTx->count ;
		  if (jTx->bestScore > iTx->bestScore) iTx->bestScore = jTx->bestScore ;
		}
	    }
	  arrayMax(aTx) = iTx - arrp(aTx,0,TaxInfo) + 1 ;
	}
      int  it ;
      TaxInfo *tx =  arrp(aTx,0,TaxInfo) ;
      for (it = 0 ; it < arrayMax(aTx) ; ++it, ++tx)
	{ oneInt(of,0) = tx->txid ; oneInt(of,1) = tx->bestScore ; oneInt(of,2) = tx->count ;
	  oneWriteLine (of, 'T', 0, 0) ;
	}
    }
  
  oneFileClose (of) ;
  newFree (a, rsize*arrayMax(bLoc), U8) ;

  printf ("  done write1read: ") ; timeUpdate (stdout) ;
}

bool bam21readSorted (char *bamFileName, char *outFileName, char *accTaxFileName)
{
  if (!outFileName) outFileName = derivedName(bamFileName, "1read") ;
  
  BamFile *bf = bamFileOpenRead(bamFileName);
  if (!bf) return false;
  int nTargets = bf->h->n_targets;
  printf("read %d references: ", nTargets); timeUpdate(stdout);

  if (sizeof(Exception) != sizeof(TaxInfo))
    die ("mismatch sizeof(Exception) %d != sizeof(TaxInfo) %d",
	 sizeof(Exception), sizeof(TaxInfo)) ;
  
  // deal with targets - first sort them
  int i, *revMap = new(nTargets, int);
  accTaxName = bf->h->target_name ;
  accTaxN = new(nTargets, I32) ;
  I64 totTargetName = 0 ;
  for (i = 0; i < nTargets; ++i)
    { revMap[i] = i ;
      totTargetName += strlen (accTaxName[i]) ;
      accTaxN[i] = accParse (accTaxName[i]) ; // also truncates to stem as side effect
    }
  qsort(revMap, nTargets, sizeof(int), accTaxOrder);
  printf ("sorted %d targets total length %lld\n", nTargets, (long long)totTargetName) ;
 
  // and next map to taxids by merging with acc2taxid information in .1acctax file
  int  *taxid = new0(nTargets+1, int) ;
  OneSchema *schema = oneSchemaCreateFromText(schemaText);
  OneFile *oa = oneFileOpenRead (accTaxFileName, schema, "acctax", 1) ;
  if (!oa) die ("failed to open .1acctax file %s", accTaxFileName) ;
  if (!oneGoto (oa, 'A', 1)) die ("no A lines in .1acctax file %s", accTaxFileName) ;
  oneReadLine (oa) ; // must be the first A line
  char accBuf[64] ; if (oneLen(oa) > 63) die ("acc root %s too long > 63", oneString(oa)) ;
  strcpy (accBuf, oneString(oa)) ;
  int  nFound = 0, nBlock = 0, nameComp ;
  bool isNewBlock = false ;
  for (i = 0; i < nTargets; ++i)
    { char *acc = accTaxName[revMap[i]] ;
      I32   k   = accTaxN[revMap[i]] ;
      while ((nameComp = strcmp(accBuf, acc)) < 0)
	{ while (oneReadLine (oa) && !oneLen(oa)) { ; } // oneLen(oa) == 0 implies same string
	  isNewBlock = true ;
	  if (!oa->lineType) break ; // end of file
	  strcpy (accBuf + oneInt(oa,3), oneString(oa)) ; // replace suffix
	}
      if (!nameComp) // a match
	{ while (oneInt(oa,1) + oneInt(oa,2) <= k)
	    { if (!oneReadLine (oa)) { nameComp = 2 ; break ; } // use nameComp as flag
	      isNewBlock = true ;
	      if (oneLen(oa))
		{ strcpy (accBuf + oneInt(oa,3), oneString(oa)) ; // replace suffix
		  nameComp = 1 ;
		  break ;
		}
	    }
	  if (nameComp == 2) break ; // finish main i loop
	  if (!nameComp && oneInt(oa,1) <= k)
	    { taxid[revMap[i]+1] = oneInt(oa,0) ;
	      ++nFound ;
	      if (isNewBlock) { ++nBlock ; isNewBlock = false ; }
	    }
	  //	  else
	  //	    printf ("failed to match %s %d to block %s %d-%d\n",
	  //		    acc, k, accBuf, (int)oneInt(oa,1), (int)oneInt(oa,2)) ;
	}
    }
  oneFileClose (oa) ;
  newFree(revMap, nTargets, int) ;
  newFree(accTaxN, nTargets, I32) ;
  printf("found %d taxids in %d acctax blocks: ", nFound, nBlock); timeUpdate(stdout) ;

  // strategy is to write to a series of .1read files (2GB in memory) then merge
  Array ofNames  = arrayCreate (16, char*) ;

  // buffer in which we will accumulate the objects as we process them, then sort when full
  I64   maxBuf   = (I64)1 << 31 ;
  U8   *bufStart = new (maxBuf, U8) ;
  U8   *bufEnd   = bufStart ;
  Array bLoc     = arrayCreate (maxBuf>>7, U8*) ; // start of each object

  // buffers for processing
  char lastqName[256] ; *lastqName = 0 ;
  char firstName[256] ; *firstName = 0 ;
  int  minNameShare = 256 ;
  int  maxNameLen   = 0 ;
  int  seqBufSize   = 1024; char  *seq   = new(seqBufSize, char);
  int  mLineSize    = 1024; char  *mLine = new(mLineSize, char);
  // int readGroupSize = 128;  char *readGroup = new(readGroupSize, char); *readGroup = 0;
  
  // hash table and array for taxonomic info per sequence
  Hash  hTx = hashCreate(8192);
  Array aTx = arrayCreate(2048, TaxInfo);

  // other counters etc.
  I64  nRecord      = 0;
  I64  seqLen       = 0;
  I32  maxScore     = -(1<<30);
  bool isMaxReverse = false;
  int  maxCigarSize = 1024; I32  *maxCigar = new(maxCigarSize, I32);
  int  maxMdSize    = 1024; char *maxMd    = new(maxMdSize, char);

  while (true) // main loop to read sequences
    { ++nRecord;
      int res = sam_read1 (bf->f, bf->h, bf->b) ;
      if (res < -1) die ("bamProcess failed to read bam record %lld", (long long)nRecord) ;
      if (res == -1) break ; // end of file

      I64   flag  = bf->b->core.flag ;
      char *qName = bam_get_qname(bf->b) ;
    
      if (strcmp(lastqName, qName)) // new query sequence
	{ int nameLen = strlen (qName) ;
	  if (nameLen > maxNameLen)
	    { maxNameLen = nameLen ;
	      if (nameLen > 255) die ("query name %s longer than 255", qName) ;
	    }

	  // process previous sequence if we have one
	  if (!*firstName)
	    strcpy (firstName, qName) ;
	  else if (maxScore > -(1<<30))
	    {
	      // if will overflow buffer, write out oneFile and restart buffer
	      if (bufEnd - bufStart
		  + 1 + sizeof(I32) + seqLen + // 'M' + maxScore + mLine
		  (1 + 3*arrayMax(aTx))*sizeof(I32) // tax info
		  > maxBuf)
		{ --arrayMax(bLoc) ;
		  U8 *bufLast = arr(bLoc, arrayMax(bLoc), U8*) ;
		  write1read (ofNames, bLoc, firstName, minNameShare, maxNameLen) ;
		  memcpy (bufStart, bufLast, bufEnd - bufLast) ;
		  arrayMax(bLoc) = 1 ; bufEnd -= (bufLast-bufStart) ;
		}
	      
	      // store max score and mismatch line from stored best alignment
	      bufEnd = bufPushChar (bufEnd, 'T') ;
	      bufEnd = bufPushI32 (bufEnd, maxScore) ;
	      ENSURE_BUF_SIZE(mLine, mLineSize, seqLen+1, char) ;
	      generateMline (mLine, seqLen, maxCigar, maxMd, isMaxReverse) ;
	      memcpy (bufEnd, mLine, seqLen) ; bufEnd += seqLen ;
	      
	      // store taxonomic information
	      bufEnd = bufPushI32 (bufEnd, (I32)arrayMax(aTx)) ;
	      memcpy (bufEnd, arrp(aTx,0,TaxInfo), arrayMax(aTx)*sizeof(TaxInfo)) ;
	      bufEnd += arrayMax(aTx)*sizeof(TaxInfo) ;
	    }
      
	  // Reset for new sequence
	  hashClear(hTx) ;
	  arrayMax(aTx) = 0 ;
	  maxScore      = -(1<<30) ;
      
	  /* Ignore read groups for now
	  // The correct way to deal with them would be with a DICT, written with counts at the end
	  // Then code reading .1read files would have to read in the DICT first
	  U8 *aux;
	  if ((aux = bam_aux_get(bf->b, "RG")) && (s = bam_aux2Z(aux)) && strcmp(s, readGroup))
	    { I64 len = strlen(s);
	      oneWriteLine(of, 'G', len, s);
	      ENSURE_BUF_SIZE(readGroup, readGroupSize, len+1, char);
	      strcpy(readGroup, s);
	    }
	  */

	  // Get new sequence
	  seqLen = bf->b->core.l_qseq;
	  ENSURE_BUF_SIZE(seq, seqBufSize, seqLen+1, char);
	  char *s = seq, *bs = (char*)bam_get_seq(bf->b);
	  if (flag & BAM_FREVERSE)
	    for (i = seqLen; i--; )
	      *s++ = binaryAmbig2text[(int)binaryAmbigComplement[(int)bam_seqi(bs,i)]];
	  else
	    for (i = 0; i < seqLen; ++i)
	      *s++ = binaryAmbig2text[bam_seqi(bs,i)];
	  *s = 0; // null terminate

	  // calculate how much of qName to store
	  char *p = firstName, *q = qName;
	  while (*p && *p == *q) { ++p; ++q; }
	  int nameShare = (q - qName) ;
	  if (nameShare < minNameShare) minNameShare = nameShare ;
      
	  // store non-ACGT sequence exceptions, borrowing aTx since it is available and right-sized
	  for (i = 0; i < seqLen; ++i)
	    if (!acgtCheck[(int)seq[i]])
	      { Exception *e = arrayp (aTx, arrayMax(aTx), Exception) ;
		e->i = i ;
		e->c = seq[i] ;
		I64 n = 1 ;
		while (++i < seqLen && seq[i] == e->c) n++ ;
		e->n = n ;
	      }
      
	  // get qualities, reversing them if necesssary
	  q = (char*)bam_get_qual(bf->b); // reuse p, q as utility char* from J line code
	  if (*q != '\xff')
	    { if (flag & BAM_FREVERSE)
		{ char t ;
		  for (p = q + seqLen ; p-- > q ; ++q) { t = *p + 33 ; *p = *q + 33 ; *q = t ; }
		}
	      else
		for (i = 0; i < seqLen; ++i) q[i] += 33;
	    }

	  if (bufEnd - bufStart
	      + 1                                   // new sequence marker 'S'
	      + 2 + nameLen - nameShare             // name
	      + sizeof(I32) + seqLen                // seqLen + seq
	      + sizeof(I32) + arrayMax(aTx)*sizeof(Exception) // exceptions
	      + ((*q != '\xff')?(1+seqLen):0)         // qualities
	      > maxBuf)
	    { write1read (ofNames, bLoc, firstName, minNameShare, maxNameLen) ;
	      arrayMax(bLoc) = 0 ; bufEnd = bufStart ; 
	    }

	  // now write the location of the start of this object into the bLoc array
	  array (bLoc, arrayMax(bLoc), U8*) = bufEnd ;

	  // now write all the required information into the buffer
	  bufEnd = bufPushChar (bufEnd, 'S') ;
	  *bufEnd++ = nameLen ;
	  *bufEnd++ = nameShare ;
	  memcpy (bufEnd, qName + nameShare, nameLen-nameShare) ; bufEnd += nameLen - nameShare ;
	  bufEnd = bufPushI32 (bufEnd, seqLen) ;
	  memcpy (bufEnd, seq, seqLen) ; bufEnd += seqLen ;
	  bufEnd = bufPushI32 (bufEnd, (I32)arrayMax(aTx)) ;
	  memcpy (bufEnd, arrp(aTx,0,Exception), arrayMax(aTx)*sizeof(Exception)) ;
	  bufEnd += arrayMax(aTx)*sizeof(Exception) ;
	  arrayMax(aTx) = 0 ; // need to reset this so it can accumuate taxon information
	  if (*q == '\xff')
	    { bufEnd = bufPushChar (bufEnd, 'Q') ;
	      memcpy (bufEnd, q, seqLen) ; bufEnd += seqLen ;
	    }

	  strcpy(lastqName, qName); // first don't forget to update lastqName
	} // end of new sequence block

      // Process current alignment - first find the txid
      TaxInfo *tx ;
      int j, txid = taxid[bf->b->core.tid+1] ;
      if (hashAdd(hTx, hashInt(txid), &j))
	{ tx = arrayp(aTx, j, TaxInfo);
	  tx->txid = txid;
	  tx->count = 0;
	  tx->bestScore = -(1<<14);
	}
      else 
	tx = arrp(aTx, j, TaxInfo);
      ++tx->count;
    
      // Get alignment score
      U8 *aux;
      if (!(aux = bam_aux_get(bf->b, "AS"))) continue; // skip if no score
      int score = bam_aux2i(aux);
      if (score > tx->bestScore)
	{ tx->bestScore = score;
	  // Check if this is the best alignment for mLine calculation
	  if (score > maxScore)
	    { maxScore = score;
	      isMaxReverse = (flag & BAM_FREVERSE) ? true : false;
	      // Store CIGAR
	      I64 nCigar = bf->b->core.n_cigar;
	      ENSURE_BUF_SIZE(maxCigar, maxCigarSize, nCigar+1, I32);
	      memcpy (maxCigar, bam_get_cigar(bf->b), nCigar*sizeof(I32)) ;
	      maxCigar[nCigar] = 0 ; // use this in generateMline()

	      // Store MD string
	      char *s ;
	      if ((aux = bam_aux_get(bf->b, "MD")) && (s = bam_aux2Z(aux)))
		{ ENSURE_BUF_SIZE(maxMd, maxMdSize, strlen(s)+1, char);
		  strcpy(maxMd, s);
		}
	      else
		*maxMd = 0 ;
	    }
	}
    }
    
  // Process final sequence
  if (*firstName)
    { if (maxScore > -(1<<30)) // same code as above
	{ if (bufEnd - bufStart
	      + 1 + sizeof(I32) + seqLen + // 'M' + maxScore + mLine
	      (1 + 3*arrayMax(aTx))*sizeof(I32) // tax info
	      > maxBuf)
	    { --arrayMax(bLoc) ;
	      U8 *bufLast = arr(bLoc, arrayMax(bLoc), U8*) ;
	      write1read (ofNames, bLoc, firstName, minNameShare, maxNameLen) ;
	      memcpy (bufStart, bufLast, bufEnd - bufLast) ;
	      arrayMax(bLoc) = 1 ; bufEnd -= (bufLast-bufStart) ;
	    }

	  bufEnd = bufPushChar (bufEnd, 'T') ;
	  bufEnd = bufPushI32 (bufEnd, maxScore) ;
	  ENSURE_BUF_SIZE(mLine, mLineSize, seqLen+1, char);
	  generateMline (mLine, seqLen, maxCigar, maxMd, isMaxReverse);
	  memcpy(bufEnd,mLine,seqLen) ; bufEnd += seqLen ;
	      
	  // store taxonomic information
	  // if (arrayMax(aTx) > 1) arraySort(aTx, txCompare) ; // don't need, will sort later
	  bufEnd = bufPushI32 (bufEnd, (I32)arrayMax(aTx)) ;
	  memcpy (bufEnd, arrp(aTx,0,TaxInfo), arrayMax(aTx)*sizeof(TaxInfo)) ;
	  bufEnd += arrayMax(aTx)*sizeof(TaxInfo) ;
	}
      write1read (ofNames, bLoc, firstName, minNameShare, maxNameLen) ;
    }

  printf ("processed %lld records\n", (long long)nRecord) ;

  if (arrayMax(ofNames) == 1) // move it
    rename (arr(ofNames,0,char*), outFileName) ;
  else if (merge1read (outFileName, arrayMax(ofNames), arrp(ofNames,0,char*)))
    { int i ;
      for (i = 0 ; i < arrayMax(ofNames) ; ++i) remove (arr(ofNames,i,char*)) ;
    }
  else
    die ("failed to merge %d files", (int)arrayMax(ofNames)) ;
  
  // Cleanup
  newFree(taxid, nTargets+1, int);
  newFree(seq, seqBufSize, char);
  newFree(mLine, mLineSize, char);
  newFree(maxCigar, maxCigarSize, I32);
  newFree(maxMd, maxMdSize, char);
  hashDestroy(hTx);
  arrayDestroy(aTx);
  { int i ; for (i = 0 ; i < arrayMax(ofNames) ; ++i) free(arr(ofNames,i,char*)) ; }
  arrayDestroy (ofNames) ;
  arrayDestroy (bLoc) ;
  bamFileClose(bf);
  
  return true;
}

/*********************** end of file ***********************/

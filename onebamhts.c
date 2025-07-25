/*  File: onebamhts.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul 17 10:53 2025 (rd109)
 * Created: Wed Jul  2 13:39:53 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"
#include "sam.h"
#include "array.h"

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

static int accOrder (void *names, const void *a, const void *b)
{ return strcmp(((char**)names)[*(int*)a], ((char**)names)[*(int*)b]) ; }
#if (!defined(_GNU_SOURCE) && FALSE)
static void *qsort_arg ;
static int accOrderBare  (const void *a, const void *b)
{ return strcmp(((char**)qsort_arg)[*(int*)a], ((char**)qsort_arg)[*(int*)b]) ; }
static void qsort_r(void *base, size_t n, size_t size, void *arg, int (*comp)(const void *, const void *, void *))
{ qsort_arg = arg ; qsort (base, n, size, accOrderBare) ; }
#endif

static const char binaryAmbigComplement[16] =
  { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 } ;
static const char binaryAmbig2text[] = "=ACMGRSVTWYHKDBN" ;
static const int  acgtCheck[] = {
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  1,  0,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,
  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
} ;

#define ENSURE_BUF_SIZE(buf,size,newSize,Type)  if (newSize>size) { Type* newBuf = new0(newSize,Type) ; \
    memcpy(newBuf,buf,size*sizeof(Type)) ; newFree(buf,size,Type) ; size = newSize ; buf = newBuf ; }

bool bam21bam (char *bamFileName, char *outFileName, char *taxidFileName, bool isNames)
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
  if (!outFileName)
    { for (s = bamFileName + strlen(bamFileName) ; *s != '.' && s > bamFileName ; --s) ;
      if (s > bamFileName)
	{ outFileName = new0 ((s-bamFileName)+6,char) ; strncpy (outFileName, bamFileName, (s-bamFileName)) ; }
      else
	{ outFileName = new0 (strlen(bamFileName)+6,char) ; strcpy (outFileName, bamFileName) ; }
      strcat (outFileName, ".1bam") ;
    }
  of = oneFileOpenWriteNew (outFileName, schema, "bam", true, 1) ;
  if (!of) die ("failed to open ONEfile %s to write", outFileName) ;
  oneAddProvenance (of, "onebam", "0.1", getCommandLine()) ;
  if (auxFields && arrayMax(auxFields))
    for (i = 0 ; i < arrayMax(auxFields) ; ++i)
      { AuxField *a = arrp(auxFields,i,AuxField) ;
	oneChar(of,0) = a->oneCode ;
	oneChar(of,1) = a->fmt ;
	oneWriteLine (of, 'Z', 2, a->tag) ;
      }
  
  // deal with the targets - first sort them, allowing for tid 0 = * (no target)
  int *tidMap = new(nTargets+1,int), *revMap = new(nTargets,int) ;
  for (i = 0 ; i < nTargets ; ++i) revMap[i] = i ;
  qsort_r (revMap, nTargets, sizeof(int), bf->h->target_name, accOrder) ;
  int *taxid = new0(nTargets+1,int) ;
  if (taxidFileName) // merge the header targets to this to identify the taxids
    { FILE *tf = fopen (taxidFileName, "r") ; if (!tf) die ("failed to open taxid file %s", taxidFileName) ;
      char accBuf[32] ; *accBuf = 0 ;
      int  tid, nFound = 0, comp ;
      for (i = 0 ; i < nTargets ; ++i)
	{ char *acc = bf->h->target_name[revMap[i]] ;
	  while ((comp = strcmp (accBuf, acc)) < 0) if (fscanf (tf, "%s\t%d\n", accBuf, &tid) != 2) break ;
	  if (!comp) { taxid[revMap[i]+1] = tid ; ++nFound ; }
	}
      fclose (tf) ;
      printf ("found %d taxids: ", nFound) ; timeUpdate (stdout) ;
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
	    { newFree (seq, seqBufSize, char) ; seqBufSize = seqLen * 1.5 ; seq = new(seqBufSize,char) ; }
	  char *bseq = (char*) bam_get_seq (bf->b) ;
	  if (flag & BAM_FREVERSE)
	    for (i = seqLen, s = seq ; i-- ; )
	      *s++ = binaryAmbig2text[(int)binaryAmbigComplement[(int)bam_seqi(bseq,i)]] ;
	  else
	    for (i = 0, s = seq ; i < seqLen ; ++i)
	      *s++ = binaryAmbig2text[bam_seqi(bseq,i)] ;
	  oneWriteLine (of, 'S', seqLen, seq) ;
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
	  // now we can write the name change, if requested
	  ENSURE_BUF_SIZE (lastqName, lastqNameSize, bf->b->core.l_qname, char) ; // needed whether isNames or not
	  if (isNames) // write the new name suffix in a J line - saves space over I lines with full name
	    { char *p = lastqName, *q = qName ;
	      while (*p && *p == *q) { ++p ; ++q ; }
	      oneInt(of,0) =(I64)(q - qName) ;
	      oneWriteLine (of, 'J', strlen(q), q) ;
	    }
	  strcpy (lastqName, qName) ;
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

/********** threaded code to read in the reference ************/

typedef struct {
  OneFile   *ofIn, *ofOut ;
  I64        start, end ;
  int       *taxid ;
} B2rThread ;

typedef struct {
  int txid ;
  int count ;
  int bestScore ;
} TaxInfo ; // info per taxid for this thread

#include <ctype.h>	// for tolower(), toupper()

static void *b2rThread (void *arg)
{
  B2rThread *ti = (B2rThread*) arg ;
  Array      aTx = arrayCreate (2048, TaxInfo) ;

  I64   seqLenMax ; oneStats (ti->ofIn, 'S', 0, &seqLenMax, 0) ;
  char *mLine = new (seqLenMax, char) ;
  I32   mScore = -(1<<30) ;
  char  mMap[128] ; { char  *s = ".-acgtnACGTN", *t = ".-tgcanTGCAN" ; while (*s) mMap[*s++] = *t++ ;  }

  int i, n = 0 ;
  oneGoto (ti->ofIn, 'S', ti->start) ;
  oneReadLine (ti->ofIn) ;
  while (ti->start <= ti->end)
    { ++ti->start ;
      I64      seqLen = oneLen(ti->ofIn) ;
      char    *seq = oneDNAchar(ti->ofIn), *qual ;
      oneWriteLine (ti->ofOut, 'S', seqLen, seq) ; // write back out the sequence
      Hash     h = hashCreate (8192) ;
      int      flags, txid, score ;
      TaxInfo *tx ;
      I64      nCigar, *i64cigar ;
      arrayMax(aTx) = 0 ;
      while (oneReadLine (ti->ofIn) && ti->ofIn->lineType != 'S')
	switch (ti->ofIn->lineType)
	  {
	  case 'Q': // save it for mismatch record
	    qual = oneString(ti->ofIn) ;
	    break ;
	  case 'J': // copy it
	    oneInt(ti->ofOut,0) = oneInt(ti->ofIn,0) ;
	    oneWriteLine (ti->ofOut,'J',oneLen(ti->ofIn),oneString(ti->ofIn)) ;
	    break ;
	  case 'B':
	    flags  = oneInt(ti->ofIn,0) ;
	    txid   = ti->taxid[oneInt(ti->ofIn,1)] ;
	    if (hashAdd (h, hashInt(txid), &i))
	      { tx = arrayp(aTx,i,TaxInfo) ;
		tx->txid = txid ; tx->count = 0 ; tx->bestScore = -(1<<30) ;
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
		int   nM = oneLen(ti->ofIn) ;
		char *md = oneString (ti->ofIn) ;
		// we actually want the edits in the read space, not the reference space
		// first reconstruct reference string, then map that to the read with cigar
		// see  https://vincebuffalo.com/notes/2014/01/17/md-tags-in-bam-files.html
		int iS = 0 ; // position in sequence, reference
		while (nCigar)
		  { int nS = bam_cigar_oplen(*i64cigar), op = bam_cigar_op(*i64cigar++) ; --nCigar ;
		    int nR = 0 ;
		    while (nM && (*md >= '0' && *md <= '9')) { nR = nR*10 + (*md++ - '0') ; --nM ; }
		    switch (bam_cigar_type(op))
		      {
		      case 0: break ; // do nothing
		      case 1: // DELETE (or REF_SKIP)
			while (nS--)
			  { if (nR || nM < 2 || *md++ != '^') die ("bad refmatch") ;
			    ++md ; nM -= 2 ; // eat up the deleted base
			  }
			break ;
		      case 2: // INSERT (or SOFT_CLIP)
			for (i = 0 ; i < nS ; ++i) mLine[iS++] = '-' ;
			break ;
		      case 3: // MATCH (or EQUAL or DIFF)
			while (nR < nS)
			  { while (nR) { mLine[iS++] = '.' ; --nR ; --nS ; }
			    // encode the edit
			    mLine[iS] = (qual[iS]<20)?tolower(*md++):toupper(*md++) ; ++iS ; --nS ;
			    
			    // 4*dna2indexConv[seq[iS++]] + dna2indexConv[*md++] ; --nS ;
			    while (nM && (*md >= '0' && *md <= '9')) {nR = nR*10 +(*md++ - '0'); --nM; }
			  }
			break ;
		      }
		  }
		if (iS != seqLen) die ("iS %d != seqLen %d", iS, seqLen) ;
		if (flags & BAM_FREVERSE) // need to reverse-complement mLine
		  { int i, j ;
		    for (i = 0, j = seqLen ; i < --j ; ++i)
		      { char t = mLine[i] ; mLine[i] = mMap[mLine[j]] ; mLine[j] = mMap[t] ; }
		    if (i == j) mLine[i] = mMap[mLine[i]] ;
		  }
	      }
	    break ;
	  default: break ;
	  }

      oneInt(ti->ofOut,0) = mScore ;
      oneWriteLine (ti->ofOut, 'M', seqLen, mLine) ;
      for (i = 0 ; i < arrayMax(aTx) ; ++i)
	{ tx = arrp(aTx,i,TaxInfo) ;
	  oneInt(ti->ofOut,0) = tx->txid ;
	  oneInt(ti->ofOut,1) = tx->bestScore ;
	  oneInt(ti->ofOut,2) = tx->count ;
	  oneWriteLine (ti->ofOut, 'T', 0, 0) ;
	}
      hashDestroy (h) ;
    }

  newFree (mLine, seqLenMax, char) ;
  arrayDestroy (aTx) ;
  return 0 ;
}

bool bamMake1read (char *bamOneFileName, char *outFileName)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *ofIn = oneFileOpenRead (bamOneFileName, 0, "bam", NTHREAD) ;
  if (!ofIn) { warn ("failed to open .1bam file %s", bamOneFileName) ; return false ; }
  if (!oneFileCheckSchemaText (ofIn, "P 3 seq\nD s 1 3 INT\nD m 1 6 STRING\n"))
    die ("input .1bam file must have s (score) and m (MD) record types") ;
  OneFile  *ofOut = oneFileOpenWriteNew (outFileName, schema, "read", true, NTHREAD) ;
  if (!ofOut) { oneFileClose(ofIn) ; warn ("failed to open %s", outFileName) ; return false ; }
  oneAddProvenance (ofOut, "onebam", "0.1", getCommandLine()) ;

  I64 nRef ; oneStats (ofIn, 'R', &nRef, 0, 0) ;
  int *taxid = new (nRef, int) ;
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

/*********************** end of file ***********************/

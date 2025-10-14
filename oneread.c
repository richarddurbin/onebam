/*  File: oneread.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: onebam functionality involving .1read files
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 14 01:41 2025 (rd109)
 * * Sep 30 04:08 2025 (rd109): added oneReader package, and dust()
 * Created: Wed Aug 13 14:10:49 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"
#include "seqio.h"  // for acgtCheck()

/********* OneReader API **********/

typedef struct {
  int      nThreads ;
  OneFile *ofIn, *ofOut ;
  int      prefixLen ;
  int      wCount[64] ; // for dust calculation
} OneReaderPrivate ;

OneReader *oneReaderCreate (char *fileName, int nThreads)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (fileName, schema, "read", nThreads) ;
  if (!of) return 0 ;

  OneReader *or = new0 (nThreads, OneReader) ;

  OneReaderPrivate *op = new0 (nThreads, OneReaderPrivate) ;
  or->private = (void*) op ;
  op->nThreads = nThreads ;
  op->ofIn = of ;
  
  if (oneGoto (of, 'P', 1) && oneReadLine (of))
    { or->namePrefix = oneString (of) ;
      op->prefixLen = oneLen (of) ;
    }

  or->taxonomy = taxonomyRead (of) ;
  
  if (!oneReadLine (of) || of->lineType != 'S')
    { oneGoto (of, 'S', 1)  ; oneReadLine (of) ; }

  int maxSeqLen, maxNameLen, maxTax ;
  oneReaderStats (or, 0, &maxSeqLen, 0, &maxNameLen, &maxTax) ; // must come after namePrefix
  
  int i ;
  for (i = 0 ; i < nThreads ; ++ i)
    { OneReader *ori = &or[i] ;
      ori->name = new0 (maxNameLen+1, char) ;
      if (op->prefixLen) strcpy (ori->name, or->namePrefix) ; // namePrefix is fixed and shared
      ori->taxid = new (maxTax, TaxID) ;
      ori->taxCount = new (maxTax, int) ;
      ori->taxBestScore = new (maxTax, int) ;
      if (i > 0)
	{ OneReaderPrivate *opi = op + i ;
	  ori->private = (void*) opi ;
	  opi->ofIn = of + i ;
	  opi->prefixLen = op->prefixLen ;
	  oneGoto (opi->ofIn, 'S', 1) ; oneReadLine (opi->ofIn) ; // sets up at start of first sequence
	}
    }

  return or ;
}

bool oneReaderGoto (OneReader *or, U64 i)
{ OneFile *of = ((OneReaderPrivate*)or->private)->ofIn ;
  if (i < 0) return false ;
  else if (i == 0) return oneGoto (of, 'S', 1) && oneReadLine (of) ;
  else if (oneGoto (of, 'S', i) && oneReadLine (of) && oneReaderNext (or)) return true ;
  else return false ;
}

void oneReaderStats (OneReader *or,
		     U64 *nReads,       // number of reads in the file
		     int *maxSeqLen,    // maximum sequence length
		     U64 *totSeqLen,    // total sequence length
		     int *maxNameLen,   // maximum name length
		     int *maxTax)       // max value of nTaxid
{ if (!or || !or->private) return ;
  OneFile *of = ((OneReaderPrivate*)or->private)->ofIn ;
  I64 a, b, c ;
  if (nReads || maxSeqLen || totSeqLen)
    { oneStats (of, 'S', &a, &b, &c) ;
      if (nReads)    *nReads = (U64) a ;
      if (maxSeqLen) *maxSeqLen = (int) b ;
      if (totSeqLen) *totSeqLen = (U64) c ;
    }
  if (maxNameLen)
    { oneStats (of, 'I', 0, &a, 0) ;
      *maxNameLen  = (int) a ;
      *maxNameLen += ((OneReaderPrivate*)or->private)->prefixLen ;
    }
  if (maxTax)
    { oneStatsContains (of, 'S', 'T', &a, 0) ;
      *maxTax = (int) a ;
    }
}

bool oneReaderNext (OneReader *or)
{
  OneReaderPrivate *op = (OneReaderPrivate*)or->private ;
  OneFile *of = op->ofIn ;
  if (!of->lineType || of->lineType == 'U') return false ; // at end of file or reads
  if (of->lineType != 'S') { printf ("**lineType %c\n", of->lineType) ; return false ; }
  or->seqLen = oneLen(of) ;
  or->seq = oneDNAchar(of) ; // NB we use the ONElib 'S' buffer
  or->nTax = 0 ;
  or->dustScore = 0 ;
  while (oneReadLine (of) && of->lineType != 'S' && of->lineType != 'U')
    switch (of->lineType)
      {
      case 'I':
	or->nameLen = op->prefixLen + oneLen(of) ;
	or->nameEnd = oneString(of) ;
	memcpy (or->name + op->prefixLen, or->nameEnd, oneLen(of)) ;
	break ;
      case 'N': // NB we are overwriting the ONElib 'S' buffer - OK
	{ int i = oneInt(of,0) ; char c = oneChar(of,1) ; int n = oneInt(of,2) ;
	  while (n--) or->seq[i++] = c ;
	}
	break ;
      case 'D':
	or->dustScore = oneInt(of,0) ;
	break ;
      case 'L':
	or->lca = oneInt(of,0) ;
	break ;
      case 'Q': // NB we are using and overwriting the ONElib 'Q' buffer - OK
	or->qual = (U8*)oneString(of) ;
	{ int i ; for (i = 0 ; i < or->seqLen ; ++i) or->qual[i] -= 33 ; }
	break ;
      case 'M':
	or->maxScore = oneInt(of,0) ;
	or->mLine = oneString(of) ;
	break ;
      case 'T':
	or->taxid[or->nTax] = oneInt(of,0) ;
	or->taxBestScore[or->nTax] = oneInt(of,1) ;
	or->taxCount[or->nTax] = oneInt(of,2) ;
	++or->nTax ;
	break ;
      default:
	die ("unknown linetype %c at line %lld sequence %d in .1read file %s",
	     of->lineType, of->line, oneObject(of,'S'), oneFileName (of)) ;
      }
  if (!or->dustScore) or->dustScore = (int)(0.5 + dust (or->seq, or->seqLen, 64, op->wCount)) ;
  return true ;
}

bool oneReaderWriteFile (OneReader *or, const char *outFileName)
{
  OneReaderPrivate *op = (OneReaderPrivate*)or->private ;
  if (!(op->ofOut = oneFileOpenWriteFrom (outFileName, op->ofIn, true, op->nThreads)))
    return false ;
  // need to change to binary, and to copy across the provenance
  oneAddProvenance (op->ofOut, "OneReader", VERSION, getCommandLine()) ;
  oneWriteLine (op->ofOut, 'P', strlen(or->namePrefix), or->namePrefix) ;
  int i ;
  for (i = 1 ; i < op->nThreads ; ++i) op[i].ofOut = op->ofOut + i ;
  return true ;
}

void oneReaderWrite (OneReader *or)
{
  OneFile *of = ((OneReaderPrivate*)or->private)->ofOut ;
  oneWriteLine (of, 'S', or->seqLen, or->seq) ;
  oneWriteLine (of, 'I', strlen(or->nameEnd), or->nameEnd) ;
  int i ;
  for (i = 0 ; i < or->seqLen ; ++i) // write exceptions for non-ACGT characters
    if (!acgtCheck[(int)or->seq[i]])
      { oneInt(of,0) = i ;
	char c = oneChar(of,1) = or->seq[i] ;
	I64 n = 1 ; while (++i < or->seqLen && or->seq[i] == c) n++ ; oneInt(of,2) = n ;
	oneWriteLine (of, 'N', 0, 0) ;
      }
  if (or->qual)
    { for (i = 0 ; i < or->seqLen ; ++i) or->qual[i] += 33 ;
      oneWriteLine (of, 'Q', or->seqLen, or->qual) ;
      for (i = 0 ; i < or->seqLen ; ++i) or->qual[i] -= 33 ;
    }
  if (or->dustScore) { oneInt(of,0) = or->dustScore ; oneWriteLine (of, 'D', 0, 0) ; }
  if (or->lca) { oneInt(of,0) = or->lca ; oneWriteLine (of, 'L', 0, 0) ; }
  oneInt(of,0) = or->maxScore ; oneWriteLine (of, 'M', or->seqLen, or->mLine) ;
  for (i = 0 ; i < or->nTax ; ++i)
    { oneInt(of,0) = or->taxid[i] ;
      oneInt(of,1) = or->taxBestScore[i] ;
      oneInt(of,2) = or->taxCount[i] ;
      oneWriteLine (of, 'T', 0, 0) ;
    }
}

void oneReaderWriteTx (OneReader *or, Taxonomy *tx, bool *txUsed)
{
  OneReaderPrivate *op = (OneReaderPrivate*)or->private ;
  taxonomyWrite (tx, op->ofOut+op->nThreads-1, txUsed) ; // LAST thread's handle for end of file
}

void oneReaderDestroy (OneReader *or)
{
  if (!or) return ;
  OneReaderPrivate *op = (OneReaderPrivate*)or->private ;
  int nThreads = op->nThreads ;

  int maxSeqLen, maxNameLen, maxTax ;
  oneReaderStats (or, 0, &maxSeqLen, 0, &maxNameLen, &maxTax) ; // must come here

  int i ;
  for (i = 0 ; i < nThreads ; ++i)
    { OneReader *ori = or + i ;
      newFree (ori->name, maxNameLen+1, char) ;
      newFree (ori->taxid, maxTax, int) ;
      newFree (ori->taxCount, maxTax, int) ;
      newFree (ori->taxBestScore, maxTax, int) ;
    }

  oneFileClose (op->ofIn) ; // closing the master closes all the slaves and releases line buffers
  if (op->ofOut) oneFileClose (op->ofOut) ;
  newFree (op, nThreads, OneReaderPrivate) ;
  newFree (or, nThreads, OneReader) ;
}

#ifdef READER_TEST1 // simple example

// compile with: gcc -DREADER_TEST1 -o reader oneread.c taxonomy.c ONElib.c merge.c seqio.c array.c dict.c utils.c -lz

int main (int argc, char *argv[])
{
  --argc ; ++argv ; // skip program name
  if (argc != 2) die ("usage: readReport <.1read file>  <i> // reports on i'th read in file") ;
  OneReader *or = oneReaderCreate (*argv, 1) ;
  if (!or) die ("failed to open .1read file %s", *argv) ;
  U64 i = atoi(*++argv) ;
  if (i <= 0) die ("second argument %s must be a positive integer", *argv) ;
  if (!oneReaderGoto (or, i)) die ("failed to locate to %llu", (long long unsigned)i) ;
  printf ("read %llu %s seqLen %d seq %s\n", (long long unsigned)i, or->name, or->seqLen, or->seq) ;
  printf ("maxScore %d dustScore %d mLine %s\n", or->maxScore, or->dustScore, or->mLine) ;
  printf   ("nTax %3d taxid   count  bestScore\n", or->nTax) ;
  int j ;
  for (j = 0 ; j < or->nTax ; ++j)
    printf ("        %8d   %5d  %8d\n", or->taxid[j], or->taxCount[j], or->taxBestScore[j]) ;
  oneReaderDestroy (or) ;
}

#endif

#ifdef READER_TEST2 // example using threading

// compile with: gcc -DREADER_TEST2 -o dustbin oneread.c taxonomy.c ONElib.c merge.c seqio.c array.c dict.c utils.c -lz

typedef struct {
  OneReader *or ;
  U64        start, end ;
  U64        dustBin[101] ;
} ThreadArg ;

static void *threadProcess (void *arg)
{
  ThreadArg *ta = (ThreadArg*)arg ;
  if (!oneReaderGoto (ta->or, ta->start)) die ("failed oneReaderGoto") ;
  while (ta->start++ < ta->end && oneReaderNext (ta->or))
    ++ta->dustBin[(int)(ta->or->dustScore)] ;
  return 0 ;
}

int main (int argc, char *argv[])
{
  int i, j, n ;
  --argc ; ++argv ; // skip program name
  if (argc != 2)
    die ("usage: dustBin <.1read file> <nThread> // report distribution of dust values") ;
  int nThread = atoi(argv[1]) ;
  if (nThread <= 0) die ("second argument %s must be a positive integer", argv[1]) ;
  OneReader *or = oneReaderCreate (argv[0], nThread) ;
  if (!or) die ("failed to open .1read file %s", argv[0]) ;

  pthread_t *threads = new (nThread, pthread_t) ;
  ThreadArg *ta      = new0 (nThread, ThreadArg) ;
  U64 nRead ; oneReaderStats (or, &nRead, 0, 0, 0, 0) ;
  if (!nRead) die ("no reads found in file %s", *argv) ;
  for (i = 0 ; i < nThread ; ++i)
    { ta[i].or  = or + i ;                     // equivalent to ta[i].or = &or[i]
      ta[i].start = (nRead * i) / nThread ;
      ta[i].end = (nRead * (i+1)) / nThread ;
    }
  for (i = 0 ; i < nThread ; ++i)              // create threads
    pthread_create (&threads[i], 0, threadProcess, &ta[i]) ;
  for (i = 0 ; i < nThread ; ++i)
    pthread_join (threads[i], 0) ;             // wait for threads to complete
  newFree (threads, nThread, pthread_t) ;

  // accumulate scores into the 0 level
  for (i = 1 ; i < nThread ; ++i) for (j = 0 ; j < 101 ; ++j) ta->dustBin[j] += ta[i].dustBin[j] ;

  // find the largest index in the histogram table with data
  for (n = 101 ; n-- ;) if (ta->dustBin[n]) break ;

  printf ("file %s containing %llu reads has dust score distribution:\n", *argv, nRead) ;
  for (j = 0 ; j <= n ; ++j)
    printf ("%d\t%.1f%%\t%llu\n", j, (100.0*ta->dustBin[j])/nRead, ta->dustBin[j]) ;

  newFree (ta, nThread, ThreadArg) ;
  oneReaderDestroy (or) ;
}
#endif

/*************************************************************************/

bool extractReads (char *inFileName, char *outFileName, TaxID lca)
{
  if (!lca) die ("extractReads can only extract based on LCA for now") ;
  
  OneReader *or = oneReaderCreate (inFileName, 1) ;
  if (!or) { warn ("failed to open .1read file %s", inFileName) ; return false ; }
  if (!outFileName)
    { int k = 0 ; while (inFileName[k] && inFileName[k] != '.') ++k ;
      inFileName[k] = 0 ;
      outFileName = new0 (k + 24, char) ;
      sprintf (outFileName, "%s-%d.1read", inFileName, lca) ;
    }
  if (!oneReaderWriteFile (or, outFileName))
    { warn ("failed to open .1read file %s to write to", outFileName) ; return false ; }

  int nIn = 0, nOut = 0 ;
  while (oneReaderNext (or))
    { ++nIn ;
      if (or->lca == lca)
	{ ++nOut ;
	  oneReaderWrite (or) ;
	}
    }
  oneReaderDestroy (or) ;

  printf ("extracted %d reads into %s from %d in %s\n", nOut, outFileName, nIn, inFileName) ;
  return true ;
}

/******* merge1read - merges arbitrarily many 1read files sorted on read name ********/
// uses Gene Myers merge package twice

#include "merge.h"

typedef struct {
  char    *fileName ;
  OneFile *of ;
  I64      seqLen ;
  U8      *dna ;      // keep pointer to 2bit form - requires less work
  char     name[1024] ;
} Input ;

bool yieldName (int t, void *arg, char **v, int *lcp) // here arg is an array of Input
{
  Input *in = &(((Input*)arg)[t]) ;
  if (in->of) { *v = in->name ; return true ; }
  else return false ;
}

bool yieldTx (int t, void *arg, I32 *v) // here arg is an array of OneFile
{
  OneFile *of = (((OneFile**)arg)[t]) ;
  if (of) { *v = oneInt(of,0) ; return true ; }
  else return false ;
}

void loadSequence (Input *in)
{
  if (in->of->lineType != 'S')
    die ("line type %c instead of expected S in %s line %lld",
	 in->of->lineType, in->fileName, in->of->line) ;
  in->seqLen = oneLen(in->of) ;
  in->dna = oneDNA2bit(in->of) ;
  while (oneReadLine (in->of) && (in->of->lineType != 'I')) ; // really should be next line...
  if (!in->of->lineType) die ("no I line found in %s", in->fileName) ;
  if (oneLen(in->of) > 1024)
    die ("read name length %d > 1024 file %s line %d", oneLen(in->of), in->fileName, in->of->line) ;
  strcpy (in->name, oneString(in->of)) ;
}

bool merge1read (char *outfile, int nIn, char **infiles)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;

  if (!outfile) outfile = "merged.1read" ;
  OneFile *ofOut = oneFileOpenWriteNew (outfile, schema, "read", true, 1) ;
  if (!ofOut) { warn ("failed to open %s to write", outfile) ; return false ; }
  oneAddProvenance (ofOut, "onebam", VERSION, getCommandLine()) ;

  int  i ;
  char t ;
  U64  nSin = 0, nTin = 0 ; // counters for input S and T lines

  // next build the inputs array and create the read merge object
  Input *inputs = new0 (nIn, Input) ;
  char  *prefix = 0 ;
  for (i = 0 ; i < nIn ; ++i)
    { Input *in = &inputs[i] ;
      in->fileName = infiles[i] ;
      if (!(in->of = oneFileOpenRead (in->fileName, 0, "read", 1)))
	die ("failed to open .1read file %s", in->fileName) ;
      if (!oneFileCheckSchemaText (in->of, "P 3 seq\nO S 1 3 DNA\nD I 1 6 STRING\n"
				   "D M 2 3 INT 6 STRING\nD T 3 3 INT 3 INT 3 INT\n"))
	die ("schema mismatch for input file %s", in->fileName) ;

      // for now just check that if there are prefixes then they are all the same
      if (oneGoto (in->of, 'P', 1))
	{ oneReadLine (in->of) ; // don't forget this!
	  if (!i) prefix = oneString (in->of) ;
	  else if (strcmp (prefix, oneString (in->of)))
	    die ("prefix mismatch %s in % versus %s in %s",
		 prefix, infiles[0], oneString(in->of), infiles[i]) ;
	}
      else if (prefix) die ("prefix %s in %s but none in %s", prefix, infiles[0], infiles[i]) ;

      if (oneGoto (in->of, 'S', 1) && oneReadLine(in->of)) // oneReadLine required to read S line
	{ loadSequence (in) ; ++nSin ; }
      else
	{ warn ("input .1read file %s appears empty", in->fileName) ; // it could actually be empty
	  in->of = 0 ;
	}
    }
  if (prefix) oneWriteLine (ofOut, 'P', strlen(prefix), prefix) ;
  
  Merge   *readMerge = mergeCreateString (nIn, inputs, yieldName) ;

  // do this now because inside the loop we will mergeRecreate(txMerge)
  OneFile **txOf      = new0 (nIn, OneFile*) ;
  Merge    *txMerge   = mergeCreateI32 (nIn, txOf, yieldTx) ;

  // now start the merge process
  int    nName, *nameList, nTx, *txList ;
//  char   lastOutName[1024] ; *lastOutName = 0 ;
  while ((nName = mergeNext (readMerge, &nameList)))
    { Input *in = &inputs[nameList[0]] ; // this is an example to copy the shared information from

      // first write the sequence record 'S'
      oneWriteLineDNA2bit (ofOut, 'S', in->seqLen, in->dna) ; // write the sequence record

      // then the name
      oneWriteLine (ofOut, 'I', strlen(in->name), in->name) ;
      
      if (nName == 1) // all very simple - go to next S line
	{ while ((t = oneReadLine(in->of)) && t != 'S')
	    { oneWriteLineFrom (ofOut, in->of) ; // copy line
	      if (t == 'T') ++nTin ;
	    }
	  if (t) { loadSequence (in) ; ++nSin ; }
	  else in->of = 0 ; // termination index for this input for readMerge
	}
      else // multiple input streams for the same read
	{ // first copy all lines up to M, which should be the same across all inputs
	  while ((t = oneReadLine(in->of)) && t != 'M')  // there should be an M line for each object
	    oneWriteLineFrom (ofOut, in->of) ; // copy line
	  if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[0]].fileName) ;
	  for (i = 1 ; i < nName ; ++i)
	    { txOf[i] = inputs[nameList[i]].of ; // will need this later, and it is more direct
	      while ((t = oneReadLine(txOf[i])) && t != 'M') ; // read these lines also
	      if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[i]].fileName) ;
	    }
	  txOf[0] = in->of ;
	  
	  // now find the highest scoring M line and write it out
	  I64 mMax = oneInt(in->of,0) ;
	  OneFile *ofM = in->of ;
	  for (i = 1 ; i < nName ; ++i)
	    if (oneInt(txOf[i],0) > mMax) { ofM = txOf[i] ; mMax = oneInt(ofM,0) ; }
	  oneWriteLineFrom (ofOut, ofM) ; // copy the highest scoring M line
	  
	  // move each txOf[] to the next line, which should be the first T line, and set up txMerge
	  for (i = 0 ; i < nName ; ++i)
	    if (oneReadLine(txOf[i]) == 'T') ++nTin ;
	    else die ("missing T line after M line input %s line %lld",
		      inputs[nameList[i]].fileName, (long long)txOf[i]->line) ;
	  mergeRecreate (txMerge, nName, txOf) ;

	  // and merge the taxids
	  while ((nTx = mergeNext (txMerge, &txList)))
	    { OneFile *ofTx = txOf[txList[0]] ;
	      oneInt(ofOut,0) = oneInt(ofTx,0) ; // copy the taxid
	      oneInt(ofOut,1) = oneInt(ofTx,1) ; // take the max of the best scores
	      oneInt(ofOut,2) = oneInt(ofTx,2) ; // sum the counts
	      for (i = 1 ; i < nTx ; ++i)
		{ ofTx = txOf[txList[i]] ;
		  if (oneInt(ofTx,1) > oneInt(ofOut,1)) oneInt(ofOut,1) = oneInt(ofTx,1) ;
		  oneInt(ofOut,2) += oneInt(ofTx,2) ;
		}
	      oneWriteLine (ofOut, 'T', 0, 0) ;
	      // now advance each line, and set txOf[] to 0 if not 'T'
	      for (i = 0 ; i < nTx ; ++i)
		if (oneReadLine (txOf[txList[i]]) == 'T') ++nTin ;
		else txOf[txList[i]] = 0 ; // terminator for the txMerge
	    }

	  // finally update the sequence and name in the inputs[] for the nameList 
	  for (i = 0 ; i < nName ; ++i)
	    { Input *in = &inputs[nameList[i]] ;
	      if (in->of->lineType) { loadSequence (in) ; ++nSin ; }
	      else { oneFileClose (in->of) ; in->of = 0 ; } // terminator for the nameMerge
	    }
	}
    }

  printf ("merged %d .1read files with %lld S lines and %lld T lines",
	  nIn, (long long) nSin, (long long) nTin) ;
  I64 nSout ; oneStats (ofOut, 'S', &nSout, 0, 0) ; 
  I64 nTout ; oneStats (ofOut, 'T', &nTout, 0, 0) ; 
  printf (" into %s with %lld S lines and %lld T lines\n",
	  oneFileName(ofOut), (long long)nSout, (long long)nTout) ;

  for (i = 0 ; i < nIn ; ++i) oneFileClose (inputs[i].of) ;
  oneFileClose (ofOut) ;

  mergeDestroy (readMerge) ;
  mergeDestroy (txMerge) ;

  newFree (inputs, nIn, Input) ;
  newFree (txOf, nIn, OneFile) ;

  return true ;
}

/******** report1read - not really functional/useful yet **********/

bool report1read (char *outFileName, char *readFileName)
{
  int i, j, k ;
  
  // open the output file
  FILE *fOut ;
  if (outFileName)
    { fOut = fopen (outFileName, "w") ;
      if (!fOut) die ("failed to open output file %s", outFileName) ;
    }
  else
    fOut = stdout ;

  // open the 1read file
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *of = oneFileOpenRead (readFileName, 0, "read", 1) ;
  if (!of) die ("failed to open .1read file %s", readFileName) ;

  // set up accumulation arrays - initialise to 0 by setting first element 0
  U64  lengths[256] = {0} ;
  U64  start[15][4][4] = {0}, end[15][4][4] = {0}, mid[4][4] = {0} ;
  U64 *edits = new0 (256, U64) ;
  U64 *scores = new0 (256, U64) ;
  I32  minScore = I32MAX, maxScore = -I32MAX ;
  int  maxEdit = 0 ;

  int  map[256] ; for (i = 0 ; i < 256 ; ++i) map[i] = -1 ;
  map['A']=map['a']=0 ; map['C']=map['c']=1 ; map['G']=map['g']=2 ; map['T']=map['t']=3 ;
  char index2char[] = "acgt" ;
  
  // process records
  U64 nRecord = 0 ;
  if (!oneGoto (of, 'S', 1)) die ("no objects in 1read file %s", readFileName) ;
  oneReadLine (of) ; // read the 'S' line
  while (of->lineType == 'S')
    { ++nRecord ;
      int seqLen = oneLen(of) ;
      ++lengths[seqLen] ;
      char *seq = oneDNAchar(of) ;

      // process the M line
      while (oneReadLine (of) && of->lineType != 'M') if (of->lineType == 'S') continue ;
      int score = oneInt(of,0) ;
      if (score > maxScore) maxScore = score ;
      if (score < minScore) minScore = score ;
      if (score > 127 || score < -128) warn ("score %d record %llu", score, nRecord) ;
      else ++scores[score+128] ;
      char *m = oneString(of) ;
      int nEdit = 0 ;
      for (i = 0 ; i < seqLen ; ++i)
	{ int from, to = map[seq[i]] ;
	  if (m[i] == '.') from = to ;
	  else { from = map[m[i]] ; ++nEdit ; }
	  if (from >= 0 && to >= 0)
	    { if (i < 15) ++start[i][from][to] ;
	      else if (seqLen - 1 - i < 15) ++end[seqLen-1-i][from][to] ;
	      else ++mid[from][to] ;
	    }
	}
      ++edits[nEdit] ;
      if (nEdit > maxEdit) maxEdit = nEdit ;

      // for now skip all the T lines - should process them
      while (oneReadLine (of) && of->lineType != 'S') ;
    }
  oneFileClose (of) ;

  // report
  fprintf (fOut, "RECORDS %llu\n", nRecord) ;
  if (nRecord)
    { int min = 0, max = 0, n50 ;
      U64 sum = 0, psum = 0 ;
      for (i = 0 ; i < 256 ; ++i)
	{ sum += i*lengths[i] ;
	  if (lengths[i] && !min) min = i ;
	  if (lengths[i]) max = i ;
	}
      for (i = 0 ; i < 256 ; ++i) { psum += i*lengths[i] ; if (psum*2 < sum) n50 = i+1 ; }
      double mean = sum/(double)nRecord ;
      fprintf (fOut, "LENGTHS min %d max %d mean %.1f n50 %d\n", min, max, mean, n50) ;
      for (i = 255 ; i >= 0 ; --i)
	if (scores[i]) fprintf (fOut, "SCORES %4d %llu\n", i-128, scores[i]) ;
      fprintf (fOut, "SCORE_MIN_MAX %d %d\n", minScore, maxScore) ;
      fprintf (fOut, "EDITS") ;
      for (i = 0 ; i <= maxEdit ; ++i) fprintf (fOut, " %llu", edits[i]) ;
      fprintf (fOut, "\n") ;
      double ct, ga, other ;
      for (i = 0 ; i < 15 ; ++i)
	{ ct = ga = other = 0. ;
	  double nC = start[i][1][0] + start[i][1][1] + start[i][1][2] + start[i][1][3] ;
	  if (nC) ct = start[i][1][3] / nC ;
	  double nG = start[i][2][0] + start[i][2][1] + start[i][2][2] + start[i][2][3] ;
	  if (nG) ga = start[i][2][0] / nG ;
	  double nX = start[i][0][0] + start[i][0][1] + start[i][0][2] + start[i][0][3] +
	    nC + nG + start[i][3][0] + start[i][3][1] + start[i][3][2] + start[i][3][3]
	    - start[i][1][3] - start[i][2][0] ;
	  if (nX) other = (nX-start[i][0][0]-start[i][1][1]-start[i][2][2]-start[i][3][3]) / nX ;
	  fprintf (fOut, "START %2d C>T %.3f G>A %.3f OTHER %.3f\n", i+1, ct, ga, other) ;
	}
      { ct = ga = other = 0. ;
	double nC = mid[1][0] + mid[1][1] + mid[1][2] + mid[1][3] ;
	if (nC) ct = mid[1][3] / nC ;
	double nG = mid[2][0] + mid[2][1] + mid[2][2] + mid[2][3] ;
	if (nG) ga = mid[2][0] / nG ;
	double nX = mid[0][0] + mid[0][1] + mid[0][2] + mid[0][3] +
	  nC + nG + mid[3][0] + mid[3][1] + mid[3][2] + mid[3][3]
	  - mid[1][3] - mid[2][0] ;
	if (nX) other = (nX-mid[0][0]-mid[1][1]-mid[2][2]-mid[3][3]) / nX ;
	fprintf (fOut, "MID      C>T %.3f G>A %.3f OTHER %.3f\n", ct, ga, other) ;
      }
      for (i = 0 ; i < 15 ; ++i)
	{ ct = ga = other = 0. ;
	  double nC = end[i][1][0] + end[i][1][1] + end[i][1][2] + end[i][1][3] ;
	  if (nC) ct = end[i][1][3] / nC ;
	  double nG = end[i][2][0] + end[i][2][1] + end[i][2][2] + end[i][2][3] ;
	  if (nG) ga = end[i][2][0] / nG ;
	  double nX = end[i][0][0] + end[i][0][1] + end[i][0][2] + end[i][0][3] +
	    nC + nG + end[i][3][0] + end[i][3][1] + end[i][3][2] + end[i][3][3]
	    - end[i][1][3] - end[i][2][0] ;
	  if (nX) other = (nX-end[i][0][0]-end[i][1][1]-end[i][2][2]-end[i][3][3]) / nX ;
	  fprintf (fOut, "END   %2d C>T %.3f G>A %.3f OTHER %.3f\n", i+1, ct, ga, other) ;
	}
      for (i = 0 ; i < 4 ; ++i)
	for (j = 0 ; j < 4 ; ++j)
	  { fprintf (fOut, " COUNTS %c>%c", index2char[i], index2char[j]) ;
	    fprintf (fOut, " MID %llu START", mid[i][j]) ;
	    for (k = 0 ; k < 15 ; ++k) fprintf (fOut, " %llu", start[k][i][j]) ;
	    fprintf (fOut, " END") ;
	    for (k = 15 ; k-- ; ) fprintf (fOut, " %llu", end[k][i][j]) ;
	    fprintf (fOut, "\n") ;
	  }
    }
  
  if (fOut) fclose (fOut) ;
  newFree (edits,256,U64) ;
  newFree (scores,256,U64) ;
  return true ;
}

/******** dust() reimplemented from Chenxi Zhou's reimplementation of published sdust ********/
// NB my version differs slightly - my window is the number of words (triplets) considered
//    whereas Chenxi 

static U8 dustMap[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

#define WLEN  3			// word length
#define WTOT  (1<<(WLEN<<1))
#define WMASK (WTOT - 1)

double dust (const char *seq, int seqLen, int window, int *wCount)
{
  static int  lastWindow = 0 ;
  static int  wCount0[WTOT], *wSeq ;

  if (window < WLEN) die ("dust window %d must be >= word length %d", window, WLEN) ;
  if (window != lastWindow)
    { if (lastWindow) newFree (wSeq, lastWindow, int) ;
      wSeq = new0 (window, int) ;
      lastWindow = window ;
    }
  if (!wCount) wCount = wCount0 ;
  memset(wCount, 0, WTOT*sizeof(int)) ;

  I64 score = 0, maxScore = 0 ;
  int i, j, t, n = -WLEN ;
  for (i = 0 ; i < seqLen ; ++i)
    { int b = dustMap[(int)seq[i]] ;
      if (b > 3) continue ; // ignore Ns
      t = (t << 2 | b) & WMASK ;
      if (++n >= 0)
	{ int k = n % window ;
	  if (n >= window)
	    { int x = wSeq[k] ;
	      if (wCount[x]) score -= --wCount[x] ;
	      score += wCount[t]++ ;
	      if (score > maxScore) maxScore = score ;
	    }
	  else
	    score += wCount[t]++ ;
	  wSeq[k] = t ;
	}
    }
  if (n >= window) return (200.0 * maxScore) / (window * (window-1)) ;
  else return (200.0 * score) / (n * (n+1)) ;
}

#ifdef DUST_TEST

// compile with gcc -DDUST_TEST -o dust oneread.c seqio.c merge.c utils.c ONElib.c -lz

#include "seqio.h"

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (argc != 2) die ("usage: dust <window> <seqfile>") ;
  int window = atoi (*argv++) ;
  SeqIO *sio = seqIOopenRead (*argv, 0, false) ;
  if (!sio) die ("failed to open sequence file %s", *argv) ;
  while (seqIOread (sio))
    printf ("%s\t%d\t%.3f\n",
	    sqioId(sio), (int)sio->seqLen, dust(sqioSeq(sio), sio->seqLen, window, 0)) ;
  seqIOclose (sio) ;
}
 
#endif // DUST_TEST

/********************* end of file ***********************/

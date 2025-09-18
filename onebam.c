/*  File: onebam.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  2 18:36 2025 (rd109)
 * Created: Wed Jul  2 10:18:19 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"

// globals
int        NTHREAD = 8 ;

// routine to convert sequence file to numbers
static bool numberSeq (char* inName, char *out1seqName, char *outFqName, bool isNames) ;

// reference
ExternalReference *reference = 0 ;
static bool makeAccTax (char *refname, char *tsvname) ;
static bool readReference (char *refname) ;

// extract stem from input filename
char *derivedName (char *inName, char *tag)
{
  int   nameLen = strlen(inName) ;
  char *stem = new (nameLen+1, char) ;
  strcpy (stem, inName) ;
  char *s = stem + nameLen ;
  if (s - stem > 3 && !strcmp (s-3, ".gz")) { s = s-3 ; *s = 0 ; }
  while (s > stem)
    if (*--s == '.') { *s = 0 ; break ; }
    else if (*s == '/') break ; // just add the tag
  char *ret = fnameTag(stem,tag) ;
  newFree (stem, nameLen+1, char) ;
  return ret ;
}

#define HIDE

static char usage[] =
  "Usage: onebam <command> <options>* <args>*\n"
  "  options and arguments depend on the command\n"
  "  commands with arguments, each followed by their options:\n"
  "    bam21read <XX.1acctax> <YY.bam>  convert read-sorted BAM/SAM/CRAM file to .1read\n"
  "      -o <ZZ.1read>                    output file - default is YY.1read for input YY.bam\n"
  "    merge1read [*.1read]+         merge .1read files\n"
  "      -o <ZZ.1read>                    output file - default is 'combined.1read'\n"
  "    report1read <XX.1read>        report on contents of .1read file\n"
  "      -o <ZZ.report>                   output to named file rather than stdout\n"
  "    makeAccTax <XX.tsv>           convert acc2taxid tab-separated file to .1acctax\n"
  "      -o <ZZ.1acctax>                  output file - default is XX.1acctax for input XX.*\n"
  "    bam21bam <XX.bam>             convert BAM/SAM/CRAM file to .1bam\n"
  "      -o <ZZ.1bam>                     output file - default is XX.1bam for input XX.bam\n"
  "      -taxid <XX.tsv>                  file with lines <acc>\\ttaxid\\n sorted on acc\n"
  "         NB you must give taxids if you plan to use your .1bam to make .1read\n"              
  "      -aux <tag>:<fmt> <char>          record BAM tag with OneCode <char>, e.g. -aux AS:i s\n"
  "         NB you must use '-aux AS:i s -aux MD:Z m' if you plant to use your .1bam to make .1read\n"
  "      -names                           keep sequence names [default to drop names]\n"
  "      -cramref <file_name|URL>         reference for CRAM - needed to read cram\n"
  "    1bam21read <XX.1bam>          convert .1bam file to .1read (simplified information per read)\n"
  "      -T <nthreads>                    number of threads [8]\n"
  "      -o <ZZ.1read>                    output file - default is XX.1read for input XX.1bam\n"
#ifndef HIDE
  "    1bam2bam <XX.1bam>            back-convert .1bam to .bam\n"
  "      -o <ZZ.bam>                      output file - default is XX.bam for input XX.1bam\n"
  "    makebin <taxid.tsv> <XX.bam>  make fixed-width .alb and .txb binary files from bam file\n"
  "      -oTxb <ZZ.txb>                   binary taxid output file - default XX.txb from XX.bam\n"
  "      -oAlb <ZZ.alb>                   binary alignment output file - default XX.alb from XX.bam\n"
  "      -prefixLen <n>                   ignore first n characters in names [0]\n"
  "      -maxChars <n>                    write the next n chars of the name, or all if fewer [48]\n"
  "      -maxEdit <n>                     maximum number of edits for read to be in .alb file [8]\n"
  "    albReport <XX.alb>            report edit (damage) distribution etc. from .alb file\n"
  "      -o <filename>                    write the report to filenam - default is stdout\n"
  "    numberSeq <XX.fq[.gz]>        make .1seq file from fq, plus new fq file with ints for names\n"
  "         NB will read and process fasta[.gz] or BAM/CRAM or even 1seq, as well as fastq[.gz]\n"
  "      -oSeq <ZZ.1seq>                  .1seq output file - default is XX.1seq for input XX.fq\n"
  "      -oFq <ZZ-i.fq[.gz]>              .fq output file - default is XX-i.fq.gz for input XX.fq\n"
  "      -names                           keep sequence names in .1seq [default to drop names]\n"
#endif
  ;
  
int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }
  char *command = *argv ; ++argv ; --argc ;

  char *outFileName = 0 ;
  char *taxidFileName = 0 ;
  bool  isNames = false ;

  if (!strcmp (command, "bam21read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bam21read option %s - run without args for usage", *argv) ;
      if (argc != 2)
	die ("onebam bam21read needs 2 not %d args; run without args for usage", argc) ;
      if (!bam21read (argv[1], outFileName, *argv))
	die ("failed to convert bam file %s to .1read file", *argv) ;
    }
  else if (!strcmp (command, "merge1read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam merge1read option %s - run without args for usage", *argv) ;
      if (argc < 2)
	die ("onebam merge1read needs at least 2 not %d args; run without args for usage", argc) ;
      if (!merge1read (outFileName, argc, argv))
	die ("failed to merge .1read files") ;
    }
  else if (!strcmp (command, "report1read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam report1read option %s - run without args for usage", *argv) ;
      if (argc < 1)
	die ("onebam report1read needs one not %d args; run without args for usage", argc) ;
      if (!report1read (outFileName, *argv))
	die ("failed to generate report from file %s", *argv) ;
    }
  else if (!strcmp (command, "makeAccTax"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam makeAccTax option %s - run without args for usage", *argv) ;
      if (argc < 1)
	die ("onebam makeAccTax needs one not %d args; run without args for usage", argc) ;
      if (!makeAccTax (outFileName, *argv))
	die ("failed to make .1acctax file from %s", *argv) ;
    }
  else if (!strcmp (command, "bam21bam"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-taxid") && argc > 1)
	  { taxidFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-names"))
	  { isNames = true ; ++argv ; --argc ; }
	else if (!strcmp (*argv, "-aux") && argc > 2)
	  { auxAdd (argv[2], argv[1]) ; argv += 3 ; argc -= 3 ; }
	else if (!strcmp (*argv, "-cramref") && argc > 1)
	  { setCramReference (argv[1]) ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bam21bam option %s - run without args for usage", *argv) ;
      if (argc != 1) die ("onebam bam21bam needs 1 not %d args; run without args for usage", argc) ;
      if (!bam21bam (*argv, outFileName, taxidFileName, isNames))
	die ("failed to convert bam file %s", *argv) ;
    }
  else if (!strcmp (command, "1bam21read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-T") && argc > 1) { NTHREAD = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
	else if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam make1read option %s - run without args for usage", *argv) ;
      if (argc != 1) die ("onebam make1read needs 1 not %d args; run without args for usage", argc) ;
      if (!bamMake1read (*argv, outFileName))
	die ("failed to convert .1bam file %s to .1read", *argv) ;
    }
  else if (!strcmp (command, "numberSeq"))
    { char *outFqName = 0 ;
      while (argc && **argv == '-')
	if (!strcmp (*argv, "-oSeq") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-oFq") && argc > 1)
	  { outFqName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-names"))
	  { isNames = true ; ++argv ; --argc ; }
	else die ("unknown onebam numberSeq option %s - run without args for usage", *argv) ;
      if (argc != 1) die ("onebam numberSeq needs 1 not %d args; run without args for usage", argc) ;
      if (!numberSeq (*argv, outFileName, outFqName, isNames))
	die ("failed to convert sequence file to .1seq and numbered .fq.gz files") ;
    }
  else if (!strcmp (command, "makebin"))
    { int   maxEdit = 8 ;
      int   prefixLen = 0 ;
      int   maxChars = 48 ;
      char *outAlbName = 0 ;
      while (argc && **argv == '-')
	if (!strcmp (*argv, "-oTxb") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-oAlb") && argc > 1)
	  { outAlbName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-prefixLen") && argc > 1)
	  { prefixLen = atoi (argv[1]) ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-maxChars") && argc > 1)
	  { maxChars = atoi (argv[1]) ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-maxEdit") && argc > 1)
	  { maxEdit = atoi (argv[1]) ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam makebin option %s - run without args for usage", *argv) ;
      if (argc != 2) die ("onebam makebin needs 2 not %d args; run without args for usage", argc) ;
      if (!makeBin (argv[1], outFileName, outAlbName, argv[0], maxEdit, prefixLen, maxChars))
	die ("failed to make binary file from bam file") ;
    }
  else if (!strcmp (command, "albReport"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam albReport option %s - run without args for usage", *argv) ;
      if (argc != 1) die ("onebam albReport needs 1 not %d args; run without args for usage",argc) ;
      albReport (*argv, outFileName) ;
    }
  else
    die ("unknown onebam command %s - run without arguments for usage", command) ;

  fprintf (stderr, "Total: ") ; timeTotal (stderr) ;
  exit (0) ;
}

/************************ numberSeq ***********************************/

#include "seqio.h"

static bool numberSeq (char* inName, char *out1seqName, char *outFqName, bool isNames)
{
  SeqIO *siIn = seqIOopenRead (inName, 0, true) ;
  if (!siIn) return false ;
  
  if (!outFqName) outFqName = derivedName (inName, "fqn.gz") ;
  SeqIO *siOut = seqIOopenWrite (outFqName, FASTQ, 0, 1) ;
  if (!siOut) { seqIOclose (siIn) ; return false ; }

  if (!out1seqName) out1seqName = derivedName (inName, "1seq") ;
  SeqIO *siSeq = seqIOopenWrite (out1seqName, ONE, 0, 1) ;
  if (!siSeq) { seqIOclose (siIn) ; seqIOclose (siOut) ; return false ; }

  I64  nSeq = 0, tot = 0 ;
  char ibuf[16] ;
  while (seqIOread (siIn))
    { sprintf (ibuf, "%lld", (long long) ++nSeq) ; // NB the increment here of nSeq
      seqIOwrite (siOut, ibuf, 0, siIn->seqLen, sqioSeq(siIn), sqioQual(siIn)) ;
      seqIOwrite (siSeq, isNames?sqioId(siIn):0, 0, siIn->seqLen, sqioSeq(siIn), sqioQual(siIn)) ;
      tot += siIn->seqLen ;
    }

  seqIOclose (siIn) ;
  seqIOclose (siOut) ;
  seqIOclose (siSeq) ;
  fprintf (stderr, "processed %lld sequences total length %lld\n", (long long) nSeq, (long long) tot) ;
  return true ;
}

/*********************** reference package *****************************/

#include "dict.h"
#include <ctype.h>

typedef struct {
  U64 start ;
  U32 count, taxid ;
} Block ;

static int blockCompare (const void *a, const void *b)
{ Block *ba = (Block*)a, *bb = (Block*)b ;
  I64 diff = (I64)ba->start - (I64)bb->start ;
  if (diff > 0) return 1 ; else if (diff < 1) return -1 ; else return 0 ;
}

static U64 accParse (char *acc) // truncates acc to its prefix before terminal digits
{
  char *s = acc ; while (*s) ++s ; // go to end of acc
  while (--s >= acc && *s >= '0' && *s <= '9') ; // go to last char before terminal digits (if any)
  char *term = ++s ;
  int n = 0 ; while (*s >= '0' && *s <= '9') n = n * 10 + (*s++ - '0') ;
  *term = 0 ; // truncate
  return n ;
}

DICT *compareDict ;
static int dictNameCompare (const void *a, const void *b)
{ char *sa = dictName(compareDict, *(int*)a), *sb = dictName(compareDict, *(int*)b) ;
  return strcmp (sa, sb) ;
}

// my own code to parse the acc2taxid file to deal with variants
static bool readAccLine (FILE *in, long long nLine, char *accBuf, int *tid)
{
  char c, *s = accBuf ;
  int n = 0 ;
  *tid = 0 ;
  // Skip blank lines and leading whitespace
  do {
    c = getc(in);
    if (c == EOF) return false;
  } while (isspace(c) && c != '\n');
  if (c == '\n') return false; // blank line
  // Start parsing the accession
  *s++ = c;
  n = 1;
  while ((c = getc(in)) != EOF && !isspace(c))
    { *s++ = c ;
      if (++n == 64)
	{ warn ("accession name is too long - line %lld", nLine) ;
	  while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
	}
      *s = 0 ;
    }
  if (!strcmp (accBuf, "accession")) // ignore this line
    { warn ("header line in acc2tax file line %lld", nLine) ;
      while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
    }
  if (c != '\t') die ("missing TAB after 1st field - line %lld", nLine) ;
  s = accBuf ;
  while ((c = getc(in)) != EOF && !isspace(c) && *s)
    if (c != *s++)
      { warn ("mismatch of first two fields - line %lld", nLine) ;
	while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
      }
  if (c != '.')
    { warn ("2nd field is not versioned 1st field - line %lld", nLine) ;
      while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
    }
  while ((c = getc(in)) != EOF && !isspace(c))
    if (c < '0' || c > '9')
      { warn ("version is not a number - line %lld", nLine) ;
	while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
      }
  if (c != '\t')
    { warn ("missing TAB after 2nd field - line %lld", nLine) ;
      while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
    }
  while ((c = getc(in)) != EOF && !isspace(c))
    if (c < '0' || c > '9')
      { warn ("taxid is not a number - line %lld", nLine) ;
	while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
      }
    else
      *tid = *tid * 10 + (c - '0') ;
  if (c != '\t' && c != '\n')
    { warn ("missing TAB after taxid - line %lld", nLine) ;
      while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
    }
  if (c != '\n')
    { warn ("error %d - line %lld", ferror(in), nLine) ;
      while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
    }
  return true ;
}

static bool makeAccTax (char *accTaxName, char *tsvName)
{
  if (!accTaxName) accTaxName = derivedName (tsvName, "1acctax") ;
  printf ("making reference %s from %s\n", accTaxName, tsvName) ;

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ; // check can open output
  OneFile *of = oneFileOpenWriteNew (accTaxName, schema, "acctax", true, 1) ;
  if (!of) die ("failed to open reference database ONEfile %s", accTaxName) ;

  FILE *in = fzopen (tsvName, "r") ;
  if (!in) die ("failed to open acc2tax file %s\n", tsvName) ;

  DICT  *pDict = dictCreate (1<<20) ;
  Array  pArray = arrayCreate (1<<20, Array) ;
  
  long long nLine = 0, nAcc = 0, nRawBlock = 0 ;
  char accBuf[64], accLast[64] ; accBuf[63] = 0 ; accLast[0] = 0 ;
  U64  index, indexLast = 0 ;
  U32  k, kLast = -1 ;
  U32  tid, tidLast ;
  while (!feof (in))
    { int tid ;
      if (!readAccLine (in, ++nLine, accBuf, &tid)) continue ;
      ++nAcc ;
      int pLen ;
      U64 index = accParse (accBuf) ;
      Block *b = 0 ;
      if (!strcmp (accBuf, accLast)) // same as previous - common so direct check is worth while
	k = kLast ;
      else
	{ if (dictAdd (pDict, accBuf, &k)) // a new prefix
	    { if (k != arrayMax(pArray)) die ("logic error in buildReference") ;
	      Array ak = array(pArray, k, Array) = arrayCreate (4, Block) ;
	      b = arrayp(ak,0,Block) ; // make a new block
	      b->taxid = tid ;
	      b->start = index ;
	      b->count = 1 ;
	      ++nRawBlock ;
	    }
	  strcpy (accLast, accBuf) ;
	}
      if (!b) // i.e. same as previous or found in the dict
	{ Array ak = arr(pArray, k ,Array) ;
	  Block *b = arrp(ak,arrayMax(ak)-1,Block) ; // pointer to existing final block
	  if (tid == b->taxid && index == b->start + b->count)
	    ++b->count ; // extend block
	  else
	    { b = arrayp(ak,arrayMax(ak),Block) ; // make a new block
	      b->taxid = tid ;
	      b->start = index ;
	      b->count = 1 ;
	      ++nRawBlock ;
	    }
	}
      kLast = k ;
      //      if (!(nLine % 1000000))
      //        { printf ("record %lld: ", nLine) ;
      //	  printf ("%d prefixes %lld raw blocks: ", dictMax(pDict), nRawBlock) ;
      //          printf ("prefix %s index %llu k %d\n", accBuf, (unsigned long long) index, k) ;
      //        }
    }
  fclose (in) ;

  fprintf (stderr, "read %lld accessions for %d prefixes %lld raw blocks: ",
	   (long long)nAcc, dictMax(pDict), (long long) nRawBlock) ;
  timeUpdate (stderr) ;

  int i, j ;
  int *iDict = new(dictMax(pDict),int) ;
  for (i = 0 ; i < dictMax(pDict) ; ++i) iDict[i] = i ;
  compareDict = pDict ; qsort (iDict, dictMax(pDict), sizeof(int), dictNameCompare) ;
  *accLast = 0 ;
  for (i = 0 ; i < dictMax(pDict) ; ++i)
    { int k = iDict[i] ;
      int d = 0 ; char *s = dictName(pDict, k), *t = accLast ;
      while (*s == *t) { ++s ; ++t ; ++d ; } // NB they can't be the same, so !(*s == *t == 0)
      Array ak = arr(pArray,k,Array) ;
      arraySort (ak, blockCompare) ;
      // should compress here - also remove duplicates at the same time
      Block *b = arrp(ak,0,Block) ;
      oneInt(of,0) = b->taxid ; oneInt(of,1) = b->start ; oneInt(of,2) = b->count ;
      oneInt(of,3) = d ;
      oneWriteLine (of, 'A', strlen(s), s) ;
      ++b ;
      oneInt(of,3) = d + strlen(s) ;
      strcpy (t, s) ; // update accLast to be full name
      for (j = 1 ; j < arrayMax(ak) ; ++j, ++b)
	{ oneInt(of,0) = b->taxid ; oneInt(of,1) = b->start ; oneInt(of,2) = b->count ;
	  oneWriteLine (of, 'A', 0, "") ;
	}
    }

  oneFileClose (of) ;
  return true ;
}
/********** threaded code to read in the reference ************/

// *** THIS CODE IS OUT OF DATE ***

typedef struct {
  OneFile   *of ;
  ExternalReference *ref ;
  I64        start, end ;
} RefReadThread ;

static void *readRefThread (void *arg)
{
  RefReadThread *ti = (RefReadThread*) arg ;

  oneReadLine (ti->of) ;
  while (ti->start <= ti->end)
    { ti->ref->len[ti->start] = oneInt(ti->of,0) ;
      ti->ref->taxid[ti->start] = (I32) oneInt(ti->of,1) ;
      ++ti->start ;
      oneReadLine (ti->of) ;
      if (ti->of->lineType != 'I')
	{ 
	}
      if (ti->of->lineType != 'A') die ("read failure at readRef %lld", (long long)ti->start) ;
    }
  return 0 ;
}

static bool readReference (char *fname)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile *of = oneFileOpenRead (fname, schema, "ref", NTHREAD) ;
  if (!of) return false ;

  fprintf (stderr, "opened oneFile: ") ; timeUpdate (stderr) ;

  reference = new (1, ExternalReference) ;
  I64 max, total ;
  oneStats (of, 'A', &reference->n, &max, &total) ;
  reference->len = new (reference->n, I64) ;
  reference->taxid = new (reference->n, I32) ;
  oneStats (of, 'I', &reference->nAcc, 0, &reference->fixBufSize) ;
  reference->acc = new (reference->nAcc, Accession) ;
  reference->fixBuf = new (reference->fixBufSize, char) ;

  pthread_t     *threads = new (NTHREAD, pthread_t) ;
  RefReadThread *ti = new0 (NTHREAD, RefReadThread) ;
  int i ;
  for (i = 0 ; i < NTHREAD ; ++i)
    { ti[i].ref = reference ;
      ti[i].of = of + i ;
      ti[i].start = 1 + (reference->n * i) / NTHREAD ;
      ti[i].end = (reference->n * (i+1)) / NTHREAD ;
      oneGoto (ti[i].of, 'A', ti->start) ;
    }
  for (i = 0 ; i < NTHREAD ; ++i) // create threads
    pthread_create (&threads[i], 0, readRefThread, &ti[i]) ;
  for (i = 0 ; i < NTHREAD ; ++i)
    pthread_join (threads[i], 0) ; // wait for threads to complete
  newFree (threads, NTHREAD, pthread_t) ;
  newFree (ti, NTHREAD, RefReadThread) ;
  
  oneFileClose (of) ;
  fprintf (stderr,"read reference %s with %lld accessions: ", fname, (long long) reference->n) ;
  timeUpdate (stderr) ;
  return true ;
}

/*************************** end of file ****************************/

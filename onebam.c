/*  File: onebam.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  3 23:31 2025 (rd109)
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
static void  buildReference (char *refname, char *tsvname) ;
static bool  readReference (char *refname) ;

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

static char usage[] =
  "Usage: onebam <command> <options>* <args>*\n"
  "  options and arguments depend on the command\n"
  "  commands with arguments, each followed by their options:\n"
  "    bam21bam <XX.bam>             convert BAM/SAM/CRAM file to .1bam\n"
  "      -o <ZZ.1bam>                     output file - default is XX.1bam for input XX.bam\n"
  "      -taxid <XX.tsv>                 file with lines <acc>\\ttaxid\\n sorted on acc\n"
  "         NB you must give taxids if you plan to use your .1bam to make .1read\n"              
  "      -aux <tag>:<fmt> <char>          record BAM tag with OneCode <char>, e.g. -aux AS:i s\n"
  "         NB you must use '-aux AS:i s -aux MD:Z m' if you plant to use your .1bam to make .1read\n"
  "      -names                           keep sequence names [default to drop names]\n"
  "      -cramref <file_name|URL>         reference for CRAM - needed to read cram\n"
  "    1bam21read <XX.1bam           convert .1bam file to .1read (simplified information per read)\n"
  "      -T <nthreads>                    number of threads [8]\n"
  "      -o <ZZ.1read>                    output file - default is XX.1read for input XX.1bam\n"
  "    makebin <taxid.tsv> <XX.bam>  make fixed-width .alb and .txb binary files from bam file\n"
  "      -oTxb <ZZ.txb>                   binary taxid output file - default XX.txb from XX.bam\n"
  "      -oAlb <ZZ.alb>                   binary alignment output file - default XX.alb from XX.bam\n"
  "      -prefixLen <n>                   ignore first n characters in names [0]\n"
  "      -maxChars <n>                    write the next n chars of the name, or all if fewer [48]\n"
  "      -maxEdit <n>                     maximum number of edits for read to be in .alb file [4]\n"
  "    bin21read <XX>                make .1read file from sorted binary XX.alb and XX.txb files\n"
  "      -T <nthreads>                    number of threads [8]\n"
  "      -o <ZZ.1read>                    output file - default is XX.1read for input XX.{alb,txb}\n"
#ifndef HIDE
  "    bin21bam <YY.1seq> <XX.bin>   make .1bam file from sorted .bin and .1seq\n"
  "      -T <nthreads>                    number of threads [8]\n"
  "      -o <ZZ.1bam>                     output file - default is XX.1bam for input XX.bin\n"
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

  if (!strcmp (command, "bam21bam"))
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
  else if (!strcmp (command, "make1read"))
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
  else if (!strcmp (command, "bin21bam"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-T") && argc > 1) { NTHREAD = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
	else if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bin21bam option %s - run without args for usage", *argv) ;
      if (argc != 2) die ("onebam bin21bam needs 2 not %d args; run without args for usage", argc) ;
      die ("bin21bam not implemented yet") ;
     }
  else if (!strcmp (command, "bin21read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-T") && argc > 1) { NTHREAD = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
	else if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bin21read option %s - run without args for usage", *argv) ;
      if (argc != 2) die ("onebam bin21read needs 2 not %d args; run without args for usage", argc) ;
      die ("bin21read not implemented yet") ;
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

static void accParse (char *acc, char *prefix, I64 *ip, int *digits, char *suffix)
{
  char *s = acc ;
  while (*s && *s < '0' || *s > '9') *prefix++ = *s++ ; *prefix = 0 ;
  *ip = 0 ; *digits = 0 ;
  while (*s >= '0' && *s <= '9') { *ip *= 10 ; *ip += *s++ - '0' ; ++*digits ; }
  while (*s) *suffix++ = *s++ ; *suffix = 0 ;
}

static void buildReference (char *dbname, char *tsvname)
{
  printf ("building reference for %s %s\n", dbname, tsvname) ;

  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *of = oneFileOpenWriteNew (dbname, schema, "ref", true, 1) ;
  if (!of) die ("failed to open reference database ONEfile %s", dbname) ;
  FILE      *in = fzopen (tsvname,"r") ;
  if (!in) die ("failed to open TSV file %s\n", tsvname) ;

  I64 nAcc = 0 ;
  char  accBuf[32], aPref[32], aSuff[32], aPref0[32], aSuff0[255] ; accBuf[31] = 0 ;
  I64 z, z0 = -2 ;
  int d, d0 ;
  while (!feof (in))
    { long long len, taxid ;
      if (fscanf (in, "%31s\t%lld\t%lld\n", accBuf, &len, &taxid) != 3)
	die ("failed to read line %lld", nAcc+1) ;
      oneInt (of,0) = (I64) len ;
      oneInt (of,1) = (I64) taxid ;
      oneWriteLine (of, 'A', 0, 0) ;
      accParse (accBuf, aPref, &z, &d, aSuff) ;
      if (z != ++z0 || d != d0 || strcmp (aPref, aPref0) || strcmp (aSuff, aSuff0))
	{ oneWriteLine (of, 'I', strlen(accBuf), accBuf) ;
	  z0 = z ; d0 = d ; strcpy (aPref0, aPref) ; strcpy (aSuff0, aSuff) ;
	}
      ++nAcc ;
    }
  fclose (in) ;
  oneFileClose (of) ;
  fprintf (stderr, "read %lld accessions: ", (long long)nAcc) ;
  timeUpdate (stderr) ;
}

/********** threaded code to read in the reference ************/

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

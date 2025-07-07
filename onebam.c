/*  File: onebam.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Jul  7 12:04 2025 (rd109)
 * Created: Wed Jul  2 10:18:19 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"

// globals
int        NTHREAD = 8 ;
OneSchema *schema ;

// reference
ExternalReference *reference = 0 ;
static void  buildReference (char *refname, char *tsvname) ;
static bool  readReference (char *refname) ;

// bam
OneFile *ofRead = 0 ;

static char oldUsage[] =
  "Usage: onebam <global_option>* <command> <option|operation>+\n"
  "\n"
  "  global options are:\n"
  "    -T <nthreads>                    number of threads [8]\n"
  "\n"
  "  commands are given below with their options/operations - options/operations are applied in order\n"
  "    reference\n"
  "      -build <XX.1ref> <tsv>         create database [default ending .1ref] from <acc,length, taxid> file\n"
  "    bam\n"
  "      -taxids <XX.tsv>               file with lines <acc>\\ttaxid\\n sorted on acc - required for writeRead\n"
  "      -names                         keep sequence names\n"
  "      -write1read <XX.1read>         write per-read information\n"
  "      -aux <tag>:<fmt> <char>        record BAM tag with OneCode <char>, e.g. -aux AS:i s\n"
  "      -write1bam <XX.1bam>           write out 1bam [default]\n"
  "      -cramref <file_name|URL>       reference for CRAM - needed to read cram - do not confuse with -ref!\n"
  "      -readBam <XX.bam>              read and process a bam (or sam or cram) file\n" ;

static char usage[] =
  "Usage: onebam <option|operation>+\n"
  "  options/operations are applied in order\n"
  "      -T <nthreads>                    number of threads [8]\n"
  "      -taxids <XX.tsv>                 file with lines <acc>\\ttaxid\\n sorted on acc - required for writeRead\n"
  "      -names                           keep sequence names\n"
  "      -aux <tag>:<fmt> <char>          record BAM tag with OneCode <char>, e.g. -aux AS:i s\n"
  "      -cramref <file_name|URL>         reference for CRAM - needed to read cram\n"
  "      -make1bam <XX.1bam> <XX.bam>     make 1bam from BAM\n"
  "      -make1read <XX.1read> <XX.1bam>  make 1read from 1bam\n" ;

int main (int argc, char *argv[])
{
  timeUpdate (0) ;
  storeCommandLine (argc, argv) ;
  argc-- ; ++argv ;
  if (!argc) { fprintf (stderr, "%s", usage) ; exit (0) ; }

  schema = oneSchemaCreateFromText (schemaText) ;

  char *taxidFileName = 0 ;
  bool  isRead = false ;
  bool  isNames = false ;

  while (argc)
    if (!strcmp (*argv, "-T") && argc > 1) { NTHREAD = atoi(argv[1]) ; argc -= 2 ; argv += 2 ; }
    else if (!strcmp (*argv, "-taxids") && argc > 1)
      { taxidFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
    else if (!strcmp (*argv, "-names"))
      { isNames = true ; ++argv ; --argc ; }
    else if (!strcmp (*argv, "-aux") && argc > 2)
      { auxAdd (argv[2], argv[1]) ; argv += 3 ; argc -= 3 ; }
    else if (!strcmp (*argv, "-cramref") && argc > 1)
      { setCramReference (argv[1]) ; argv += 2 ; argc -= 2 ; }
    else if (!strcmp (*argv, "-make1bam") && argc > 2)
      { if (!bamConvert1bam (argv[2], argv[1], taxidFileName, isNames))
	  die ("failed to convert bam file %s", argv[2]) ;
	argv += 3 ; argc -= 3 ;
      }
    else if (!strcmp (*argv, "-make1read") && argc > 2)
      { if (!bamConvert1read (argv[2], argv[1]))
	  die ("failed to convert bam file %s", argv[2]) ;
	argv += 3 ; argc -= 3 ;
      }
    else die ("unknown operation %s for onebam - run without args for usage", *argv) ;

  fprintf (stderr, "Total: ") ; timeTotal (stderr) ;
  exit (0) ;
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

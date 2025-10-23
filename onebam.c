/*  File: onebam.c
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 23 12:52 2025 (rd109)
 * Created: Wed Jul  2 10:18:19 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"

// globals
int        NTHREAD = 8 ;

// reference
static bool makeAccTax (char *refname, char *tsvname) ;

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
  "    report1read <XX.1read>        report on contents of .1read file // needs more work...\n"
  "      -o <ZZ.report>                   output to named file rather than stdout\n"
  "    makeAccTax <XX.tsv>           convert acc2taxid tab-separated file to .1acctax\n"
  "      -o <ZZ.1acctax>                  output file - default is XX.1acctax for input XX.*\n"
  "    addLCA <taxdir> <XX.1read>    add LowestCommonAncestor and dust information\n"
  "      -o <ZZ.1read>                    output file - default is to overwrite input\n"
  "      -T <nthreads>                    number of threads [8]\n"
  //  "      -scoreThresh <T>                 only include taxa with score above T\n"
  //  "      -maxDivergence <M>               maximum divergence, e.g. 0.05 on top of end deamination\n"
  "    reportLCA <XX.1read>          report aggregated by LCA\n"
  "      -o <ZZ.tsv>                      tab-separated output file - default is XX.reportLCA\n"
  "      -rank <taxonomic rank>           collect to and report at named rank, e.g. genus, family\n"
  "      -group <top level tax group>     restrict to group, eg animals, plants, fungi, bacteria\n"
  "    extractReads <XX.1read>       extract reads from XX.1read\n"
  "      -o <ZZ.1read>                    output file - default is XX-<lca>.1read\n"
  "      -lca <lca>                       extract reads with LCA at or below <lca>\n"
#ifndef HIDE
  "    bam21bam <XX.bam>             convert BAM/SAM/CRAM file to .1bam\n"
  "      -o <ZZ.1bam>                     output file - default is XX.1bam for input XX.bam\n"
  "      -accTax <YY.1acctax>             file with lines <acc>\\ttaxid\\n sorted on acc\n"
  "         NB you must give taxids if you plan to use your .1bam to make .1read\n"              
  "      -aux <tag>:<fmt> <char>          record BAM tag with OneCode <char>, e.g. -aux AS:i s\n"
  "         NB you must use '-aux AS:i s -aux MD:Z m' if you plant to use your .1bam to make .1read\n"
  "      -names                           keep sequence names [default to drop names]\n"
  "      -cramref <file_name|URL>         reference for CRAM - needed to read cram\n"
  "    1bam21read <XX.1bam>          convert .1bam file to .1read (simplified information per read)\n"
  // no - this is not sorted so not valid
  "      -T <nthreads>                    number of threads [8]\n"
  "      -o <ZZ.1read>                    output file - default is XX.1read for input XX.1bam\n"
  "    1bam2bam <XX.1bam>            back-convert .1bam to .bam\n"
  "      -o <ZZ.bam>                      output file - default is XX.bam for input XX.1bam\n"
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
  char *accTaxFileName = 0 ;
  bool  isNames = false ;

  if (!strcmp (command, "bam21read"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bam21read option %s - run without args for usage", *argv) ;
      if (argc != 2)
	die ("onebam bam21read needs 2 not %d args; run without args for usage", argc) ;
      if (!bam21readSorted (argv[1], outFileName, *argv))
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
      if (argc != 1)
	die ("onebam report1read needs one not %d args; run without args for usage", argc) ;
      if (!report1read (outFileName, *argv))
	die ("failed to generate report from file %s", *argv) ;
    }
  else if (!strcmp (command, "makeAccTax"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam makeAccTax option %s - run without args for usage", *argv) ;
      if (argc != 1)
	die ("onebam makeAccTax needs one not %d args; run without args for usage", argc) ;
      if (!makeAccTax (outFileName, *argv))
	die ("failed to make .1acctax file from %s", *argv) ;
    }
  else if (!strcmp (command, "addLCA"))
    { int scoreThresh = - (1<<30) ;
      double maxDivergence = 1.0 ;
      int nThreads = 8 ;
      while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-T") && argc > 1)
	  { nThreads = atoi(argv[1]) ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-scoreThresh") && argc > 1)
	  { scoreThresh = atoi(argv[1]) ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-maxDivergence") && argc > 1)
	  { maxDivergence = atof(argv[1]) ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam addLCA option %s - run without args for usage", *argv) ;
      if (argc != 2)
	die ("onebam addLCA needs two not %d args; run without args for usage", argc) ;
      if (!addLCA (outFileName, argv[1], argv[0], scoreThresh, maxDivergence, nThreads))
      	die ("failed to add LCAs to %s using taxonomy in %s", argv[1], argv[0]) ;
    }
  else if (!strcmp (command, "reportLCA"))
    { char *rank = 0, *group = 0 ;
      while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-rank") && argc > 1)
	  { rank = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-group") && argc > 1)
	  { group = argv[1] ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam reportLCA option %s - run without args for usage", *argv) ;
      if (argc != 1)
	die ("onebam reportLCA needs one not %d args; run without args for usage", argc) ;
      if (!reportLCA (argv[0], outFileName, rank, group))
      	die ("failed to report LCAs for %s", *argv) ;
    }
  else if (!strcmp (command, "extractReads"))
    { int lca = 0 ;
      while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-lca") && argc > 1)
	  { lca = atoi(argv[1]) ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam extractReads option %s - run without args for usage", *argv) ;
      if (argc != 1)
	die ("onebam extractReads needs one not %d args; run without args for usage", argc) ;
      if (!extractReads (argv[0], outFileName, lca))
      	die ("failed to extract reads from %s", *argv) ;
    }
  else if (!strcmp (command, "bam21bam"))
    { while (argc && **argv == '-')
	if (!strcmp (*argv, "-o") && argc > 1)
	  { outFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-accTax") && argc > 1)
	  { accTaxFileName = argv[1] ; argv += 2 ; argc -= 2 ; }
	else if (!strcmp (*argv, "-names"))
	  { isNames = true ; ++argv ; --argc ; }
	else if (!strcmp (*argv, "-aux") && argc > 2)
	  { auxAdd (argv[2], argv[1]) ; argv += 3 ; argc -= 3 ; }
	else if (!strcmp (*argv, "-cramref") && argc > 1)
	  { setCramReference (argv[1]) ; argv += 2 ; argc -= 2 ; }
	else die ("unknown onebam bam21bam option %s - run without args for usage", *argv) ;
      if (argc != 1) die ("onebam bam21bam needs 1 not %d args; run without args for usage", argc) ;
      if (!bam21bam (*argv, outFileName, accTaxFileName, isNames))
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
  else
    die ("unknown onebam command %s - run without arguments for usage", command) ;

  fprintf (stderr, "Total: ") ; timeTotal (stderr) ;
  exit (0) ;
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
  if ((c = getc(in)) != EOF && isspace(c))
    { while (isspace(c))
	{ if (c == '\n') { warn ("empty line %lld", nLine) ; return false ; }
	  if ((c = getc(in)) == EOF)
	    { warn ("whitespace on empty final line %lld", nLine) ; return false ; }
	}
      warn ("initial whitespace on line %lld", nLine) ;
    }
  ungetc (c, in) ; // put the initial character back on the stack
  while ((c = getc(in)) != EOF && !isspace(c))
    { *s++ = c ;
      if (++n == 64)
	{ warn ("accession name is too long - line %lld", nLine) ;
	  while ((c = getc(in)) != EOF && c != '\n') { ; } return false ;
	}
      *s = 0 ;
    }
  if (feof (in)) return false ;
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
  destroyCommandLine () ;
  return true ;
}

/*************************** end of file ****************************/

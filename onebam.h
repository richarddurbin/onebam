/*  File: onebam.h
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct  2 00:20 2025 (rd109)
 * Created: Wed Jul  2 10:18:52 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "ONElib.h"

#define VERSION "0.2"

static char *schemaText =
  "1 3 def 1 0                  schema for onebam\n"
  ".\n"
  "P 6 acctax                   ACCESSION TO TAXID INFORMATION\n"
  "O A 5 3 INT 3 INT 3 INT 3 INT 6 STRING  taxid start count shared_prefix_len new_suffix\n"
  ".\n"
  "P 3 seq                      SEQUENCE file - making filetype seq enables .1seq functionality\n"
  "S 4 read                     1read file\n"
  "S 3 bam                      1bam file\n"
  ".\n"
  "O Z 3 4 CHAR 4 CHAR 6 STRING"
  " aux table (for 1bam): ONEcode (lower case), tag, format [ACcSsIifHZ]|[dDtTjJg]"
  " using [dDtTjJg] for B[cCsSiIf] for arrays\n"
  ".\n"
  "O R 3 3 INT 3 INT 6 STRING   reference (for 1bam): length, taxid, accession\n"
  ".\n"
  "O G 1 6 STRING               read group: name - only do this when it changes\n"
  "G S                          groups sequences\n"
  ".\n"
  "O P 1 6 STRING               identifier prefix - for 1read files\n"
  ".\n"
  "O S 1 3 DNA                  the sequence\n"
  "D I 1 6 STRING               identifier\n"
  "D N 3 3 INT 4 CHAR 3 INT     non-acgt base: pos (0-indexed), base, number\n"
  "D Q 1 6 STRING               quality: Q values (ascii string = q+33)\n"
  ". // the following lines for 1read\n"
  "D M 2 3 INT 6 STRING         match string: . for match, ACGT for HiQ sub, acgt for lowQ sub, - for delete\n"
  "D T 3 3 INT 3 INT 3 INT      TaxID, best score(AS), count of records to this taxid\n"
  ". // the following lines for 1bam\n"
  "D B 5 3 INT 3 INT 3 INT 3 INT 8 INT_LIST"
  " BAM line: flags, target, pos, mapq, cigar (stored as BAM I32s) - can have multiple per sequence\n"
  "D X 3 3 INT 3 INT 3 INT      next fragment line for last B line: RNEXT, PNEXT, TLEN\n" ;
// NB the schema must end inside the sequence ("O S 1 DNA") type, so we can add auxiliary linetypes dynamically

// onebam.c
extern int NTHREAD ;
char *derivedName (char *inName, char *tag) ;

// onebamhts.c
void setCramReference (char *cramRef) ;
void auxAdd (char *oneCode, char *spec) ;
bool bam21readSorted (char *bamFileName, char *outFileName, char *accTaxName) ;
bool bam21bam (char *bamFileName, char *outFileName, char *accTaxFileName, bool isNames) ;
bool bamMake1read (char *bamOneFileName, char *outFileName) ;

// oneread.c
// first the user API for reading .1read files

// the OneReader handle gives access to objects in the file - read-only please
typedef struct {
  // properties of the current read, filled by oneReaderNext()
  int   seqLen ;
  char *seq ;
  U8   *qual ;          // 0-based; length seqLen (NB may contain internal 0s - can use str* funcs)
  int   nameLen ;
  char *name ;
  char *namePrefix ;    // initial segment of the name that is shared by all reads in the file
  char *nameEnd ;       // the segment of the current name after the prefix
  int   maxScore ;      // the best score (AS field) of any alignment
  float dustScore ;     // the dust score of this sequence, in range 0 to 100 (max dusty)
  char *mLine ;         // string containing reference substitutions for an alignment with maxScore
                        //     length seqLen ; positions matching seq are '.', deletions '-'
  int   nTax ;          // how many taxids this read has alignments to
  int  *taxid ;         // list of nTax taxids that have hits from this read
  int  *taxCount ;      // list of nTax counts of how many hits to the respective taxid
  int  *taxBestScore ;  // list of nTax best scores for the respective taxid
  int   lca ;           // NCBI taxid of least common ancestor of taxids hit - NOT YET
  
  void *private ;       // a handle for private information for the oneRead code
} OneReader ;

OneReader *oneReaderCreate (char *fileName, int nThreads) ;
// returns an array of nThreads OneReader objects; just a simple pointer if nThreads = 1

void oneReaderDestroy (OneReader *or) ;
// cleans up memory etc.; call on the base pointer (return value of oneReaderCreate)

bool oneReaderNext (OneReader *or) ;
// moves on to the next read, updating the contents of *or; returns false at end of file
// call on or+i (or equivalently &or[i] to access the i'th thread if threaded

bool oneReaderGoto (OneReader *or, U64 i) ;
// goes to the requested read
// goto 0 places before the first read, like when you have just started - oneReaderNext() to read it
// goto 1 loads the first read - oneReaderNext() will read the second read, same for following reads
// call on or+i to access the i'th thread if threaded
// returns true if successful, false if i < 0 or i > number of reads

void oneReaderStats (OneReader *or,
		     U64 *nReads,       // number of reads in the file
		     int *maxSeqLen,    // maximum sequence length
		     U64 *totSeqLen,    // total sequence length
		     int *maxNameLen,   // maximum name length
		     int *maxTax        // max value of nTax ;
		     ) ;
// any of the arguments after the first (or) can be 0 (NULL), in which case they won't be filled
// call on the base pointer if you have a threaded OneReader array

bool   merge1read (char *outfile, int nIn, char **infiles) ;
bool   report1read (char *readFileName, char *outFileName) ;
double dust (const char *seq, int seqLen, int window, int *wCount) ;

// MSDsort.c
void msd_sort (U8 *array, I64 nels, int rsize, int ksize, int depth, int mark, int nthreads) ;

typedef struct {
  I32 prefix, suffix ;  // offsets in fixbuf of prefix and suffix
  I32 count ;
} Accession ;

typedef struct {
  I64        n ;       // number of entries
  I64       *len ;     // n * length
  I32       *taxid ;   // n * taxid
  I32       *k ;       // n * index into original accession order
  I64        nAcc, fixBufSize ;
  Accession *acc ;
  char      *fixBuf ;  // holds the prefixes and suffixes
} ExternalReference ;
extern ExternalReference *reference ;

/*********************** end of file *************************/

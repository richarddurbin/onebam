/*  File: onebam.h
 *  Author: Richard Durbin (rd109@sanger.ac.uk)
 *  Copyright (C) Richard Durbin, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  4 23:33 2025 (rd109)
 * Created: Wed Jul  2 10:18:52 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "ONElib.h"

static char *schemaText =
  "1 3 def 1 0                  schema for onebam\n"
  ".\n"
  "P 3 ref                      REFERENCE DATABASE\n"
  "O A 2 3 INT 3 INT            accession: length, taxid\n"
  "D I 1 6 STRING               accession identifier\n"
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
  "O S 1 3 DNA                  the sequence\n"
  "D N 3 3 INT 4 CHAR 3 INT     non-acgt base: pos (0-indexed), base, number\n"
  "D Q 1 6 STRING               quality: Q values (ascii string = q+33)\n"
  ". // the following lines for 1read\n"
  "D T 3 3 INT 3 INT 3 INT      TaxID, best score(AS), count of records to this taxid\n"
  "D M 2 3 INT 6 STRING         match string: . for match, ACGT for HiQ sub, acgt for lowQ sub, - for delete\n"
  ". // the following lines for 1bam\n"
  "D J 2 3 INT 6 STRING         identifier: updates previous with <STRING> from character <INT>\n"
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
bool bam21bam (char *bamFileName, char *outFileName, char *taxidFileName, bool isNames) ;
bool bamMake1read (char *bamOneFileName, char *outFileName) ;
bool makeBin (char *bamFileName, char *outTxbName, char *outAlbName, char *taxidFileName,
	      int maxEdit, int prefix, int maxChars) ;

// albcode.c
void albReport (char *fileName, char *outFileName) ;

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

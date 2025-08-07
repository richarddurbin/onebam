/*  File: txbview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 08:34 2025 (rd109)
 * Created: Thu Aug  7 08:02:54 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;

typedef struct {
  int txid ;
  int count ;
  int bestScore ;
} TaxInfo ; // info per taxid for this thread

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (!argc) die ("Usage: txbview XXX.txb [start [count]]") ;
  FILE *f = fopen (*argv++, "r") ; --argc ;
  if (!f) die ("failed to open .alb file %s", *--argv) ;
  int start = 0, count = 20 ;
  if (argc--) start = atoi (*argv++) ;
  if (argc--) count = atoi (*argv++) ;
  U8 b ;
  if ((b = fgetc(f)) != 0) die ("first byte %d of input is non-zero", b) ;
  int prefixLen = (U8)fgetc(f), maxChars = (U8)fgetc(f) ;
  int recordLen = maxChars + (int)sizeof(TaxInfo) ;
  char *prefix = malloc (prefixLen+1) ;
  int z = (prefixLen <= recordLen - 3) ? prefixLen : recordLen - 3 ; // read 3 fgetc's
  if (fread (prefix, z, 1, f) != 1) die ("failed to read prefix") ;
  prefix[z] = 0 ;
  int i ; for (i = 0 ; i < recordLen - 3 - z ; ++i) fgetc (f) ; // ignore rest of record
  printf ("recordLen %d prefixLen %d maxChars %d prefix %s\n",
	  recordLen, prefixLen, maxChars, prefix) ;
  char *buf = malloc(recordLen) ;
  U64 nRecord = 0 ;
  while ((fread(buf,recordLen,1,f) == 1))
    { ++nRecord ;
      if (nRecord >= start && nRecord < start + count)
	{ TaxInfo *tx = (TaxInfo*)(buf+maxChars) ;
	  printf ("%-*.*s ", prefixLen, prefixLen, buf) ;
	  printf ("%6d %4d %4d\n", tx->txid, tx->count, tx->bestScore) ;
	}
    }
  if (!feof(f)) die ("failed to read after record %llu while not at end of file", nRecord) ;
  printf ("read %llu records\n", nRecord) ;
  fclose (f) ;
}

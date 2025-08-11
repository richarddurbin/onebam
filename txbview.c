/*  File: txbview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 11 12:55 2025 (rd109)
 * Created: Thu Aug  7 08:02:54 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;

typedef struct {
  U32 txid ;
  U16 count ;
  I16 bestScore ;
} TaxInfo ; // info per taxid for this thread

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (!argc) die ("Usage: txbview XXX.txb [start [count]]") ;
  FILE *f = fopen (*argv++, "rb") ; --argc ;
  if (!f) die ("failed to open .txb file %s", *--argv) ;
  int start = 0, count = 10 ;
  if (argc-- > 0) start = atoi (*argv++) ;
  if (argc-- > 0) count = atoi (*argv++) ;
  U8 b ;
  if ((b = fgetc(f)) != 0) die ("first byte %d of input is non-zero", b) ;
  int prefixLen = (U8)fgetc(f), maxChars = (U8)fgetc(f) ;
  int recordLen = maxChars + (int)sizeof(I32) ;
  char *prefix = calloc (prefixLen+1,1) ;
  if (fread (prefix, prefixLen, 1, f) != 1) die ("failed to read prefix") ;
  int i ; for (i = 0 ; i < recordLen - 3 - prefixLen ; ++i) fgetc (f) ; // ignore rest of record
  printf ("recordLen %d prefixLen %d maxChars %d prefix %s\n",
	  recordLen, prefixLen, maxChars, prefix) ;
  char *buf = malloc(recordLen) ;
  U64 nRecord = 0, nHits = 0 ;
  TaxInfo tx ;
  TaxInfo *txBuf = new(1<<20,TaxInfo) ;
  while ((fread(buf,recordLen,1,f) == 1))
    { ++nRecord ;
      int n = *(I32*)(buf+maxChars) ;
      nHits += n ;
      if (nRecord >= start && nRecord < start + count)
	{ printf ("%.*s %d\n", prefixLen, buf, n) ;
	  while (n--)
	    { if ((fread(&tx,sizeof(TaxInfo),1,f) != 1)) die ("failed hit read") ;
	      printf ("  %6u %4u %4d\n", tx.txid, tx.count, tx.bestScore) ;
	    }
	}
      else
	{ if (n >= 1<<20) die ("too many hits %d for record %llu", n, nRecord) ;
	  if ((fread(txBuf, sizeof(TaxInfo), n, f) != n)) die ("failed hit block read") ;
	}
    }
  if (!feof(f)) die ("failed to read after record %llu while not at end of file", nRecord) ;
  printf ("read %llu records with %llu hits\n", nRecord, nHits) ;
  fclose (f) ;
}

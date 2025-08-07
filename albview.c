/*  File: albview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: simple utility for checking .alb files
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  7 08:35 2025 (rd109)
 * Created: Mon Aug  4 13:43:40 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

static const char binary2char[] = "-ACMGRSVTWYHKDBN" ;

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (!argc) die ("Usage: albview XXX.alb") ;
  FILE *f = fopen (*argv, "r") ;
  if (!f) die ("failed to open .alb file %s", *argv) ;
  U8 b ;
  if ((b = fgetc(f)) != 0) die ("first byte %d of input is non-zero", b) ;
  int prefixLen = (U8)fgetc(f), maxChars = (U8)fgetc(f), maxEdit = (U8)fgetc(f) ;
  int recordLen = maxChars + 2*maxEdit + 13 ;
  char *prefix = malloc (prefixLen+1) ;
  int z = (prefixLen <= recordLen - 4) ? prefixLen : recordLen - 4 ;
  if (fread (prefix, z, 1, f) != 1) die ("failed to read prefix") ;
  prefix[z] = 0 ;
  int i ; for (i = 0 ; i < recordLen - 4 - z ; ++i) fgetc (f) ; // ignore rest of record
  printf ("recordLen %d prefixLen %d maxChars %d maxEdit %d prefix %s\n",
	  recordLen, prefixLen, maxChars, maxEdit, prefix) ;
  char *buf = malloc(recordLen) ;
  U64 nRecord = 0 ;
  while ((fread(buf,recordLen,1,f) == 1))
    { ++nRecord ;
      if (nRecord < 10)
	{ int seqLen = ((U8*)buf)[maxChars] ;
	  int score = (int)*(I32*)(buf+maxChars+1+2*maxEdit) ;
	  printf ("%-*.*s len %3d score %3d", prefixLen, prefixLen, buf, seqLen, score) ;
	  for (i = 0 ; i < maxEdit ; ++i)
	    { U8 pos = (U8)buf[prefixLen+1+2*i] ;
	      U8 x = (U8)buf[prefixLen+1+2*i+1] ;
	      if (!pos && !x) break ;
	      printf (" %d:%c>%c", pos, binary2char[x >> 4], binary2char[x & 0xf]) ;
	    }
	  putchar ('\n') ;
	}
    }
  if (!feof(f)) die ("failed to read after record %llu while not at end of file", nRecord) ;
  printf ("read %llu records\n", nRecord) ;
  fclose (f) ;
}

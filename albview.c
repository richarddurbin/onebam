/*  File: albview.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  4 15:39 2025 (rd109)
 * Created: Mon Aug  4 13:43:40 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  if (!argc) die ("Usage: albview XXX.alb") ;
  FILE *f = fopen (*argv, "r") ;
  if (!f) die ("failed to open .alb file %s", *argv) ;
  U8 b ;
  if ((b = fgetc(f)) != 0) die ("first byte %d of input is non-zero", b) ;
  int prefixLen = (U8)fgetc(f), maxChars = (U8)fgetc(f), maxEdit = (U8)fgetc(f) ;
  int recordLen = maxChars + 2*maxEdit + 4 ;
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
      
    }
  if (!feof(f)) die ("failed to read after record %llu while not at end of file", nRecord) ;
  printf ("read %llu records\n", nRecord) ;
  fclose (f) ;
}

/*  File: albcode.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Aug  8 23:50 2025 (rd109)
 * Created: Mon Aug  4 19:34:03 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"

static U8   bin2index[] = { -1, 0, 1, -1, 2, -1, -1, -1, 3, -1, -1, -1, -1, -1, -1, -1 } ;
static char index2char[] = "ACGT" ;

static void readHeaderLine (FILE *f, int *prefixLen, int *maxChars, int *maxEdit, int *recordLen,
			    char **prefix)
{
  *prefixLen = (U8)fgetc(f) ; // second byte - the first byte was 0
  *maxChars = (U8)fgetc(f) ;
  *maxEdit = (U8)fgetc(f) ;
  *recordLen = *maxChars + 2 * *maxEdit + 13 ;
  *prefix = malloc (*prefixLen+1) ;
  int z = (*prefixLen <= *recordLen - 4) ? *prefixLen : *recordLen - 4 ;
  if (fread (*prefix, z, 1, f) != 1) die ("failed to read prefix") ;
  (*prefix)[z] = 0 ;
  int i ; for (i = 0 ; i < *recordLen - 4 - z ; ++i) fgetc (f) ; // ignore rest of record
}

void albReport (char *albFileName, char *outFileName)
{
  int i, j, k ;
  
  // open the output file
  FILE *fOut ;
  if (outFileName)
    { fOut = fopen (outFileName, "w") ;
      if (!fOut) die ("failed to open output file %s", outFileName) ;
    }
  else
    fOut = stdout ;

  // open the alb file
  FILE *fAlb = fopen (albFileName, "r") ;
  if (!fAlb) die ("failed to open .alb file %s", albFileName) ;
  U8 b ; if ((b = fgetc(fAlb)) != 0) die ("first byte %d of input is non-zero", b) ;

  // read the header line
  int prefixLen, maxChars, maxEdit, recordLen ;
  char *prefix ;
  readHeaderLine (fAlb, &prefixLen, &maxChars, &maxEdit, &recordLen, &prefix) ;
  fprintf (fOut, "ALB %s recordLen %d prefixLen %d maxChars %d maxEdit %d prefix %s\n",
	   albFileName, recordLen, prefixLen, maxChars, maxEdit, prefix) ;

  // set up accumulation arrays - initialise to 0 by setting first element 0
  U64  lengths[256] = {0} ;
  U64  start[15][4][4] = {0}, end[15][4][4] = {0}, mid[4][4] = {0} ;
  U64 *edits = new0 (maxEdit+1,U64) ;
  U64 *scores = new0 (256, U64) ;
  I32  minScore = I32MAX, maxScore = -I32MAX ;
  
  // process records
  char *buf = malloc(recordLen) ;
  U64 nRecord = 0 ;
  while ((fread(buf,recordLen,1,fAlb) == 1))
    { ++nRecord ;
      int seqLen = ((U8*)buf)[maxChars] ;
      ++lengths[seqLen] ;
      int score = (int)*(I32*)(buf+maxChars+1+2*maxEdit) ;
      if (score > maxScore) maxScore = score ;
      if (score < minScore) minScore = score ;
      // if (score > 127 || score < -128) warn ("score %d record %llu", score, nRecord) ;
      else ++scores[score+128] ;
      // if (scores[1]) die ("something wrote to scores[1] at record %llu", nRecord) ;
      for (i = 0 ; i < maxEdit ; ++i)
	{ U8 pos = (U8)buf[prefixLen+1+2*i] ;
	  if (!pos) break ; else --pos ;  // because pos was 1-based
	  U8 x = (U8)buf[prefixLen+1+2*i+1] ;
	  U8 to = bin2index[x >> 4], from = bin2index[x & 0xf] ;
	  if (from >= 0 && to >= 0)
	    { if (pos < 15) ++start[pos][from][to] ;
	      else if (seqLen - 1 - pos < 15) ++end[seqLen-1-pos][from][to] ;
	      else ++mid[from][to] ;
	    }
	}
      ++edits[i] ;
    }
  if (!feof(fAlb)) die ("failed to read after record %llu while not at end of file", nRecord) ;
  fclose (fAlb) ;

  // report
  fprintf (fOut, "RECORDS %llu\n", nRecord) ;
  if (nRecord)
    { int min = 0, max = 0, n50 ;
      U64 sum = 0, psum = 0 ;
      for (i = 0 ; i < 256 ; ++i)
	{ sum += i*lengths[i] ;
	  if (lengths[i] && !min) min = i ;
	  if (lengths[i]) max = i ;
	}
      for (i = 0 ; i < 256 ; ++i) { psum += i*lengths[i] ; if (psum*2 < sum) n50 = i+1 ; }
      double mean = sum/(double)nRecord ;
      fprintf (fOut, "LENGTHS min %d max %d mean %.1f n50 %d\n", min, max, mean, n50) ;
      for (i = 255 ; i >= 0 ; --i)
	if (scores[i]) fprintf (fOut, "SCORES %4d %llu\n", i-128, scores[i]) ;
      fprintf (fOut, "SCORE_MIN_MAX %d %d\n", minScore, maxScore) ;
      fprintf (fOut, "EDITS") ;
      for (i = 0 ; i <= maxEdit ; ++i) fprintf (fOut, " %llu", edits[i]) ;
      fprintf (fOut, "\n") ;
      double ct, ga, other ;
      for (i = 0 ; i < 15 ; ++i)
	{ ct = ga = other = 0. ;
	  double nC = start[i][1][0] + start[i][1][1] + start[i][1][2] + start[i][1][3] ;
	  if (nC) ct = start[i][1][3] / nC ;
	  double nG = start[i][2][0] + start[i][2][1] + start[i][2][2] + start[i][2][3] ;
	  if (nG) ga = start[i][2][0] / nG ;
	  double nX = start[i][0][0] + start[i][0][1] + start[i][0][2] + start[i][0][3] +
	    nC + nG + start[i][3][0] + start[i][3][1] + start[i][3][2] + start[i][3][3]
	    - start[i][1][3] - start[i][2][0] ;
	  if (nX) other = (nX-start[i][0][0]-start[i][1][1]-start[i][2][2]-start[i][3][3]) / nX ;
	  fprintf (fOut, "START %2d C>T %.3f G>A %.3f OTHER %.3f\n", i+1, ct, ga, other) ;
	}
      { ct = ga = other = 0. ;
	double nC = mid[1][0] + mid[1][1] + mid[1][2] + mid[1][3] ;
	if (nC) ct = mid[1][3] / nC ;
	double nG = mid[2][0] + mid[2][1] + mid[2][2] + mid[2][3] ;
	if (nG) ga = mid[2][0] / nG ;
	double nX = mid[0][0] + mid[0][1] + mid[0][2] + mid[0][3] +
	  nC + nG + mid[3][0] + mid[3][1] + mid[3][2] + mid[3][3]
	  - mid[1][3] - mid[2][0] ;
	if (nX) other = (nX-mid[0][0]-mid[1][1]-mid[2][2]-mid[3][3]) / nX ;
	fprintf (fOut, "MID      C>T %.3f G>A %.3f OTHER %.3f\n", ct, ga, other) ;
      }
      for (i = 0 ; i < 15 ; ++i)
	{ ct = ga = other = 0. ;
	  double nC = end[i][1][0] + end[i][1][1] + end[i][1][2] + end[i][1][3] ;
	  if (nC) ct = end[i][1][3] / nC ;
	  double nG = end[i][2][0] + end[i][2][1] + end[i][2][2] + end[i][2][3] ;
	  if (nG) ga = end[i][2][0] / nG ;
	  double nX = end[i][0][0] + end[i][0][1] + end[i][0][2] + end[i][0][3] +
	    nC + nG + end[i][3][0] + end[i][3][1] + end[i][3][2] + end[i][3][3]
	    - end[i][1][3] - end[i][2][0] ;
	  if (nX) other = (nX-end[i][0][0]-end[i][1][1]-end[i][2][2]-end[i][3][3]) / nX ;
	  fprintf (fOut, "END   %2d C>T %.3f G>A %.3f OTHER %.3f\n", i+1, ct, ga, other) ;
	}
      for (i = 0 ; i < 4 ; ++i)
	for (j = 0 ; j < 4 ; ++j)
	  { fprintf (fOut, " COUNTS %c>%c", index2char[i], index2char[j]) ;
	    fprintf (fOut, " MID %llu START", mid[i][j]) ;
	    for (k = 0 ; k < 15 ; ++k) fprintf (fOut, " %llu", start[k][i][j]) ;
	    fprintf (fOut, " END") ;
	    for (k = 15 ; k-- ; ) fprintf (fOut, " %llu", end[k][i][j]) ;
	    fprintf (fOut, "\n") ;
	  }
    }
  
  if (fOut) fclose (fOut) ;
  newFree (edits,maxEdit+1,U64) ;
}

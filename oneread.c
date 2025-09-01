/*  File: oneread.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: onebam functionality involving .1read files
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 31 23:39 2025 (rd109)
 * Created: Wed Aug 13 14:10:49 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"
#include "merge.h"

typedef struct {
  char    *fileName ;
  OneFile *of ;
  I64      seqLen ;
  U8      *dna ;      // keep pointer to 2bit form - requires less work
  int      lcp ;
  char     name[1024] ;
} Input ;

bool yieldName (int t, void *arg, char **v, int *lcp) // here arg is an array of Input
{
  Input *in = &(((Input*)arg)[t]) ;
  if (in->of) { *lcp = in->lcp ; *v = in->name ; return true ; }
  else return false ;
}

bool yieldTx (int t, void *arg, I32 *v) // here arg is an array of OneFile
{
  OneFile *of = (((OneFile**)arg)[t]) ;
  if (of) { *v = oneInt(of,0) ; return true ; }
  else return false ;
}

void loadSequence (Input *in)
{
  if (in->of->lineType != 'S')
    die ("line type %c instead of expected S in %s line %lld",
	 in->of->lineType, in->fileName, in->of->line) ;
  in->seqLen = oneLen(in->of) ;
  in->dna = oneDNA2bit(in->of) ;
  while (oneReadLine (in->of) && (in->of->lineType != 'J')) ; // really should be next line...
  if (!in->of->lineType) die ("no J line found in %s", in->fileName) ;
  in->lcp = oneInt(in->of, 0) ;
  if (in->lcp + oneLen(in->of) > 1024)
    die ("read name length %d > 1024 file %s line %d",
	 in->lcp + oneLen(in->of), in->fileName, in->of->line) ;
  strcpy (in->name + in->lcp, oneString(in->of)) ;
}

bool merge1read (char *outfile, int nIn, char **infiles)
{
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;

  if (!outfile) outfile = "merged.1read" ;
  OneFile *ofOut = oneFileOpenWriteNew (outfile, schema, "read", true, 1) ;
  if (!ofOut) { warn ("failed to open %s to write", outfile) ; return false ; }
  oneAddProvenance (ofOut, "onebam", VERSION, getCommandLine()) ;

  int  i ;
  char t ;
  U64  nSin = 0, nTin = 0 ; // counters for input S and T lines

  // next build the inputs array and create the read merge object
  Input *inputs = new0 (nIn, Input) ;
  for (i = 0 ; i < nIn ; ++i)
    { Input *in = &inputs[i] ;
      in->fileName = infiles[i] ;
      if (!(in->of = oneFileOpenRead (in->fileName, 0, "read", 1)))
	die ("failed to open .1read file %s", in->fileName) ;
      if (!oneFileCheckSchemaText (in->of, "P 3 seq\nO S 1 3 DNA\nD J 2 3 INT 6 STRING\n"
				   "D M 2 3 INT 6 STRING\nD T 3 3 INT 3 INT 3 INT\n"))
	die ("schema mismatch for input file %s", in->fileName) ;
      if (oneGoto (in->of, 'S', 1) && oneReadLine(in->of)) // oneReadLine required to read S line
	{ loadSequence (in) ; ++nSin ; }
      else
	{ warn ("input .1read file %s appears empty", in->fileName) ; // it could actually be empty
	  in->of = 0 ;
	}
    }
  Merge   *readMerge = mergeCreateString (nIn, inputs, yieldName) ;

  // do this now because inside the loop we will mergeRecreate(txMerge)
  OneFile **txOf      = new0 (nIn, OneFile*) ;
  Merge    *txMerge   = mergeCreateI32 (nIn, txOf, yieldTx) ;

  // now start the merge process
  int    nName, *nameList, nTx, *txList ;
  char   lastOutName[1024] ; *lastOutName = 0 ;
  while ((nName = mergeNext (readMerge, &nameList)))
    { Input *in = &inputs[nameList[0]] ; // this is an example to copy the shared information from

      // first write the sequence record 'S'
      oneWriteLineDNA2bit (ofOut, 'S', in->seqLen, in->dna) ; // write the sequence record

      // and the name in lcp/suffix form 'J'
      char *p = lastOutName, *q = in->name ; while (*p && *p == *q) { ++p ; ++q ; }
      oneInt(ofOut,0) = (I64)(q - in->name) ;
      oneWriteLine (ofOut, 'J', strlen(q), q) ;
      strcpy (p, q) ; // copy the distinct suffix onto the shared prefix to update lastOutName

      if (nName == 1) // all very simple - go to next S line
	{ while ((t = oneReadLine(in->of)) && t != 'S')
	    { oneWriteLineFrom (ofOut, in->of) ; // copy line
	      if (t == 'T') ++nTin ;
	    }
	  if (t) { loadSequence (in) ; ++nSin ; }
	  else in->of = 0 ; // termination index for this input for readMerge
	}
      else // multiple input streams for the same read
	{ // first copy all lines up to M, which should be the same across all inputs
	  while ((t = oneReadLine(in->of)) && t != 'M')  // there should be an M line for each object
	    oneWriteLineFrom (ofOut, in->of) ; // copy line
	  if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[0]].fileName) ;
	  for (i = 1 ; i < nName ; ++i)
	    { txOf[i] = inputs[nameList[i]].of ; // will need this later, and it is more direct
	      while ((t = oneReadLine(txOf[i])) && t != 'M') ; // read these lines also
	      if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[i]].fileName) ;
	    }
	  txOf[0] = in->of ;
	  
	  // now find the highest scoring M line and write it out
	  I64 mMax = oneInt(in->of,0) ;
	  OneFile *ofM = in->of ;
	  for (i = 1 ; i < nName ; ++i)
	    if (oneInt(txOf[i],0) > mMax) { ofM = txOf[i] ; mMax = oneInt(ofM,0) ; }
	  oneWriteLineFrom (ofOut, ofM) ; // copy the highest scoring M line
	  
	  // move each txOf[] to the next line, which should be the first T line, and set up txMerge
	  for (i = 0 ; i < nName ; ++i)
	    if (oneReadLine(txOf[i]) == 'T') ++nTin ;
	    else die ("missing T line after M line input %s line %lld",
		      inputs[nameList[i]].fileName, (long long)txOf[i]->line) ;
	  mergeRecreate (txMerge, nName, txOf) ;

	  // and merge the taxids
	  while ((nTx = mergeNext (txMerge, &txList)))
	    { OneFile *ofTx = txOf[txList[0]] ;
	      oneInt(ofOut,0) = oneInt(ofTx,0) ; // copy the taxid
	      oneInt(ofOut,1) = oneInt(ofTx,1) ; // take the max of the best scores
	      oneInt(ofOut,2) = oneInt(ofTx,2) ; // sum the counts
	      for (i = 1 ; i < nTx ; ++i)
		{ ofTx = txOf[txList[i]] ;
		  if (oneInt(ofTx,1) > oneInt(ofOut,1)) oneInt(ofOut,1) = oneInt(ofTx,1) ;
		  oneInt(ofOut,2) += oneInt(ofTx,2) ;
		}
	      oneWriteLine (ofOut, 'T', 0, 0) ;
	      // now advance each line, and set txOf[] to 0 if not 'T'
	      for (i = 0 ; i < nTx ; ++i)
		if (oneReadLine (txOf[txList[i]]) == 'T') ++nTin ;
		else txOf[txList[i]] = 0 ; // terminator for the txMerge
	    }

	  // finally update the sequence and name in the inputs[] for the nameList 
	  for (i = 0 ; i < nName ; ++i)
	    { Input *in = &inputs[nameList[i]] ;
	      if (in->of->lineType) { loadSequence (in) ; ++nSin ; }
	      else { oneFileClose (in->of) ; in->of = 0 ; } // terminator for the nameMerge
	    }
	}
    }

  printf ("merged %d .1read files with %lld S lines and %lld T lines",
	  nIn, (long long) nSin, (long long) nTin) ;
  I64 nSout ; oneStats (ofOut, 'S', &nSout, 0, 0) ; 
  I64 nTout ; oneStats (ofOut, 'T', &nTout, 0, 0) ; 
  printf (" into %s with %lld S lines and %lld T lines\n",
	  oneFileName(ofOut), (long long)nSout, (long long)nTout) ;

  oneFileClose (ofOut) ;

  mergeDestroy (readMerge) ;
  mergeDestroy (txMerge) ;

  newFree (inputs, nIn, Input) ;
  newFree (txOf, nIn, OneFile) ;

  return true ;
}

/*********************************************************/

bool report1read (char *readFileName, char *outFileName)
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

  // open the 1read file
  OneSchema *schema = oneSchemaCreateFromText (schemaText) ;
  OneFile   *of = oneFileOpenRead (readFileName, 0, "read", 1) ;
  if (!of) die ("failed to open .1read file %s", readFileName) ;

  // set up accumulation arrays - initialise to 0 by setting first element 0
  U64  lengths[256] = {0} ;
  U64  start[15][4][4] = {0}, end[15][4][4] = {0}, mid[4][4] = {0} ;
  U64 *edits = new0 (256, U64) ;
  U64 *scores = new0 (256, U64) ;
  I32  minScore = I32MAX, maxScore = -I32MAX ;
  int  maxEdit = 0 ;

  int  map[256] ; for (i = 0 ; i < 256 ; ++i) map[i] = -1 ;
  map['A']=map['a']=0 ; map['C']=map['c']=1 ; map['G']=map['g']=2 ; map['T']=map['t']=3 ;
  char index2char[] = "acgt" ;
  
  // process records
  U64 nRecord = 0 ;
  if (!oneGoto (of, 'S', 1)) die ("no objects in 1read file %s", readFileName) ;
  oneReadLine (of) ; // read the 'S' line
  while (of->lineType == 'S')
    { ++nRecord ;
      int seqLen = oneLen(of) ;
      ++lengths[seqLen] ;
      char *seq = oneDNAchar(of) ;

      // process the M line
      while (oneReadLine (of) && of->lineType != 'M') if (of->lineType == 'S') continue ;
      int score = oneInt(of,0) ;
      if (score > maxScore) maxScore = score ;
      if (score < minScore) minScore = score ;
      if (score > 127 || score < -128) warn ("score %d record %llu", score, nRecord) ;
      else ++scores[score+128] ;
      char *m = oneString(of) ;
      int nEdit = 0 ;
      for (i = 0 ; i < seqLen ; ++i)
	if (m[i] != '.') // an edit
	  { ++nEdit ;
	    int from = map[m[i]], to = map[seq[i]] ;
	    if (from >= 0 && to >= 0)
	      { if (i < 15) ++start[i][from][to] ;
		else if (seqLen - 1 - i < 15) ++end[seqLen-1-i][from][to] ;
		else ++mid[from][to] ;
	      }
	  }
      ++edits[nEdit] ;
      if (nEdit > maxEdit) maxEdit = nEdit ;

      // for now skip all the T lines - should process them
      while (oneReadLine (of) && of->lineType != 'S') ;
    }
  oneFileClose (of) ;

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
  newFree (edits,256,U64) ;
  newFree (scores,256,U64) ;
  return true ;
}

/********************* end of file ***********************/

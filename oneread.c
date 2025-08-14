/*  File: oneread.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: onebam functionality involving .1read files
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 14 01:46 2025 (rd109)
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
  OneFile *of = &(((OneFile*)arg)[t]) ;
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
  OneFile *ofOut = oneFileOpenWriteNew (outfile, schema, "read", true, NTHREAD) ;
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
	    transferLine (ofOut, in->of) ; // copy line
	  if (t) { loadSequence (in) ; ++nSin ; }
	  else in->of = 0 ; // termination index for this input for readMerge
	}
      else
	{ // first copy all lines up to M, which should be the same across all inputs
	  while ((t = oneReadLine(in->of)) && t != 'M')  // there should be an M line for each object
	    transferLine (ofOut, in->of) ; // copy line
	  if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[0]].fileName) ;
	  for (i = 1 ; i < nName ; ++i)
	    { txOf[i] = inputs[nameList[i]].of ; // will need this later, and it is more direct
	      while ((t = oneReadLine(txOf[i]) && t != 'M')) ; // read these lines also
	      if (t != 'M') die ("failed to find M line on input %s", inputs[nameList[i]].fileName) ;
	    }
	  txOf[0] = in->of ;
	  
	  // now find the highest scoring M line and write it out
	  I64 mMax = oneInt(in->of,0) ;
	  OneFile *ofM = in->of ;
	  for (i = 1 ; i < nName ; ++i)
	    if (oneInt(txOf[i],0) > mMax) { ofM = txOf[i] ; mMax = oneInt(ofM,0) ; }
	  transferLine (ofOut, ofM) ; // copy the highest scoring M line
	  
	  // move each txOf[] to the next line, which should be the first T line, and set up txMerge
	  for (i = 0 ; i < nName ; ++i)
	    if (oneReadLine(txOf[i]) == 'T') ++nTin ;
	    else die ("missing T line after M line input %s line %lld",
		      inputs[nameList[i]].fileName, (long long)txOf[i]->line) ;
	  mergeRecreate (txMerge, nName, txOf) ;

	  // and merge the taxids
	  while ((nTx = mergeNext (txMerge, &txList)))
	    { oneInt(ofOut,0) = oneInt(txOf[0],0) ; // copy the taxid
	      oneInt(ofOut,1) = oneInt(txOf[0],1) ; // take the max of the best scores
	      oneInt(ofOut,2) = oneInt(txOf[0],2) ; // sum the counts
	      for (i = 1 ; i < nTx ; ++i)
		{ if (oneInt(txOf[txList[i]],1) > oneInt(ofOut,1))
		    oneInt(ofOut,1) = oneInt(txOf[txList[i]],1) ;
		  oneInt(ofOut,2) += oneInt(txOf[txList[i]],2) ;
		}
	      oneWriteLine (ofOut, 'T', 0, 0) ;
	      // now advance each line, and set txOf[] to 0 if not 'T'
	      for (i = 0 ; i < nTx ; ++i)
		if (oneReadLine (txOf[i]) == 'T') ++nTin ;
		else txOf[i] = 0 ; // terminator for the txMerge
	    }

	  // finally update the sequence and name in the inputs[] for the nameList 
	  for (i = 0 ; i < nName ; ++i)
	    { Input *in = &inputs[nameList[i]] ;
	      if (in->of->lineType)	{ loadSequence (in) ; ++nSin ; }
	      else in->of = 0 ; // terminator for the nameMerge
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
  for (i = 0 ; i < nIn ; ++i) oneFileClose (inputs[i].of) ;

  newFree (inputs, nIn, Input) ;
  newFree (txOf, nIn, OneFile) ;

  return true ;
}

/********************* end of file ***********************/

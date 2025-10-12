/*  File: onebamtax.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 12 22:59 2025 (rd109)
 * * Oct  7 17:07 2025 (rd109): 
 * Created: Tue Oct  7 17:06:32 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "onebam.h"

typedef struct {
  OneReader *or ;
  U64        start, end ;
  TaxLCA    *txLCA ;		// shared - readonly
  bool      *txUsed ;		// private
} ThreadArg ;

static void *threadProcess (void *arg)
{
  ThreadArg *ta = (ThreadArg*)arg ;
  if (!oneReaderGoto (ta->or, ta->start)) die ("failed oneReaderGoto") ;
  int n = 0 ;
  while (ta->start++ < ta->end && oneReaderNext (ta->or))
    { ta->or->lca = taxLCAfind (ta->txLCA, ta->or->nTax, ta->or->taxid) ;
      oneReaderWrite (ta->or) ;
      int i ; for (i = 0 ; i < ta->or->nTax ; ++i) ta->txUsed[ta->or->taxid[i]] = true ;
    }
  return 0 ;
}

#include <sys/errno.h>
#include <unistd.h>
#include <sys/stat.h> // for umask() and fchmod()

bool addLCA (char *outFileName, char *inFileName, char *taxPath,
	     int scoreThresh, double maxDivergence, int nThreads)
{
  OneReader *or = oneReaderCreate (inFileName, nThreads) ;
  if (!or) die ("failed to open .1read file %s", inFileName) ;

  Taxonomy *tx = taxonomyFromNCBIfiles (taxPath) ;
  if (!tx) die ("failed to create taxonomy from directory %s", taxPath) ;
  printf ("read taxonomy from %s with %d entries: ", taxPath, tx->nTx) ; timeUpdate (stdout) ;
  
  char outName[] = "./onebam-LCA-XXXXXX" ;
  if (outFileName)
    { if (!oneReaderWriteFile (or, outFileName)) die ("failed to open %s", outFileName) ; }
  else
    { int fd = mkstemp (outName) ;
      if (errno) die ("failed to open temporary file %s", outName) ;
      // Set permissions respecting umask - code from Claude
      mode_t mask = umask (0) ;   // Get current umask
      umask (mask) ;              // Restore it immediately
      fchmod (fd, 0666 & ~mask) ; // Apply desired perms with umask
      close (fd) ; // just did this to get unique name and confirm we could open write to it
      oneReaderWriteFile (or, outName) ;
    }

  pthread_t *threads = new (nThreads, pthread_t) ;
  ThreadArg *ta      = new0 (nThreads, ThreadArg) ;
  U64 nRead ; oneReaderStats (or, &nRead, 0, 0, 0, 0) ;
  if (!nRead) die ("no reads found in file %s", inFileName) ;
  int i, j ;
  for (i = 0 ; i < nThreads ; ++i)
    { ta[i].or     = or + i ;                     // equivalent to ta[i].or = &or[i]
      ta[i].start  = (nRead * i) / nThreads ;
      ta[i].end    = (nRead * (i+1)) / nThreads ;
      ta[i].txLCA  = taxLCAcreate (tx) ;
      ta[i].txUsed = new0 (arrayMax(tx->nodes), bool) ;
    }
  for (i = 0 ; i < nThreads ; ++i)              // create threads
    pthread_create (&threads[i], 0, threadProcess, &ta[i]) ;
  for (i = 0 ; i < nThreads ; ++i)
    pthread_join (threads[i], 0) ;             // wait for threads to complete
  newFree (threads, nThreads, pthread_t) ;

  for (i = 1 ; i < nThreads ; ++i)
    for (j = 0 ; j < arrayMax(tx->nodes) ; ++j)
      ta->txUsed[j] |= ta[i].txUsed[j] ;
  oneReaderWriteTx (or, tx, ta->txUsed) ;

  for (i = 1 ; i < nThreads ; ++i)
    { taxLCAdestroy (ta[i].txLCA) ;
      newFree (ta[i].txUsed, arrayMax(tx->nodes), bool) ;
    }

  oneReaderDestroy (or) ;
  
  if (!outFileName)
    if (rename (outName, inFileName))  // rename() returns 0 on success, -1 on fail
      die ("failed to rename temp file %s to %s", outName, inFileName) ;

  return true ;
}

/*  File: onebamtax.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 16 00:26 2025 (rd109)
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

/*********************** reporting *************************/

typedef struct {
  int      n ;
  TaxID    txid ;
  I64      dustSum ;
  I64      totLen ;
  I64      scoreSum ;
  int      nACGT[5] ; 		// reference base count: use 4 for all non-ACGT
  int      nC1, nC1toT ;	// reference C in position 1, C->T in position 1
  int      counts[1<<14] ;
} ReportStats ;

#include <math.h>

bool reportLCA (char *readFileName, char *outFileName, char *rankText, char *groupText)
{
  OneReader *or = oneReaderCreate (readFileName, 1) ;
  if (!or) die ("failed to open .1read file %s", readFileName) ;
  Taxonomy *tax = or->taxonomy ;
  if (!tax) die ("missing taxonomy - you must run addLCA first before reportLCA\n") ;
  if (!outFileName) outFileName = derivedName (readFileName, "reportLCA") ;
  FILE *out = fopen (outFileName, "w") ;
  if (!out) die ("failed to open output file %s", outFileName) ;

  TaxID *rankID = 0 ;
  if (rankText)
    { U32 rank ;
      if (!dictFind (tax->rankDict, rankText, &rank))
	die ("unrecognised taxonomic rank %s", rankText) ;
      rankID = new0 (arrayMax(tax->nodes), TaxID) ;
      TaxID    t = tax->root ;
      TaxNode *x ;
      while (true)
	{ x = arrp(tax->nodes, t, TaxNode) ;
	  rankID[t] = (x->rank == rank) ? t : rankID[x->parent] ; // inherit if not at rank
	  // now recurse
	  if (x->left)
	    t = x->left ; // next level down
	  else            // go right (next) or up (parent) then right
	    { while (!x->next && t != tax->root)
		{ t = x->parent ; x = arrp(tax->nodes, t, TaxNode) ; }
	      if (t == tax->root) break ; // we are done - leave main loop
	      t = x->next ;
	    }
	}
    }

  U32 group = 0 ;
  if (groupText && !dictFind (tax->groupDict, groupText, &group))
    die ("group %s not one of animals, plants, fungi, bacteria, archaea", groupText) ;
  
  Hash  taxHash = hashCreate (4096) ;
  Array stats   = arrayCreate (4096, ReportStats) ;

  int i, baseMap[256] ;
  for (i = 0 ; i < 256 ; ++i) baseMap[i] = 4 ;
  baseMap['a'] = 0 ; baseMap['c'] = 1 ; baseMap['g'] = 2 ; baseMap['t'] = 3 ;
  baseMap['A'] = 0 ; baseMap['C'] = 1 ; baseMap['G'] = 2 ; baseMap['T'] = 3 ;

  double log2small[1<<16] ;
  for (i = 0 ; i < 1<<16 ; ++i) log2small[i] = log2(1.0*i) ;

  int j, flatten[1<<14] ;
  for (i = 0 ; i < 1<<14 ; ++i)
    { j = ((3-(i & 3)) << 12) | ((3-((i>>2) & 3)) << 10) | ((3-((i>>4) & 3)) << 8) |
	((3-((i>>6) & 3)) << 6) | ((3-((i>>8) & 3)) << 4) | ((3-((i>>10) & 3)) << 2) | (3-(i>>12)) ;
      if (i < j) flatten[i] = i ;
      else flatten[i] = j ;
    }

  int mask = (1<<14) - 1 ;
  while (oneReaderNext (or))
    { TaxID txid = or->lca ;
      if (rankID) txid = rankID[txid] ;
      if (!txid) continue ;
      if (group && arrp(tax->nodes, txid, TaxNode)->group != group) continue ;
      int k ;
      hashAdd (taxHash, hashInt(txid), &k) ;
      ReportStats *rs = arrayp(stats, k, ReportStats) ;
      ++rs->n ;
      rs->txid = txid ; // only needed the first time, but faster to set than test
      rs->totLen += or->seqLen ;
      rs->dustSum += or->dustScore ;
      rs->scoreSum += or->maxScore ;
      if (*or->mLine == 'C' || (*or->mLine == '.' && *or->seq == 'c')) ++rs->nC1 ;
      if (*or->mLine == 'C' && *or->seq == 't') ++rs->nC1toT ;
      k = 0 ;
      for (i = 0 ; i < or->seqLen ; ++i)
	{ if (or->mLine[i] == '.') or->mLine[i] = or->seq[i] ; // map up to ref
	  ++rs->nACGT[baseMap[or->mLine[i]]] ;
	  k = ((k << 2) & mask) | baseMap[or->mLine[i]] ; // not right at -, but so what
	  if (i > 5) ++rs->counts[k] ;
	}
    }
  hashDestroy (taxHash) ;
  if (rankID) newFree (rankID, arrayMax(tax->nodes), TaxID) ;
  oneReaderDestroy (or) ;

  fprintf (out, "#tax_id\ttax_name\ttax_rank\tcount\t") ;
  fprintf (out, "av_len\tav_score\tav_dust\tav_GC\t5'CtoT\tHscore\tfamily\n") ;
  for (i = 0 ; i < arrayMax(stats) ; ++i)
    { ReportStats *rs = arrp(stats,i,ReportStats) ;
      TaxNode *tx = arrp(tax->nodes, rs->txid, TaxNode) ;
      fprintf (out, "%d\t", rs->txid) ;
      fprintf (out, "%s\t", dictName (tax->nameDict, rs->txid)) ;
      fprintf (out, "%s\t", dictName (tax->rankDict, tx->rank)) ;
      fprintf (out, "%d\t", rs->n) ;
      fprintf (out, "%.1f\t", rs->totLen / (1.0*rs->n)) ;
      fprintf (out, "%.1f\t", rs->scoreSum / (1.0*rs->n)) ;
      fprintf (out, "%.2f\t", rs->dustSum / (1.0*rs->n)) ;
      fprintf (out, "%.3f\t", (rs->nACGT[1]+rs->nACGT[2]) /
	       (1.0*(rs->nACGT[0]+rs->nACGT[1]+rs->nACGT[2]+rs->nACGT[3]))) ;
      fprintf (out, "%.2f %d/%d\t", rs->nC1?rs->nC1toT/(1.0*rs->nC1):0.0, rs->nC1toT, rs->nC1) ;
      int sum = 0, *count = rs->counts ;
      for (j = 0 ; j < 1<<14 ; ++j)
	if (flatten[j] < j) { count[flatten[j]] += count[j] ; count[j] = 0 ; }
      double info = 0.0 ;
      for (j = 0 ; j < 1<<14 ; ++j, ++count)
	if (*count)
	  { sum += *count ;
	    if (*count < 1<<16) info += *count * log2small[*count] ;
	    else info += *count * log2(*count) ;
	  }
      if (sum > 1<<16) fprintf (out, "%3d\t", (int)(1000*(log2(sum) - info/sum)/13)) ;
      else if (sum) fprintf (out, "%3d\t", (int)(1000*(log2small[sum] - info/sum)/13)) ;
      else fprintf (out, "0.00\t") ;
      fprintf (out, "%s\n", tx->family ? dictName (tax->nameDict, tx->family) : "-") ;
    }
  arrayDestroy (stats) ;
  fclose (out) ;
  return true ;
}

// sed -e 's/\t/,/g' ERR10493337.reportLCA | sort -t , -k4nr | column -s , -t | less -S
// cat Z3.reportLCA | sed -e 's/\t/,/g' | awk 'BEGIN{FS=","}($4 > 10)' | sort -t , -k8nr | column -s , -t | less -S



/*  File: merge.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: implementations of Gene Myers heap merging algorithms
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 10 15:00 2025 (rd109)
 * * Aug 10 11:51 2025 (rd109): see TEST in merge.c for example usage
 * Created: Sat Aug  9 23:19:03 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

#ifndef MERGE_DEFINED
typedef void Merge ;
#endif

Merge *mergeCreateI32 (int T, void *x, bool (*yield)(int t, void* x, I32 *v)) ;
Merge *mergeCreateString (int T, void *x, bool (*yield)(int t, void* x, char **v)) ;
void   mergeDestroy (Merge *m) ;
void   mergeRecreate (Merge *m, int T, void *x) ; // reset using existing allocated structures
// the puropose of mergeRecreate is to allow lightweight use for many merges of small sets
int    mergeNext (Merge *m, int **tList) ; // returns number of items with next smallest value
                                           // *tList contains a list of the input indices

// see TEST code at end of merge.c for example usage

/************** end of file *****************/

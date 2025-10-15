/*  File: merge.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: implementations of Gene Myers heap merging algorithms
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 15 12:09 2025 (rd109)
 * * Oct  6 14:29 2025 (rd109): see https://drops.dagstuhl.de/storage/00lipics/lipics-vol259-cpm2023/LIPIcs.CPM.2023.22/LIPIcs.CPM.2023.22.pdf
 * * Aug 10 11:51 2025 (rd109): see TEST in merge.c for example usage
 * Created: Sat Aug  9 23:19:03 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"

#ifndef MERGE_DEFINED
typedef void Merge ;
#endif

Merge *mergeCreateI32 (int T, void *arg, bool (*yield)(int t, void* arg, I32 *v)) ;
// yield is a callback to return true with the next value on input stream t in *v, or false at end
Merge *mergeCreateString (int T, void *arg, bool (*yield)(int t, void* arg, char **v, int *lcp)) ;
// yield returns the next full string in **v, with optionally the lcp to the last value in *lcp
void   mergeDestroy (Merge *m) ;
void   mergeRecreate (Merge *m, int T, void *arg) ; // reset using existing allocated structures
// the puropose of mergeRecreate is to allow lightweight use for many merges of small sets
int    mergeNext (Merge *m, int **tList) ; // returns number of items with next smallest value
                                           // *tList contains a list of the input indices
bool   mergeUpdateString (Merge *m, int t, char *sNew) ;
// this allows the user to move the location of the last yielded string for stream t to sNew
// if the string contents at the new and old location are the same it accepts and returns true 
// if they are different it returns false and does not update - strdup your string first!

// see TEST code at end of merge.c for example usage

/************** end of file *****************/

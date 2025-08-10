/*  File: merge.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: implementations of Gene Myers heap merging algorithms
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 10 21:34 2025 (rd109)
 * Created: Sat Aug  9 18:05:50 2025 (rd109)
 *-------------------------------------------------------------------
 */


#include "utils.h"

typedef struct {
  int     T ;     // the number of inputs
  int     Tmax ;  // allocated 
  void   *x ;     // opaque argument to pass to callback
  int    *H ;     // 1..T - Gene's H: the heap - values are input indices
  I32    *V ;     // 0..T - Gene's V: current value of the t'th input
  bool   *L, *R ; // 1..T - Gene's L and R: is equal to left child, right child
  int     nG ;    // Gene's len: number in next cohort of equal-valued samples
  int    *G ;     // 0..T-1 - Gene's G offset by 1: heap indices for next cohort of samples
  int    *HG ;    // 0..T-1 - H[G[i]: used to return tList in mergeNext
  bool  (*yieldI32)(int t, void *x, I32 *v) ;
                  // return true with next value of input t in *v, or false if end of input t
} Merge ;
#define MERGE_DEFINED   // prevents typedef of Merge to void

#include "merge.h"

static void printHeap (Merge *m) // for debug
{
  int t ;
  printf ("  HEAP") ;
  for (t = 0 ; t <= m->T ; ++t) printf (" (%d %d:%d)", t, m->H[t], m->V[m->H[t]]) ;
  putchar ('\n') ;
}

// Heapify functions insert x in stream t starting at i

static void heapifyI32 (Merge *m, int i, I32 x, int t) // collision heap from Figure 3
{
  int    c = i, l ;
  int    hl, hr ; // left and right indexes
  I32    vl, vr ; // current left and right values

  //  printHeap (m) ;
  m->V[t] = x ;
  while ((l = 2*c) <= m->T)
    { int hl = m->H[l] ;
      int vl = m->V[hl] ;
      if (l < m->T)
	{ hr = m->H[l+1] ;
	  vr = m->V[hr] ;
	}
      else
	vr = I32MAX ;
      if (vl < vr) // so we will go down left route, if anything
	{ m->R[c] = false ;
	  if (x > vl) // Case 1L - swap down left with hl
	    { m->H[c] = hl ;
	      m->L[c] = m->L[l] || m->R[l] ;
	      c = l ;
	    }
	  else if (x == vl) // Case 2L - done
	    { m->H[c] = t ;
	      m->L[c] = true ;
	      return ;
	    }
	  else // x < vl Case 3L - done
	    { m->H[c] = t ;
	      m->L[c] = false ;
	      return ;
	    }
	}
      else if (vr < vl) // symmetric Cases 1R, 2R, 3R
	{ m->L[c] = false ;
	  if (x > vr) // Case 1R - swap down right with hr
	    { m->H[c] = hr ;
	      m->R[c] = m->L[l+1] || m->R[l+1] ;
	      c = l+1 ;
	    }
	  else if (x == vr) // Case 2R - done
	    { m->H[c] = t ;
	      m->R[c] = true ;
	      return ;
	    }
	  else // x < vr - case 3R - done
	    { m->H[c] = t ;
	      m->R[c] = false ;
	      return ;
	    }
	}
      else // vr == vl
	{ if (x > vl) // Case 4 - swap down left
	    { m->H[c] = hl ;
	      m->L[c] = m->L[l] || m->R[l] ;
	      m->R[c] = true ;
	      c = l ;
	    }
	  else if (x == vl) // Case 5 - done
	    { m->H[c] = t ;
	      m->L[c] = m->R[c] = true ;
	      return ;
	    }
	  else // x < vl - Case 3M - done
	    { m->H[c] = t ;
	      m->L[c] = m->R[c] = false ;
	      return ;
	    }
	}
    }
  // if we get here we have hit the end
  m->H[c] = t ;
  m->L[c] = m->R[c] = false ;
} 

/** next the code to collect the next set of samples */

static int addToCohort (Merge *m, int c, int len)
{ if (m->L[c]) len = addToCohort (m, 2*c, len) ;
  if (m->R[c]) len = addToCohort (m, 2*c+1, len) ;
  m->G[len] = c ;
  m->HG[len] = m->H[c] - 1 ; // we return values as 0..T-1
  return ++len ;
}

/** and initialisation - set up so the heap fills from bottom to top on first call to mergeNext() */

static void initialise (Merge *m)
{
  int t, T = m->T ;
  for (t = 0 ; t < T ; ++t)
    { m->H[t+1] = t+1 ;
      m->G[t]   = T - t ; // last element first
      m->HG[t]  = T - t - 1 ;
    }
  m->nG = T ;
}

/********** public interface ************/

Merge *mergeCreateI32 (int T, void *x, bool (*yield)(int t, void *x, I32 *v))
{
  Merge *m = new (1, Merge) ;
  m->T  = T ;
  m->x  = x ;
  m->Tmax = T ;
  m->H  = new0 (T+1, int) ;
  m->V  = new0 (T+1, I32) ;
  m->L  = new0 (T+1, bool) ;
  m->R  = new0 (T+1, bool) ;
  m->G  = new0 (T, int) ;
  m->HG  = new0 (T, int) ;
  m->yieldI32 = yield ;
  initialise (m) ;
  return m ;
}

void mergeRecreate (Merge *m, int T, void *x)
{
  if (T > m->Tmax) // must reallocate
    { newFree (m->H, m->Tmax+1, int) ;
      if (m->V) newFree (m->V, m->Tmax+1, I32) ;
      newFree (m->L, m->Tmax+1, bool) ;
      newFree (m->R, m->Tmax+1, bool) ;
      newFree (m->G, m->Tmax+1, int) ;
      newFree (m->HG, m->Tmax, int) ;
      m->Tmax = T ;
      m->H  = new0 (T+1, int) ;
      if (m->V) m->V = new0 (T+1, I32) ;
      m->L  = new0 (T+1, bool) ;
      m->R  = new0 (T+1, bool) ;
      m->G  = new0 (T, int) ;
      m->HG  = new0 (T, int) ;
    }
  m->T = T ;
  m->x = x ;
  initialise (m) ;
}

void mergeDestroy (Merge *m)
{
  newFree (m->H, m->Tmax+1, int) ;
  if (m->V) newFree (m->V, m->Tmax+1, I32) ;
  newFree (m->L, m->Tmax+1, bool) ;
  newFree (m->R, m->Tmax+1, bool) ;
  newFree (m->G, m->Tmax, int) ;
  newFree (m->HG, m->Tmax, int) ;
  newFree (m, 1, Merge) ;
}

int mergeNext (Merge *m, int **tList) // returns the number of entries in m->G, 0 at finish
{
  int k ;
  if (m->V)
    { I32 v ;
      // printf ("entering mergeNext nG %d\n", m->nG) ;
      for (k = 0 ; k < m->nG ; ++k)
	{ bool res = (m->yieldI32)(m->HG[k], m->x, &v) ;
	  //	  printf ("  yield k %d HG[k] %d G[k] %d v %d\n", k, m->HG[k], m->G[k], v) ;
	  if (res)
	    heapifyI32 (m, m->G[k], v, m->HG[k]+1) ; // recall that HG is one less than H[G]
	  else
	    heapifyI32 (m, m->G[k], I32MAX, 0) ;
	}
    }
  if (!m->H[1]) return 0 ; // we are done
  m->nG = addToCohort (m, 1, 0) ; // add 1 and anything below with equal value to return list
  if (tList) *tList = m->HG ;
  return m->nG ;
}

/************************* test package ************************/

#ifdef TESTI32

typedef struct {
  int   T ;   // number of lists
  int  *n ;   // length of each list
  int  *i ;   // current position in each list (next item to be yielded, starts at 0)
  int **val ; // the lists themselves
} ListSet ;

static bool yield (int t, void *x, I32 *v)
{ ListSet *ls = (ListSet*)x ;
  //  printf ("  yield called on input %d\n", t) ;
  //  printf ("    ls->i[t] %d ls->n[t] %d\n", ls->i[t], ls->n[t]) ;
  if (ls->i[t] < ls->n[t])
    { *v = ls->val[t][(ls->i[t])++] ; return true ; }
  else
    return false ;
}

static int intOrder (const void *a, const void *b)
{ int ia = *(int*)a, ib = *(int*)b ; return (ia-ib) ; }

int main (int argc, char *argv[])
{
  if (argc != 4) die ("usage: %s <nList> <nItemPerList> <valueMax>", argv[0]) ;
  --argc ; ++argv ; // swallow the program name

  // set up the list set
  ListSet *ls = new0 (1, ListSet) ;
  int T = atoi(*argv++), N = atoi(*argv++), max = atoi(*argv++) ;
  ls->T   = T ;
  ls->n   = new (T, int) ;
  ls->i = new0 (T, int) ; // initialise to 0
  ls->val = new (T, int*) ;
  int t, i ;
  for (t = 0 ; t < T ; ++t)
    { ls->val[t] = new (N, int) ;
      for (i = 0 ; i < N ; ++i) ls->val[t][i] = random() % max ;
      qsort (ls->val[t], N, sizeof(int), intOrder) ; // sort
      int n = 0 ; // and make unique
      for (i = 1 ; i < N ; ++i) if (ls->val[t][i] != ls->val[t][n]) ls->val[t][++n] = ls->val[t][i] ;
      ls->n[t] = n+1 ;
    }
  printf ("INPUTS\n") ;
  for (i = 0 ; i < N ; ++i)
    { for (t = 0 ; t < T ; ++t)
	if (i < ls->n[t]) printf (" %2d", ls->val[t][i]) ;
	else printf ("   ") ;
      putchar ('\n') ;
    }

  Merge *m = mergeCreateI32 (T, ls, yield) ;
  int n, *tList ;
  printf ("OUTPUT\n") ;
  while ((n = mergeNext (m, &tList)))
    { // printHeap (m) ;
      int t = tList[0] ;
      printf (" output value %2d count %d inputs", ls->val[t][ls->i[t]-1], n) ;
      int i ; for (i = 0 ; i < n ; ++i) printf (" %d", tList[i]) ;
      putchar ('\n') ;
    }
}

#endif

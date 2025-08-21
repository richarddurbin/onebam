/*  File: merge.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: implementations of Gene Myers heap merging algorithms
 * Exported functions:
 * HISTORY:
 * Last edited: Aug 21 23:36 2025 (rd109)
 * Created: Sat Aug  9 18:05:50 2025 (rd109)
 *-------------------------------------------------------------------
 */


#include "utils.h"

typedef struct {
  int     T ;     // the number of inputs
  int     Tmax ;  // allocated 
  void   *arg ;   // opaque argument to pass to callback
  int    *H ;     // 1..T - Gene's H: the heap - values are input indices
  int     nG ;    // Gene's len: number in next cohort of equal-valued samples
  int    *G ;     // 0..T-1 - Gene's G offset by 1: heap indices for next cohort of samples
  int    *HG ;    // 0..T-1 - H[G[i]: used to return tList in mergeNext
  bool  (*yieldI32)(int t, void *arg, I32 *v) ; // client function to give next value of input t
  I32    *V ;     // 0..T - Gene's collision (int) V: current value of the t'th input
  bool   *L, *R ; // 1..T - Gene's L and R: is equal to left child, right child
  bool  (*yieldString)(int t, void *arg, char **v, int *p) ; // same for strings
  char  **S ;     // 1..T - Gene's string V
  int    *P ;     // 1..T - Gene's P: the LCP to the parent
} Merge ;
#define MERGE_DEFINED   // prevents typedef of Merge to void

#include "merge.h"

static void printHeapI32 (Merge *m) // for debug
{
  int t ;
  printf ("  HEAP") ;
  for (t = 0 ; t <= m->T ; ++t) printf (" (%d %d:%d)", t, m->H[t], m->V[m->H[t]]) ;
  putchar ('\n') ;
}

// Heapify functions insert x in stream t starting at i

static void heapifyI32 (Merge *m, int i, I32 x, int t) // string heap from Figure 2 with V->S, v->s
{
  int    c = i ;
  int    l, r ; // left, right
  int    hl, hr ; // left and right indexes
  I32    vl, vr ; // current left and right values

  //  printHeapI32 (m) ;
  m->V[t] = x ;
  while ((l = 2*c) <= m->T)
    { r = l+1 ;
      hl = m->H[l] ; vl = m->V[hl] ;
      if (l < m->T) { hr = m->H[r] ; vr = m->V[hr] ; } else vr = I32MAX ;
      if (vl < vr) // so we will go down left route, if anything
	{ m->R[c] = false ;
	  if (x > vl) // Case 1L - swap down left with hl
	    { m->H[c] = hl ; m->L[c] = m->L[l] || m->R[l] ; c = l ; }
	  else if (x == vl) // Case 2L - done
	    { m->H[c] = t ; m->L[c] = true ; return ; }
	  else // x < vl Case 3L - done
	    { m->H[c] = t ; m->L[c] = false ; return ; }
	}
      else if (vr < vl) // symmetric Cases 1R, 2R, 3R
	{ m->L[c] = false ;
	  if (x > vr) // Case 1R - swap down right with hr
	    { m->H[c] = hr ; m->R[c] = m->L[r] || m->R[r] ; c = r ; }
	  else if (x == vr) // Case 2R - done
	    { m->H[c] = t ; m->R[c] = true ; return ; }
	  else // x < vr - case 3R - done
	    { m->H[c] = t ; m->R[c] = false ; return ; }
	}
      else // vr == vl
	{ if (x > vl) // Case 4 - swap down left
	    { m->H[c] = hl ; m->L[c] = m->L[l] || m->R[l] ; m->R[c] = true ; c = l ; }
	  else if (x == vl) // Case 5 - done
	    { m->H[c] = t ; m->L[c] = m->R[c] = true ; return ; }
	  else // x < vl - Case 3M - done
	    { m->H[c] = t ; m->L[c] = m->R[c] = false ; return ; }
	}
    }
  // if we get here we have hit the end
  m->H[c] = t ;
  m->L[c] = m->R[c] = false ;
} 

/** next the code to collect the next set of samples */

static int addToCohortI32 (Merge *m, int c, int len)
{ if (m->L[c]) len = addToCohortI32 (m, 2*c, len) ;
  if (m->R[c]) len = addToCohortI32 (m, 2*c+1, len) ;
  m->G[len] = c ;
  m->HG[len] = m->H[c] - 1 ; // we return values as 0..T-1
  return ++len ;
}

/** now the same routines for strings - first a couple of LCP functions */

static inline int lcp2 (char *x, char *y, int n)
{ while (true)
    if (!x[n] || x[n] != y[n]) return n ;
    else ++n ;
}

static inline int lcp3 (char *x, char *y, char *z, int n)
{ while (true)
    if (!x[n] || x[n] != y[n] || x[n] != z[n]) return n ;
    else ++n ;
} 

static void printHeapString (Merge *m) // for debug
{
  int t ;
  printf ("  HEAP") ;
  for (t = 1 ; t <= m->T ; ++t) printf (" (%d:%s)", m->H[t], m->S[m->H[t]]) ;
  putchar ('\n') ;
}

#define DEBUG(x) // printf(x)

static void heapifyString (Merge *m, int i, char *x, int t, int p) // collision heap from Figure 3
{
  int    l,   r ;  // left and right
  int    hl,  hr ; // left and right indexes
  char  *sl, *sr ; // current left and right strings
  int    pl,  pr ; // LCP of current to left and right children

  //  printf ("heapify %s t %d at i %d", x, t, i) ;
  
  int    c = i ;

  if (p == -1) p = lcp2 (m->S[m->H[i]],x,0) ; // the yieldString call did not pass the lcp
  
  while ((l = 2*c) <= m->T) // we will break once we find where to place x
    { hl = m->H[l] ; pl = m->P[l] ;
      r = l+1 ;
      if (l < m->T) { hr = m->H[r] ; pr = m->P[r] ; } else pr = -1 ; // don't need hr in this case
      if (pr < pl)  // we will go down left route, if anything
	{ if (p < pl)                                   // case 1L - swap with left and iterate
	    { m->H[c] = hl ; m->P[c] = pl ; c = l ; DEBUG(" case 1L\n") ; }
	  else if (p > pl)                              // case 2L - done
	    { DEBUG(" case 2L\n") ; break ; }
	  else // p == pl
	    { sl = m->S[hl] ;
	      int px = lcp2 (sl, x, pl) ;
	      if (sl[px] < x[px])                       // Case 3La - swap with left and iterate
		{ m->H[c] = hl ; m->P[c] = pl ; c = l ; p = px ; DEBUG(" case 3La\n") ; }
	      else                                      // Case 3Lb - set m->P[l] = px then done
		{ m->P[l] = px ; DEBUG(" case 3Lb\n") ; break ; }
	    }
	}
      else if (pr > pl) // symmetric - explore going down right route (can't happen if pr == -1)
	{ if (p < pr)                                   // case 1R
	    { m->H[c] = hr ; m->P[c] = pr ; c = r ; DEBUG(" case 1R\n") ; }
	  else if (p > pr)                              // case 2R
	    { DEBUG(" case 2R\n") ; break ; }
	  else // p == pr
	    { sr = m->S[hr] ;
	      int px = lcp2 (sr, x, pr) ;
	      if (sr[px] < x[px])                       // case 3Ra
		{ m->H[c] = hr ; m->P[c] = pr ; c = r ; p = px ; DEBUG(" case 3Ra\n") ; }
	      else                                      // case 3Rb
		{ m->P[r] = px ; DEBUG(" case 3Rb\n") ; break ; }
	    }
	}
      else             // pl == pr
	{ if (p > pl)                                   // case 2 - done
	    { DEBUG(" case 2.3\n") ; break ; }
	  else
	    { sl = m->S[hl] ; sr = m->S[hr] ;
	      if (p < pl)                               // case 4
		{ int px = lcp2 (sr, sl, pl) ;
		  if (sl[px] <= sr[px]) // go left
		    { m->H[c] = hl ; m->P[c] = pl ; m->P[r] = px ; c = l ; DEBUG(" case 4L\n") ; }
		  else // go right
		    { m->H[c] = hr ; m->P[c] = pr ; m->P[l] = px ; c = r ; DEBUG(" case 4R\n") ; }
		}
	      else // p == pl == pr                    // case 5 - lots of subcases
		{ int px = lcp3 (sl, sr, x, p) ;
		  if (sr[px] > sl[px])
		    { if (x[px] > sl[px])              // case 5.1L
			{ m->H[c] = hl ; m->P[c] = pl ; m->P[r] = px ; c = l ; p = px ; DEBUG(" case 5.1L\n") ; }
		      else if (x[px] < sl[px])         // case 5.2-1
			{ m->P[l] = m->P[r] = px ; DEBUG(" case 5.2-1\n") ; break ; }
		      else
			{ int py = lcp2 (sl, x, px) ;
			  if (sl[py] < x[py])          // case 5.3La
			    { m->H[c] = hl ; m->P[c] = pl ; m->P[r] = px ; c = l ; p = py ; DEBUG(" case 5.3La\n") ; }
			  else                         // case 5.3Lb
			    { m->P[l] = py ; m->P[r] = px ; DEBUG(" case 5.3Lb\n") ; break ; }
			}
		    }
		  else if (sr[px] < sl[px])
		    { if (x[px] > sr[px])              // case 5.1R
			{ m->H[c] = hr ; m->P[c] = pr ; m->P[l] = px ; c = r ; p = px ; DEBUG(" case 5.1R\n") ; }
		      else if (x[px] < sr[px])         // case 5.2 again
			{ m->P[l] = m->P[r] = px ; DEBUG(" case 5.2-2\n") ; break ; }
		      else
			{ int py = lcp2(sr,x,px) ;
			  if (sr[py] < x[py])          // case 5.3Ra
			    { m->H[c] = hr ; m->P[c] = pr ; m->P[l] = px ; c = r ; p = py ; DEBUG(" case 5.3Ra\n") ; }
			  else                         // case 5.3Rb
			    { m->P[r] = py ; m->P[l] = px ; DEBUG(" case 5.3Rb\n") ; break ; }
			}
		    }
		  else // sr[px] == sl[px]
		    { if (x[px] <= sl[px])              // case 5.2 again
			{ m->P[l] = m->P[r] = px ; DEBUG(" case 5.2-3\n") ; break ; }
		      else                             // case 5.4
			{ int py = lcp2 (sl, sr, px) ;
			  if (sl[py] < sr[py])         // case 5.4L
			    { m->H[c] = hl ; m->P[c] = pl ; m->P[r] = py ; c = l ; p = px ; DEBUG(" case 5.4L\n") ; }
			  else                         // case 5.4R
			    { m->H[c] = hr ; m->P[c] = pr ; m->P[l] = py ; c = r ; p = px ; DEBUG(" case 5.4R\n") ; }
			}
		    }
		}
	    }
	}
    }
  m->H[c] = t ; m->P[c] = p ; m->S[t] = x ;  // finally place x at c
  //  printHeapString (m) ;
} 

/** next the code to collect the next set of samples */

static int addToCohortString (Merge *m, int c, int len)
{ if (2*c <= m->T && !m->S[m->H[2*c]][m->P[2*c]]) len = addToCohortString (m, 2*c, len) ; // full string match
  if (2*c+1 <= m->T && !m->S[m->H[2*c+1]][m->P[2*c+1]]) len = addToCohortString (m, 2*c+1, len) ;
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

static Merge *mCreate (int T, void *arg)
{
  Merge *m = new0 (1, Merge) ;
  m->T    = T ;
  m->arg  = arg ;
  m->Tmax = T ;
  m->H    = new0 (T+1, int) ;
  m->G    = new0 (T, int) ;
  m->HG   = new0 (T, int) ;
  initialise (m) ;
  return m ;
}


/********** public interface ************/

Merge *mergeCreateI32 (int T, void *arg, bool (*yield)(int t, void *x, I32 *v))
{
  Merge *m = mCreate (T, arg) ;
  m->V = new0 (T+1, I32) ;
  m->L = new0 (T+1, bool) ;
  m->R = new0 (T+1, bool) ;
  m->yieldI32 = yield ;
  return m ;
}

Merge *mergeCreateString (int T, void *arg, bool (*yield)(int t, void *x, char **s, int *p))
{
  Merge *m = mCreate (T, arg) ;
  m->S = new0 (T+1, char*) ;
  int t ; for (t = 0 ; t <= T ; ++t) m->S[t] = "" ; // point to empty string
  m->P = new0 (T+1, int) ;
  m->yieldString = yield ;
  return m ;
}

void mergeRecreate (Merge *m, int T, void *arg)
{
  if (T > m->Tmax) // must reallocate
    { newFree (m->H, m->Tmax+1, int) ;
      newFree (m->G, m->Tmax+1, int) ;
      newFree (m->HG, m->Tmax, int) ;
      if (m->V)
	{ newFree (m->V, m->Tmax+1, I32) ;
	  newFree (m->L, m->Tmax+1, bool) ;
	  newFree (m->R, m->Tmax+1, bool) ;
	}
      else
	{ newFree (m->S, m->Tmax+1, char*) ;
	  newFree (m->P, m->Tmax+1, int) ;
	}
      m->Tmax = T ;
      m->H  = new0 (T+1, int) ;
      if (m->V)
	{ m->V = new0 (T+1, I32) ;
	  m->L  = new0 (T+1, bool) ;
	  m->R  = new0 (T+1, bool) ;
	}
      else
	{ m->S = new0 (m->Tmax+1, char*) ;
	  m->P = new0 (m->Tmax+1, int) ;
	}
      m->G  = new0 (T, int) ;
      m->HG  = new0 (T, int) ;
    }
  m->T = T ;
  m->arg = arg ;
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
  static char LAST[2] = { 0xff, 0 } ;
  int k ;
  if (m->V)
    { I32 x ;
      // printf ("entering mergeNext I32 nG %d\n", m->nG) ;
      for (k = 0 ; k < m->nG ; ++k)
	{ bool res = (m->yieldI32)(m->HG[k], m->arg, &x) ;
	  //	  printf ("  yield k %d HG[k] %d G[k] %d x %d\n", k, m->HG[k], m->G[k], x) ;
	  if (res)
	    heapifyI32 (m, m->G[k], x, m->HG[k]+1) ; // recall that HG is one less than H[G]
	  else
	    heapifyI32 (m, m->G[k], I32MAX, 0) ;
	}
      if (!m->H[1]) return 0 ; // we are done
      m->nG = addToCohortI32 (m, 1, 0) ; // add 1 and anything below with equal value to return list
    }
  else // string
    { char *s ;
      // printf ("entering mergeNext string nG %d\n", m->nG) ;
      for (k = 0 ; k < m->nG ; ++k)
	{ int p = -1 ;
	  bool res = (m->yieldString)(m->HG[k], m->arg, &s, &p) ;
	  // printf ("  yield k %d HG[k] %d G[k] %d s %s\n", k, m->HG[k], m->G[k], s) ;
	  if (res)
	    heapifyString (m, m->G[k], s, m->HG[k]+1, p) ; // recall that HG is one less than H[G]
	  else
	    heapifyString (m, m->G[k], LAST, 0, 0) ;
	}
      if (!m->H[1]) return 0 ; // we are done
      m->nG = addToCohortString (m, 1, 0) ; // add 1 and anything below with equal value to return list
    }
  if (tList) *tList = m->HG ;
  // printf ("  leaving mergeNext\n") ;
  return m->nG ;
}

/************************* test packages ************************/

#ifdef TEST_I32

typedef struct {
  int   T ;   // number of lists
  int  *n ;   // length of each list
  int  *i ;   // current position in each list (next item to be yielded, starts at 0)
  int **val ; // the lists themselves
} ListSet ;

static bool yield (int t, void *arg, I32 *v)
{ ListSet *ls = (ListSet*)arg ;
  // printf ("  yield called on input %d\n", t) ;
  // printf ("    ls->i[t] %d ls->n[t] %d\n", ls->i[t], ls->n[t]) ;
  if (ls->i[t] < ls->n[t])
    { *v = ls->val[t][(ls->i[t])++] ; return true ; }
  else
    return false ;
}

static int intOrder (const void *a, const void *b)
{ int ia = *(int*)a, ib = *(int*)b ; return (ia-ib) ; }

int main (int argc, char *argv[])
{
  if (argc != 4) die ("usage: %s <nList> <nItemPerList> <valueMax> # e.g 10 10 50", argv[0]) ;
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
    { // printHeapI32 (m) ;
      int t = tList[0] ;
      printf (" output value %2d count %d inputs", ls->val[t][ls->i[t]-1], n) ;
      int i ; for (i = 0 ; i < n ; ++i) printf (" %d", tList[i]) ;
      putchar ('\n') ;
    }
}

#endif

#ifdef TEST_STRING

typedef struct {
  int     T ;   // number of lists
  int    *n ;   // length of each list
  int    *i ;   // current position in each list (next item to be yielded, starts at 0)
  char ***val ; // the lists themselves
} ListSet ;

static bool yield (int t, void *arg, char **s, int *p)
{ ListSet *ls = (ListSet*)arg ;
  //  printf ("  yield called on input %d\n", t) ;
  //  printf ("    ls->i[t] %d ls->n[t] %d\n", ls->i[t], ls->n[t]) ;
  if (ls->i[t] < ls->n[t])
    { *s = ls->val[t][(ls->i[t])++] ; return true ; }
  else
    return false ;
}

static int stringOrder (const void *a, const void *b)
{ char *sa = *(char**)a, *sb = *(char**)b ; return strcmp (sa, sb) ; }

// Function to generate a Poisson random variable with rate lambda
#include <math.h>
int poisson_random(double lambda)
{ int n = 0 ;
  double sum = 0.0 ;
  while (sum < lambda) { sum -= log ((double)rand() / RAND_MAX) ; n++ ; }
  return n - 1;
}

int main (int argc, char *argv[])
{
  if (argc != 4) die ("usage: %s <nList> <nItemPerList> <meanLen> # e.g 10 10 4", argv[0]) ;
  --argc ; ++argv ; // swallow the program name

  // set up the list set
  ListSet *ls = new0 (1, ListSet) ;
  int T = atoi(*argv++), N = atoi(*argv++) ;
  double lambda = atof(*argv++) ; if (lambda <= 0) die ("mean %f must be positive", lambda) ;
  ls->T   = T ;
  ls->n   = new (T, int) ;
  ls->i   = new0 (T, int) ; // initialise to 0
  ls->val = new (T, char**) ;
  int t, i ;
  for (t = 0 ; t < T ; ++t)
    { ls->val[t] = new (N, char*) ;
      for (i = 0 ; i < N ; ++i)
	{ int len = poisson_random(lambda) ;
	  char *s = ls->val[t][i] = new0(len+1,char) ; // inefficient
	  while (len--) *s++ = (random() & 0x1) ? 'a' : 'b' ;
	}
      qsort (ls->val[t], N, sizeof(char*), stringOrder) ; // sort
      int n = 0 ; // and make unique
      for (i = 1 ; i < N ; ++i)
	if (strcmp(ls->val[t][i],ls->val[t][n])) ls->val[t][++n] = ls->val[t][i] ;
      ls->n[t] = n+1 ;
    }
  printf ("INPUTS\n") ;
  for (i = 0 ; i < N ; ++i)
    { for (t = 0 ; t < T ; ++t)
	if (i < ls->n[t]) printf (" %-10s", ls->val[t][i]) ;
	else printf ("           ") ;
      putchar ('\n') ;
    }

  Merge *m = mergeCreateString (T, ls, yield) ;
  int n, *tList ;
  printf ("OUTPUT\n") ;
  while ((n = mergeNext (m, &tList)))
    { int t = tList[0] ;
      printf (" output value %s count %d inputs", ls->val[t][ls->i[t]-1], n) ;
      int i ; for (i = 0 ; i < n ; ++i) printf (" %d", tList[i]) ;
      putchar ('\n') ;
    }
}

#endif

/************************** end of file ************************/

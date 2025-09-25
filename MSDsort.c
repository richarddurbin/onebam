/*  Last edited: Sep 25 16:33 2025 (rd109) */
// Gene Myer's msd_sort()
// copyright Eugene Myers 2025 
// available under MIT licence

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <errno.h>

#include "gene_core.h"

#define IS_SORTED
#undef  REPORT_FITS
#undef  SHOW_LAUNCH

static int    NTHREADS;   //  Command line parameters
static int    RSIZE;
static int    KSIZE;
static int    MARK;
static int    DEPTH;

static pthread_t *THREADS;    //  NTHREADS thread records

#define SMAX  20
#define NMAX   3

#define THR0 15
#define THR1 15
#define THR2  8
#define GAP1  9
#define GAP2  4

static int S_thr0, S_thr1, S_thr2;
static int S_gap1, S_gap2;

static inline void mycpy(uint8 *a, uint8 *b, int n)
{ while (n--)
    *a++ = *b++;
}

static inline int mycmp(uint8 *a, uint8 *b, int n)
{ while (n--)
    { if (*a++ != *b++)
        return (a[-1] < b[-1] ? -1 : 1);
    }
  return (0);
}

#ifdef IS_SORTED

static inline void sorted(uint8 *array, int64 asize)
{ int64 p, i;

  array += DEPTH;
  for (p = RSIZE; p < asize; p += RSIZE)
    if (mycmp(array + (p-RSIZE),array + p,KSIZE-DEPTH) > 0)
      { printf("  Not sorted %12lld: ",p/RSIZE);
        for (i = 0; i < KSIZE-DEPTH; i++)
          printf(" %02x",array[(p-RSIZE)+i]);
        printf(" vs ");
        for (i = 0; i < KSIZE-DEPTH; i++)
          printf(" %02x",array[p+i]);
        printf("\n");
      }
}

#endif

static inline void gap_sort(uint8 *array, int asize, int gap, int cmp, int rem)
{ int    i, j;
  uint8  temp[RSIZE];
  uint8 *garray;

  garray = array + gap;
  for (i = gap; i < asize; i += RSIZE)
    { j = i-gap;
      if (mycmp(array+j,array+i,cmp) <= 0)
        continue;
      mycpy(temp,array+i,rem);
      mycpy(array+i,array+j,rem);
      for(j -= gap; j >= 0; j -= gap)
        { if (mycmp(array+j,temp,cmp) <= 0)
            break;
          mycpy(garray+j,array+j,rem);
        }
      mycpy(garray+j,temp,rem);
    }
}

static inline void shell_sort(uint8 *array, int asize, int digit)
{ int    cmp, rem;
  uint8 *garray;

  cmp    = KSIZE-digit;
  rem    = RSIZE-digit;
  garray = array+digit;

  // if (asize > S_thr1)
    // gap_sort(garray,asize,S_gap1,cmp,rem);
  if (asize > S_thr2)
    gap_sort(garray,asize,S_gap2,cmp,rem);
  gap_sort(garray,asize,RSIZE,cmp,rem);

  if (MARK)
    { int i, j;

      j = 0;
      for (i = RSIZE; i < asize; i += RSIZE)
        if (mycmp(garray+j,garray+i,cmp) != 0)
          { array[i] = 1;
            j = i;
          }
    }
}


/**********************************************************************************************
*
*  THREAD LAUNCHED RADIX SORT OF A SUBSECTION (IN MEMORY VERSION)
*
**********************************************************************************************/

static void radix_sort(uint8 *array, int64 asize, int depth, int64 *alive)
{ int64  n, len[256];
  int    y, q, ntop;
  int    nzero[256];

  { uint8 *end[256];
    uint8 *u, *arrow;
    int64  o;
    int    e, x, z, rem;

    uint8 *off[256];
    uint8  temp[RSIZE];
    uint8 *stack[SMAX];
    uint8 *copy;
    
    arrow = array + depth;
    while (1)
      { e = arrow[0];
        for (o = RSIZE; o < asize; o += RSIZE)
          { if (arrow[o] != e)
              break;
          }
        if (o < asize)
          break;
        depth += 1;
        if (depth >= KSIZE)
          return;
        arrow = array + depth;
      }
    
    rem  = RSIZE - depth;
    copy = temp + depth;
    
    ntop = 1;
    nzero[0] = e;
    alive[e] = o;
    for (; o < asize; o += RSIZE)
      { x = arrow[o];
        if (alive[x] == 0)
          nzero[ntop++] = x;
        alive[x] += RSIZE;
      }

    u = array;
    if (ntop <= NMAX)
      { for (y = 1; y < ntop; y++)
          { x = nzero[y];
            for (z = y-1; z >= 0; z--)
              if (nzero[z] < x)
                break;
              else
                nzero[z+1] = nzero[z];
            nzero[z+1] = x;
          }
        for (y = 0; y < ntop; y++)
          { x = nzero[y];
            len[x] = alive[x];
            alive[x] = 0;
            off[x] = u;
            end[x] = u += len[x];
          }
      }
    else
      { ntop = 0;
        for (x = 0; x < 256; x++)
          if (alive[x] > 0)
            { len[x] = alive[x];
              alive[x] = 0;
              off[x] = u;
              end[x] = u += len[x];
              nzero[ntop++] = x;
            }
      }

    for (y = 0; y < ntop; y++)
      { uint8   *p;
        int      t, s;
        int      z;

        x = nzero[y];
        while (off[x] < end[x])
          { t = off[x][depth];

            if (t == x)
              off[x] += RSIZE;
	          else
              { s = 0;
                stack[s++] = off[x];
                while (s < SMAX)
                  if (t == x)
                    { off[x] += RSIZE;
                      break;
                    }
                  else
                    { u = off[t];
                      while ((z = u[depth]) == t)
                        u += RSIZE;
                      off[t] = u+RSIZE;
                      stack[s++] = u;
                      t = z;
                    }

                u = stack[--s];
                mycpy(copy,u+depth,rem);
	              while (s > 0)
                  { p = stack[--s];
                    mycpy(u+depth,p+depth,rem);
                    u = p;
                  }
                mycpy(u+depth,copy,rem);
              }
          }
      }
  }
  
  depth += 1;
  if (MARK)
    { if (depth < KSIZE)
        for (y = 0; y < ntop; y++)
          { q = nzero[y];
            n = len[q];
            if (n > S_thr0)
              radix_sort(array, n, depth, alive);
            else if (n > RSIZE)
              shell_sort(array, n, depth);
            *array = 1;
            array += n;
          }
      else
        for (y = 0; y < ntop; y++)
          { q = nzero[y];
            n = len[q];
	          *array = 1;
            array += n;
          }
    }
  else
    { if (depth < KSIZE)
        for (y = 0; y < ntop; y++)
          { q = nzero[y];
            n = len[q];
            if (n > S_thr0)
              radix_sort(array, n, depth, alive);
            else if (n > RSIZE)
              shell_sort(array, n, depth);
            array += n;
          }
    }
}

typedef struct
  { int     tid;
    uint8  *array;
    int64   offs;
    int64   span;
  } Range_Array;

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;

static int  *Tstack;
static int   Tavail;

static void *subsort(void *arg)
{ Range_Array *param = (Range_Array *) arg;

  int      tid   = param->tid;
  uint8   *array = param->array;
  int64    offs  = param->offs;
  int64    span  = param->span;

  int      x;
  int64    alive[256];

  for (x = 0; x < 256; x++)
    alive[x] = 0;

#ifdef SHOW_LAUNCH
  printf("Bucket %d: %12lld - %12lld\n",tid,offs,offs+span); fflush(stdout);
#endif

  radix_sort(array+offs, span, DEPTH, alive);

#ifdef IS_SORTED
  sorted(array+offs,span);
#endif

  if (MARK)
    array[offs] = 1;

  pthread_mutex_lock(&TMUTEX);
    Tstack[Tavail++] = tid;
  pthread_mutex_unlock(&TMUTEX);

  pthread_cond_signal(&TCOND);

  return (NULL);
}


/**********************************************************************************************
* 
*  MSD SORT ENTRANCE, MEMORY_BASED
*
**********************************************************************************************/

void msd_sort(uint8 *array, int64 nels, int rsize, int ksize, int depth, int mark, int nthreads)
{ if (depth >= ksize) return;
  if (nthreads < 1) nthreads = 1;

  pthread_t threads[nthreads];
  int       x;
  int64     part[256];

  NTHREADS = nthreads;   //  global read-only throughout
  RSIZE    = rsize;
  KSIZE    = ksize;
  MARK     = mark;
  DEPTH    = depth;

  THREADS  = threads;
  
  for (x = 0; x < 256; x++)
    part[x] = 0;
  
  if (NTHREADS > 1)
    { if (!DEPTH)
        { KSIZE = 1; // sort by first byte only
          radix_sort(array, nels * RSIZE, 0, part);
          for (x = 0; x < 256; x++)
            part[x] = 0; //restore part array
          KSIZE = ksize; // restore original KSIZE
          DEPTH = 1; // set depth to 1 for next call
        }
      // partition by first byte
      uint64 i;
      uint8 *e = array + DEPTH - 1;
      for (i = 0; i < nels; i++, e += RSIZE)
        part[*e]++;
    }
  else part[0] = nels;
  
  if (KSIZE > DEPTH)
    { Range_Array sparm[nthreads];
      int         stack[nthreads];
      int64       offs;

      S_thr0 = THR0*RSIZE;
      S_thr1 = THR1*RSIZE;
      S_thr2 = THR2*RSIZE;
      S_gap1 = GAP1*RSIZE;
      S_gap2 = GAP2*RSIZE;

      for (x = 0; x < NTHREADS; x++)
        { sparm[x].tid   = x;
          sparm[x].array = array;
        }

      Tstack = stack;
      for (x = 0; x < NTHREADS; x++)
        Tstack[x] = x;
      Tavail = NTHREADS;

      //  For each non-empty part do

      pthread_mutex_init(&TMUTEX,NULL);
      pthread_cond_init(&TCOND,NULL);

      offs = 0;
      for (x = 0; x < 256; x++)
        { int tid;

          if (part[x] == 0)
            continue;

          pthread_mutex_lock(&TMUTEX);

          if (Tavail <= 0)
            pthread_cond_wait(&TCOND,&TMUTEX);

          tid = Tstack[--Tavail];

#ifdef SHOW_LAUNCH
          fprintf(stderr,"Launching part %d on thread %d\n",x,tid);
#endif
    
          pthread_mutex_unlock(&TMUTEX);
    
          sparm[tid].offs  = offs;
          sparm[tid].span  = part[x]*RSIZE;
          offs += sparm[tid].span;
    
          pthread_create(threads+tid,NULL,subsort,sparm+tid);
        }
    
      pthread_mutex_lock(&TMUTEX);
      while (Tavail < NTHREADS)
        pthread_cond_wait(&TCOND,&TMUTEX);
      pthread_mutex_unlock(&TMUTEX);

      if (MARK)
        array[offs] = 1;
    }
}

/*  File: dict.c
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2011
 *-------------------------------------------------------------------
 * Description: based on acedb code from Jean Thierry-Mieg and Richard Durbin 1999-2004
 * -------------------------------------------------------------------
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  1 00:33 2025 (rd109)
 * Created: July 2003 (rd)
 *-------------------------------------------------------------------
 */

#include "dict.h"

#include "array.h"

/****************************************/

static void* remap (void *old, U32 oldSize, U32 newSize)
{
  void* new = mycalloc (newSize, 1) ;
  memcpy (new, old, oldSize) ;
  free (old) ;
  return new ;
}

/****************************************/

static U32 hashString (char *cp, U32 n, bool isDiff)
{
  U32 i ;
  U32 j, x = 0 ;
  U32 rotate = isDiff ? 21 : 13 ;
  U32 leftover = 8 * sizeof(U32) - rotate ;

  while (*cp)
    x = (*cp++) ^ ((x >> leftover) | (x << rotate)) ;

  for (j = x, i = n ; i < sizeof(U32) ; i += n)
    j ^= (x >> i) ;
  j &= (1 << n) - 1 ;

  if (isDiff)
    j |= 1 ;

  return j ;
}

/*****************************/

DICT *dictCreate (U32 size)
{
  DICT *dict = (DICT*) mycalloc (1, sizeof(DICT)) ;

  for (dict->dim = 10, dict->size = 1024 ; dict->size < size ; ++dict->dim, dict->size *= 2) ;
  dict->table = new0 (dict->size, U32) ;
  dict->nameSpace = arrayCreate (((U64)dict->size)<<3, char) ;
  dict->name = new0 (dict->size/2, U32) ;
  return dict ; 
}

/*****************************/

void dictDestroy (DICT *dict)
{
  U32 i ;
  newFree (dict->table, dict->size, U32) ;
  arrayDestroy (dict->nameSpace) ;
  newFree (dict->name, ((U64)dict->size)<<3, char) ;
  free (dict) ;
}

/*****************************/

bool dictWrite (DICT *dict, FILE *f)
{
  if (fwrite (&dict->dim,sizeof(U32),1,f) != 1) return false ;
  if (fwrite (&dict->max,sizeof(U32),1,f) != 1) return false ;
  if (fwrite (dict->table,sizeof(U32),dict->size,f) != dict->size) return false ;
  if (fwrite (dict->name,sizeof(U32),dict->max+1,f) != dict->max+1) return false ;
  if (!arrayWrite (dict->nameSpace, f)) return false ;
  return true ;
}
  
DICT *dictRead (FILE *f)
{
  U32 dim ; if (fread (&dim,sizeof(U32),1,f) != 1) return 0 ;
  DICT *dict = dictCreate (1 << dim) ;
  if (fread (&dict->max,sizeof(U32),1,f) != 1) return 0 ;
  if (fread (dict->table,sizeof(U32),dict->size,f) != dict->size) return 0 ;
  if (fread (dict->name,sizeof(U32),dict->max+1,f) != dict->max+1) return 0 ;
  if (!(dict->nameSpace = arrayRead (f))) return 0 ;
  return dict ;
}

/*****************************/

bool dictFind (DICT *dict, char *s, U32 *ip)
{
  U32 i, x, d ;

  if (!dict) die ("dictFind/Add received null dict\n") ;
  if (!s) die ("dictFind/Add received null string\n") ;

  x = hashString (s, dict->dim, false) ;
  if (!(i = dict->table[x]))
    { dict->newPos = x ; 
      return false ; 
    }
  else if (!strcmp (s, arrp(dict->nameSpace,dict->name[i],char)))
    { if (ip) *ip = i-1 ; 
      return true ; 
    }
  else
    { d = hashString (s, dict->dim, true) ;
      while (1)
	{ x = (x + d) & ((1 << dict->dim) - 1) ;
	  if (!(i = dict->table[x]))
	    { dict->newPos = x ; 
	      return false ; 
	    }
	  else if (!strcmp (s, arrp(dict->nameSpace,dict->name[i],char)))
	    { if (ip) *ip = i-1 ; 
	      return true ; 
	    }
	}
    }
}

/*****************************/

bool dictAdd (DICT *dict, char *s, U32 *ip)
{
  U32 i, x ;

  if (dictFind (dict, s, ip)) return false ;

  i = ++dict->max ;
  dict->table[dict->newPos] = i ;
  dict->name[i] = arrayMax(dict->nameSpace) ;
  array(dict->nameSpace, arrayMax(dict->nameSpace)+strlen(s), char) = 0 ; // terminator for string
  strcpy (arrp(dict->nameSpace, dict->name[i], char), s) ;
  if (ip) *ip = i-1 ;

  if (dict->max > 0.3 * dict->size) /* double table size and remap */
    { U32 *newTable ;
      ++dict->dim ; dict->size *= 2 ;
      dict->name = (U32*) remap (dict->name, (dict->max+1)*sizeof(U32), (dict->size/2)*sizeof(U32)) ;
      newTable = new0 (dict->size, U32) ;
      for (i = 1 ; i <= dict->max ; ++i)
	{ s = arrp(dict->nameSpace, dict->name[i], char) ;
	  x = hashString (s, dict->dim, false) ;
	  if (!newTable[x])
	    newTable[x] = i ;
	  else
	    { U32 d = hashString (s, dict->dim, true) ;
	      while (true)
		{ x = (x + d) & ((1 << dict->dim) - 1) ;
		  if (!newTable[x])
		    { newTable[x] = i ; break ; }
		}
	    }
	}
      newFree (dict->table, dict->size/2, 32) ; dict->table = newTable ;
    }

  return true ;
}

/*********** end of file ***********/

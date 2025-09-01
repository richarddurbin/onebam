/*  File: dict.h
 *  Author: Richard Durbin (rd@sanger.ac.uk)
 *  Copyright (C) Genome Research Limited, 2003-2008
 *-------------------------------------------------------------------
 * Description: header file for DICT package, string hash tables
                developed from the corresponding functions in acedb
		Jean Thierry-Mieg and Richard Durbin 1989-
 * Exported functions:
 * HISTORY:
 * Last edited: Sep  1 00:02 2025 (rd109)
 * * Aug 28 00:23 2025 (rd109): new version using Array for all the names - reduces malloc() usage
 * Created: Sat Dec 20 09:34:14 2008 (rd)
 *-------------------------------------------------------------------
 */

#ifndef DICT_DEFINED
#define DICT_DEFINED

#include "utils.h"
#include "array.h"

typedef struct {
  Array nameSpace ;
  U32  *name ;                  /* index into nameSpace */
  U32  *table ;
  U32   max ;			/* current number of entries */
  U32   dim ;
  U32   size ;			/* 2^dim = size of table */
  U32   newPos ;
} DICT ;

DICT *dictCreate (U32 size) ;
void dictDestroy (DICT *dict) ;
bool dictWrite (DICT *dict, FILE *f) ; /* return success or failure */
DICT *dictRead (FILE *f) ;	       /* return 0 on failure */
bool dictAdd (DICT *dict, char* string, U32 *index) ; /* return TRUE if added, always fill index */
bool dictFind (DICT *dict, char *string, U32 *index) ; /* return TRUE if found */

static inline char* dictName (DICT *dict, U32 i)
{ return arrp(dict->nameSpace, dict->name[i+1], char) ; }

#define dictMax(dict)  ((dict)->max)

#endif

/*********** end of file ***********/

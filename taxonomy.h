/*  File: taxonomy.h
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 15 21:17 2025 (rd109)
 * Created: Tue Oct  7 12:08:08 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "utils.h"
#include "array.h"
#include "hash.h"
#include "dict.h"

typedef I32 TaxID ;     // need to be able to store -1

typedef enum { NONE, ANIMALS, PLANTS, FUNGI, BACTERIA, ARCHAEA } TaxGroup ;

typedef struct {
  int   nTx ;
  Array nodes ;         // of TaxNode - index is the TaxID
  DICT *nameDict ;      // node names
  DICT *rankDict ;
  DICT *divisionDict ;
  DICT *groupDict ;
  TaxID root ;
} Taxonomy ;

typedef struct {
  TaxID    parent ;
  U8       rank ;
  U8       division ;
  U8       geneticCode ;     // we could make dicts for the codes, but for now we don't use them
  U8       mitochondrialCode ;
  TaxID    left, next ;             // for traversal
  TaxID    family ; // corresponding family, or 0 if rank is above
  TaxGroup group ;
} TaxNode ;

// next data structure is working space for finding LCAs - require an object to make threadsafe
typedef struct {
  Taxonomy *tx ;
  TaxID    *mark ;       // will store TaxID or -1, so must be a signed int
  Array     markList ; 
} TaxLCA ;

Taxonomy *taxonomyCreate (void) ;
Taxonomy *taxonomyFromNCBIfiles (const char *path) ; // directory containing nodes.dmp, names.dmp
void      taxonomyDestroy (Taxonomy *tx) ;
TaxLCA   *taxLCAcreate (Taxonomy *tx) ;
TaxID     taxLCAfind (TaxLCA *tl, int n, TaxID *tid) ;
void      taxLCAdestroy (TaxLCA *tl) ;

#include "ONElib.h"

bool      taxonomyWrite (Taxonomy *tx, OneFile *of, bool *txUsed) ;
// if txUsed, then only write that part of the taxonomy, with all paths from it to the root
Taxonomy *taxonomyRead (OneFile *of) ;

/*****************************************/

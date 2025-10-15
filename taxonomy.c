/*  File: taxonomy.c
 *  Author: Richard Durbin (rd109@cam.ac.uk)
 *  Copyright (C) Richard Durbin, Cambridge University, 2025
 *-------------------------------------------------------------------
 * Description: read NCBI taxonomy *.dmp and provide services
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 15 22:07 2025 (rd109)
 * Created: Fri Oct  3 08:00:43 2025 (rd109)
 *-------------------------------------------------------------------
 */

#include "taxonomy.h"

static FILE *openDmpFile (const char *path, char *fileName)
{
  char *name = (char*) malloc (strlen(path) + strlen(fileName) + 6) ;
  strcpy (name, path) ; strcat (name, "/") ; strcat (name, fileName) ; strcat (name, ".dmp") ;
  FILE *f = fopen (name, "r") ;
  free (name) ;
  return f ;
}

static char* *parseLine (FILE *f, int n) // n is expected number of fields
{
  static char buf[1024] ;
  static char *fields[16] ;
  if (n > 16) die ("parseLine %d > 16 in parseLine", n) ;
  int k = 0 ;
  fields[k++] = buf ;
  char *s = buf ;
  while ((*s = getc (f)) != '\n')
    if (*s == EOF) break ;
    else if (*s == '\t')
      { *s++ = 0 ;
	if (getc (f) != '|') die ("file parse error |") ;
	if ((*s = getc (f)) == '\n') break ;
	if (*s != '\t') die ("file parse error \\t") ;
	fields[k++] = s ;
      }
    else ++s ;
  if (k == n) return fields ;
  if (feof(f)) return 0 ;
  die ("read %d fields not %d in parseLine", k, n) ; return 0 ;
}

/* nodes.dmp columns: NCBI documentation preceded by counts in 250926 version
 2700207 tax_id				   -- node id in GenBank taxonomy database                
  250429 parent tax_id			   -- parent node id in GenBank taxonomy database         
      48 rank				   -- rank of this node (domain, kingdom, ...)            
     677 embl code			   -- locus-name prefix; not unique                 IGNORE
      12 division id			   -- see division.dmp file                               
       2 inherited div flag  (1 or 0)	   -- 1 if node inherits division from parent             
      18 genetic code id		   -- see gencode.dmp file                                
       2 inherited GC  flag  (1 or 0)	   -- 1 if node inherits genetic code from parent         
      16 mitochondrial genetic code id	   -- see gencode.dmp file                                
       2 inherited MGC flag  (1 or 0)	   -- 1 if node inherits mitochondrial gencode from parent
       2 GenBank hidden flag (1 or 0)      -- 1 if name is suppressed in GenBank entry lineage    
       1 hidden subtree root flag (1 or 0) -- 1 if this subtree has no sequence data yet          
      10 comments			   -- free-text comments and citations
*/

/* names.dmp columns: NCBI documentation preceded by counts in 250926 version
 2700207 tax_id		-- the id of node associated with this name          
 4463358 name_txt	-- name itself                                        
  339830 unique name	-- the unique variant of this name if name not unique 
      12 name class	-- (synonym, common name, ...)                       
 we will just take the ones for which the name-class == "scientific name" - there is one per tax_id
 we take the unique_name if it exists (only if there is confounding) else name_txt
*/

static void readFiles (Taxonomy *tx, const char *path)
{
  FILE  *f ;
  char **fields ;
  U32    k ;
    
  if (!(f = openDmpFile (path, "nodes"))) die ("failed to open nodes.dmp file at %s", path) ;
  while ((fields = parseLine (f, 13)))
    { TaxID txid = atoi(fields[0]) ;
      TaxNode *node = arrayp(tx->nodes, txid, TaxNode) ;
      node->parent = atoi(fields[1]) ;
      dictAdd (tx->rankDict, fields[2], &k) ; node->rank = k ;
      // ignore fields[3] = embl code
      node->division = atoi(fields[4]) ;
      // ignore fields[5] = inherited division flag
      node->geneticCode = atoi(fields[6]) ;
      // ignore fields[7] = inherited GC flag ;
      node->mitochondrialCode = atoi(fields[8]) ;
      // ignore further fields
      ++tx->nTx ;
    }
  fclose (f) ;
		  
  if (!(f = openDmpFile (path, "names"))) die ("failed to open names.dmp file at %s", path) ;
  char buf[32] ;
  while ((fields = parseLine (f, 4)))
    if (!strcmp (fields[3], "scientific name"))
      { TaxID txid = atoi(fields[0]) ;
	if (txid < dictMax(tx->nameDict)) die ("names out of order") ;
	while (dictMax(tx->nameDict) < txid)
	  { sprintf (buf, "z%d", dictMax(tx->nameDict)) ;
	    dictAdd (tx->nameDict, buf, &k) ;
	  }
	if (*fields[2])
	  dictAdd (tx->nameDict, fields[2], &k) ;
	else
	  dictAdd (tx->nameDict, fields[1], &k) ;
	if (k != txid)
	  die ("mismatch in dictAdd for tx %d node name %s k %d", txid, fields[1], k) ;
      }
  fclose (f) ;
  if (dictMax (tx->nameDict) != arrayMax(tx->nodes)) die ("mismatch nameMax != nodeMax") ;

  if ((f = openDmpFile (path, "division"), 4))
    { tx->divisionDict = dictCreate (16) ;
      while ((fields = parseLine (f, 4))) dictAdd (tx->divisionDict, fields[2], &k) ;
      fclose (f) ;  // code here assumes they come 0, 1, 2.. as they do
    }
}

static void buildTraversal (Taxonomy *tx) // fill left and next
{
  TaxID t ;
  TaxNode *x, *p, *y ;
  TaxID *last = new0 (arrayMax(tx->nodes), TaxID) ;
  for (t = 0 ; t < arrayMax (tx->nodes) ; ++t)
    { x = arrp(tx->nodes, t, TaxNode) ;
      if (!x->parent) continue ; // not a used node
      if (x->parent == t) { tx->root = t ; continue ; }
      p = arrp(tx->nodes, x->parent, TaxNode) ;
      if (!p->left) { p->left = t ; last[t] = t ; continue ; }
      y = arrp (tx->nodes, last[p->left], TaxNode) ;
      y->next = t ;
      last[p->left] = t ;
    }
  newFree (last, arrayMax(tx->nodes), TaxID) ;
}

static inline void groupAssign (Taxonomy *tx, char *name, TaxGroup tg)
{
  U32 t ;
  if (dictFind (tx->nameDict, name, &t)) arrp(tx->nodes, t, TaxNode)->group = tg ;
}

static void annotate (Taxonomy *tx)
{
  TaxID fam = 0 ;
  U32   family ;
  if (!dictFind (tx->rankDict, "family", &family)) die ("failed to find rank family") ;

  TaxGroup group = 0 ;
  tx->groupDict = dictCreate(32) ;
  dictAdd (tx->groupDict, "none", 0) ;     
  dictAdd (tx->groupDict, "animals", 0) ;  groupAssign (tx, "Metazoa", ANIMALS) ;
  dictAdd (tx->groupDict, "plants", 0) ;   groupAssign (tx, "Embryophyta", PLANTS) ;
  dictAdd (tx->groupDict, "fungi", 0) ;    groupAssign (tx, "Fungi", FUNGI) ;
  dictAdd (tx->groupDict, "bacteria", 0) ; groupAssign (tx, "Bacteria <bacteria>", BACTERIA) ;
  dictAdd (tx->groupDict, "archaea", 0) ;  groupAssign (tx, "Archaea", ARCHAEA) ;

  // strategy is to recurse through the tree, left (down) first
  TaxNode *x ;
  TaxID t = tx->root ;
  while (true)
    { x = arrp(tx->nodes, t, TaxNode) ;
      if (x->rank == family) fam = t ;
      x->family = fam ;
      if (!x->group) x->group = arrp(tx->nodes, x->parent, TaxNode)->group ;
      if (x->left) // next level down (left) - inherit current state
	t = x->left ;
      else         // go right (next) or up
	{ if (x->rank == family) fam = 0 ;
	  if (x->next)
	    t = x->next ;
	  else
	    { while (t != tx->root) // at the root
		{ t = x->parent ;
		  x = arrp(tx->nodes, t, TaxNode) ;
		  if (x->rank == family) fam = 0 ;
		  if (x->next) { t = x->next ; break ; }
		}
	    }
	}
      if (t == tx->root) break ;
    }
}

Taxonomy *taxonomyCreate (void)
{
  Taxonomy *tx = new0 (1, Taxonomy) ;
  tx->nodes    = arrayCreate (1<<20, TaxNode) ;
  tx->nameDict = dictCreate (1<<20) ;
  tx->rankDict = dictCreate (256) ;
  return tx ;
}

Taxonomy *taxonomyFromNCBIfiles (const char *path)
{
  Taxonomy *tx = taxonomyCreate () ;
  readFiles (tx, path) ;
  buildTraversal (tx) ;
  annotate (tx) ;
  return tx ;
}

void taxonomyDestroy (Taxonomy *tx)
{
  dictDestroy (tx->nameDict) ;
  dictDestroy (tx->rankDict) ;
  if (tx->divisionDict) dictDestroy (tx->divisionDict) ;
  arrayDestroy (tx->nodes) ;
  newFree (tx, 1, TaxNode) ;
}

TaxLCA *taxLCAcreate (Taxonomy *tx)
{
  TaxLCA *tl = new (1, TaxLCA) ;
  tl->tx = tx ;
  tl->mark = new0 (arrayMax(tx->nodes), I32) ;
  tl->markList = arrayCreate (1024, TaxID) ;
  return tl ;
}

void taxLCAdestroy (TaxLCA *tl)
{ newFree (tl->mark, arrayMax(tl->tx->nodes), I32) ;
  arrayDestroy (tl->markList) ;
  newFree (tl, 1, TaxLCA) ;
}

TaxID taxLCAfind (TaxLCA *tl, int n, TaxID *tid)
{
  if (n && *tid == 0) { --n ; ++tid ; } // ignore 0 which is for unidentified TaxID
  if (!n) return 0 ;
  if (n == 1) return *tid ;

  // Strategy is to mark the tree below and including the current LCA with -1,
  // and the path above the current LCA up to the root with the node below.
  // When adding a leaf, mark with -1 up until you hit a mark.  If -1 then the LCA doesn't change.
  // If positive, then change the LCA and mark back down the path until you hit the old LCA.
  arrayMax (tl->markList) = 0 ;
  TaxID t, lca = *tid, last = -1 ;
  if (lca < 0) printf (" ZZ n %d tid %llx *tid %d\n", n, (long long) tid, *tid) ;
  for (t = *tid ; !tl->mark[t] ; t = arrp(tl->tx->nodes,t,TaxNode)->parent) // root is its own parent
    { tl->mark[t] = last ; array(tl->markList, arrayMax(tl->markList), TaxID) = t ; last = t ; }
  int i ;
  for (i = 1 ; i < n ; ++i)
    { last = -1 ;
      for (t = tid[i] ; !tl->mark[t] ; t = arrp(tl->tx->nodes,t,TaxNode)->parent)
	{ tl->mark[t] = -1 ; array(tl->markList, arrayMax(tl->markList), TaxID) = t ; }
      if (tl->mark[t] > 0) // hit above the current LCA
	{ lca = t ;
	  while (tl->mark[t] > 0) { last = tl->mark[t] ; tl->mark[t] = -1 ; t = last ; }
	}
    }
  
  // clean up the mark space
  for (i = 0 ; i < arrayMax(tl->markList) ; ++i) tl->mark[arr(tl->markList, i, TaxID)] = 0 ;

  return lca ;
}

/**************** read/write to OneFiles  *****************/

// for now ignore division, geneticCode, mitochondrialCode

bool taxonomyWrite (Taxonomy *tx, OneFile *of, bool *txUsed)
// if txUsed, then only write that part of the taxonomy, with all paths from it to the root
{
  if (!of || !oneFileCheckSchemaText (of,
			    "P 3 seq\n"   // would prefer no required primary here
			    "O U 4 3 INT 3 INT 3 INT 6 STRING\n" // taxid, rank, parent, name
			    "O V 1 6 STRING\n"))                 // taxonomy rank
    return false ;

  TaxID    i, t ;
  TaxNode *x ;

  // first build the list of taxids to write - need parental path to the root
  bool *txInclude = new0 (arrayMax(tx->nodes), bool) ;
  for (i = 0 ; i < arrayMax(tx->nodes) ; ++i)
    if (txUsed[i])
      for (t = i ; !txInclude[t] ; t = x->parent) // will break at root whose parent is self
	{ txInclude[t] = true ;
	  x = arrp(tx->nodes,t,TaxNode) ;
	}

  // now write the data for these taxids
  for (i = 0 ; i < arrayMax(tx->nodes) ; ++i)
    if (txInclude[i])
      { oneInt(of,0) = i ;
	TaxNode *x = arrp(tx->nodes,i,TaxNode) ;
	oneInt(of,1) = x->rank ;
	oneInt(of,2) = x->parent ;
	oneWriteLine (of, 'U', strlen(dictName(tx->nameDict,i)), dictName(tx->nameDict,i)) ;
      }
  // and also the rank names
  for (i = 0 ; i < dictMax(tx->rankDict) ; ++i)
    oneWriteLine (of, 'V', strlen(dictName(tx->rankDict,i)), dictName(tx->rankDict,i)) ;
  
  return true ;
}

Taxonomy *taxonomyRead (OneFile *of)
{
  if (!of || !oneFileCheckSchemaText (of,
		      "P 3 seq\n"   // would prefer no required primary here
		      "O U 4 3 INT 3 INT 3 INT 6 STRING\n" // taxid, rank, parent, name
		      "O V 1 6 STRING\n")                  // taxonomy rank
      || !oneGoto (of, 'U', 1)) return 0 ;
  Taxonomy *tx = taxonomyCreate() ;

  char buf[64] ;
  while (oneReadLine (of) && of->lineType == 'U')
    { TaxID    txid = oneInt(of,0) ;
      TaxNode *node = arrayp(tx->nodes, txid, TaxNode) ;
      node->rank    = oneInt(of,1) ;
      node->parent  = oneInt(of,2) ;
      if (txid < dictMax(tx->nameDict)) die ("names out of order") ;
      U32 k ;
      while (dictMax(tx->nameDict) < txid)
	{ sprintf (buf, "z%d", dictMax(tx->nameDict)) ;
	  dictAdd (tx->nameDict, buf, &k) ;
	}
      dictAdd (tx->nameDict, oneString(of), &k) ;
      if (k != txid) die ("mismatch for tx %d node name %s k %d", txid, oneString(of), k) ;
    }
  if (oneGoto (of, 'V', 1))
    while (oneReadLine (of) && of->lineType == 'V')
      dictAdd (tx->rankDict, oneString (of), 0) ;
  
  buildTraversal (tx) ;
  annotate (tx) ;

  return tx ;
}

/****************************************/

#ifdef TAX_TEST

// compile with: gcc -DTAX_TEST -o taxtest taxonomy.c ONElib.c array.c dict.c utils.c -lz

int main (int argc, char *argv[])
{
  --argc ; ++argv ;
  timeUpdate (0) ;
  if (argc < 1) die ("usage: taxtest <dir containing [nodes,names].dmp> <taxid>+") ;
  Taxonomy *tx = taxonomyFromNCBIfiles (*argv) ;
  printf ("read taxonomy with %d nodes\n", tx->nTx) ;
  timeUpdate (stdout) ;

  --argc ; ++argv ;
  if (argc)
    { TaxID *tid = new (argc, TaxID) ;
      int i ;
      for (i = 0 ; i < argc ; ++i)
	{ tid[i] = atoi(argv[i]) ;
	  TaxNode *x = arrp (tx->nodes, tid[i], TaxNode) ;
	  printf ("%u %s %s\n", tid[i], dictName(tx->nameDict, tid[i]),
		  dictName(tx->rankDict, x->rank)) ;
	}
      TaxLCA *tl = taxLCAcreate (tx) ;
      TaxID lca = taxLCAfind (tl, argc, tid) ;
      printf ("  LCA is %u %s %s\n", lca, dictName(tx->nameDict,lca),
	      dictName(tx->rankDict, arrp(tx->nodes, lca, TaxNode)->rank)) ;
      taxLCAdestroy (tl) ;
      newFree (tid, argc, TaxID) ;
    }
  
  taxonomyDestroy (tx) ;
  timeTotal (stdout) ;
}

#endif

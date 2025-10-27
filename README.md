# Onebam
ONEcode package to replace SAM/BAM, especially for eDNA mapping to very large reference databases containing everything.

First, `onebam bam21read` converts [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) files to a [ONEcode](https://github.com/thegenemyers/ONEcode) equivalent with default file ending `.1read`. The corresponding schema is given below. In this file the primary objects are reads, not alignments as in bam files.  For each read that has any matches it contains the name, sequence, qualities and the information required to calculate DNA damage patterns (the best match together with the bases that are different from the reference in the match, and the reference values at those positions), and the information required to assign the read to a taxon (all the taxids to which it has hits, and for each of these the number of hits and the best hit score). We don't keep the alignment coordinates. If you want to study genetic variation in a specific taxon we recommend to remap just the reads assigned to that taxon to its reference genome, and proceed with standard pipelines.

The resulting `.1read` file has general [ONEcode](https://github.com/thegenemyers/ONEcode) properties of being binary and compressed with adaptive data-specific compressors, containing its own indices, supporting multithreading, carrying summary statistics, and enabling generic ascii conversion via `ONEview`. It is sorted on read name, enabling highly efficient merging of `.1read` files with `onebam merge1read`, which makes use of two of Gene Myers' super-efficient [list-merging algorithms](https://drops.dagstuhl.de/storage/00lipics/lipics-vol259-cpm2023/LIPIcs.CPM.2023.22/LIPIcs.CPM.2023.22.pdf).

In order to convert the target identifiers in the BAM file alignment lines to NCBI taxids they need to be matched against entries in a table normally provided as a tab delimited 3 or 4 column file as available from [NCBI](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/). Because this can contain up to billions of lines and take time to parse, we convert it to another ONEcode file of type `.1acctax` with `onebam makeAccTax`.

As well as the `.1read` pacakge, there is also a `onebam bam21bam` command that generates a much more complete onecode representation of a BAM file as a `.1bam` file.  By default the conversion drops the read names, though option `-names` adds them back in.  Also, although `.1bam` supports arbitrary auxiliary information, you have to explicitly specify the tags that you want to keep using option `-aux` with their data types.  (The exception is `RG` for readgroups, which are automatically processed.) The `.1bam` file has general ONEcode file properties but is actually a little larger than the original BAM file in limited testing, and for now this has not been developed further.

## Installation

```
github clone https://github.com/richarddurbin/onebam.git
cd onebam
git submodule update --init --recursive # download required external libraries htslib and zstd
make # -j      # optionally with -j to enable parallel compilation
make install   # by default installs in ~/bin - to put elsewhere "make DESTDIR=<target> install"
```

## Command line arguments/options
The following help information is provided when you run onebam without arguments:

```
Usage: onebam <command> <options>* <args>*
  options and arguments depend on the command
  commands with arguments, each followed by their options:
    bam21read <XX.1acctax> <YY.bam>  convert read-sorted BAM/SAM/CRAM file to .1read
      -o <ZZ.1read>                    output file - default is YY.1read for input YY.bam
    merge1read [*.1read]+         merge .1read files
      -o <ZZ.1read>                    output file - default is 'combined.1read'
    report1read <XX.1read>        report on contents of .1read file // needs more work...
      -o <ZZ.report>                   output to named file rather than stdout
    makeAccTax <XX.tsv>           convert acc2taxid tab-separated file to .1acctax
      -o <ZZ.1acctax>                  output file - default is XX.1acctax for input XX.*
    addLCA <taxdir> <XX.1read>    add LowestCommonAncestor and dust information
      -o <ZZ.1read>                    output file - default is to overwrite input
      -T <nthreads>                    number of threads [8]
    reportLCA <XX.1read>          report aggregated by LCA
      -o <ZZ.tsv>                      tab-separated output file - default is XX.reportLCA
      -rank <taxonomic rank>           collect to and report at named rank, e.g. genus, family
      -group <top level tax group>     restrict to group, eg animals, plants, fungi, bacteria
    extractReads <XX.1read>       extract reads from XX.1read
      -o <ZZ.1read>                    output file - default is XX-<lca>.1read
      -lca <lca>                       extract reads with LCA at or below <lca>
    bamsort [*.bam|*.sam]+        sort BAM/SAM file by query name (ASCII)
      -o <ZZ.bam>                      output file - default is XX.[sorted|merged].bam for input XX.bam
      -T <nthreads>                    number of threads [8]
      -m <maxMem>                      maximum memory per thread; suffix K/M/G recognized [768M]
      -M                               merge only - input files must be presorted
      -zstOutput                       output ZSTD-compressed file
      -noTrimHeader                    do not trim header
```

## Synopsis

### Onebam workflow
```
onebam makeAccTax -o reference.1acctax reference.accession2taxid.gz
for file in *.bam; do                                    #  naturally parallelisable on a cluster
    onebam bam21read reference.1acctax $file             #  makes *.1read
done
onebam merge1read -o merged.1read *.1read                #  merge individual bam files
onebam addLCA ../NCBI-taxonomy/ merged.1read             #  adds the LCA and taxonomy to merged.1read file
onebam reportLCA -rank family merged.1read               #  writes out report collated at family level
sed 's/\t/,/g' merged.reportLCA | sort -t , -k4nr | column -s , -t  #  pretty prints report
onebam extractReads -lca 3688 merged.1read               #  extract reads with LCA Salicacae (taxid 3688)
seqconvert -fa merged-3688.1read                         #  convert reads to fasta (or -fq for fastq)
```

You can view/convert to ASCII the contents of any ONEfile with `ONEview`. Execute this without arguments to see options.  Useful versions are `ONEview -H xx.1read` to see just the header, `ONEview -h xx.1read` to see just the data, and `ONEview -i S 5 xx.1read` to output just the 5th read (reads are sequence objects designated by the `S` onecode).

### BAM/SAM sort

```
onebam bamsort -T 24 -o sorted.bam *.sam *.bam           # sort and merge BAM/SAM files (by read name) 
                                                         # and remove unused references
onebam bamsort -T 24 -M -o merged.bam *.bam *.sam *.zst  # merge sorted BAM/SAM/ZST files 
```


## C interface

There is also a C programming interface, implemented in `oneread.c` with function protypes in `onebam.h` as

```
OneReader *oneReaderCreate (char *fileName, int nThreads) ;
// returns an array of nThreads OneReader objects; just a simple pointer if nThreads = 1

void oneReaderDestroy (OneReader *or) ;
// cleans up memory etc.; call on the base pointer (return value of oneReaderCreate)

bool oneReaderNext (OneReader *or) ;
// moves on to the next read, updating the contents of *or; returns false at end of file
// call on or+i (or equivalently &or[i]) to access the i'th thread if threaded

bool oneReaderGoto (OneReader *or, U64 i) ;
// goes to the requested read
// goto 0 places before the first read, like when you have just started - oneReaderNext() to read it
// goto 1 loads the first read - oneReaderNext() will read the second read, same for following reads
// call on or+i to access the i'th thread if threaded
// returns true if successful, false if i < 0 or i > number of reads

void oneReaderStats (OneReader *or,
		     U64 *nReads,       // number of reads in the file
		     int *maxSeqLen,    // maximum sequence length
		     U64 *totSeqLen,    // total sequence length
		     int *maxNameLen,   // maximum name length
		     int *maxTax        // max value of nTax ;
		     ) ;
// any of the arguments after the first (or) can be 0 (NULL), in which case they won't be filled
```
Then you access the information about the current read by direct inspection of the `*or` object, with structure

```
// the OneReader handle gives access to objects in the file - read-only please
typedef struct {
  // properties of the current read, filled by oneReaderNext()
  int   seqLen ;
  char *seq ;
  U8   *qual ;          // 0-based; length seqLen (NB may contain internal 0s - can use str* funcs)
  int   nameLen ;
  char *name ;
  char *namePrefix ;    // initial segment of the name that is shared by all reads in the file
  char *nameEnd ;       // the segment of the current name after the prefix
  int   maxScore ;      // the best score (AS field) of any alignment
  float dustScore ;     // the dust score of this sequence, in range 0 to 100 (max dusty)
  char *mLine ;         // string containing reference substitutions for an alignment with maxScore
                        //     length seqLen ; positions matching seq are '.', deletions '-'
  int   nTax ;          // how many taxids this read has alignments to
  int  *taxid ;         // list of nTax taxids that have hits from this read
  int  *taxCount ;      // list of nTax counts of how many hits to the respective taxid
  int  *taxBestScore ;  // list of nTax best scores for the respective taxid
  int   lca ;           // NCBI taxid of least common ancestor of taxids hit, after addLCA
  Taxonomy *taxonomy ;  // taxonomy subset for taxids in this file, after addLCA
  
  void *private ;       // a handle for private information for the oneRead code
} OneReader ;
```

Usage is illustrated by the two examples below that are also provided in `oneread.c` and can be compiled as indicated. For compilation you also need `ONElib.h, merge.h, array.h, hash.h` and `utils.h` in the same directory or include path.

```
#ifdef READER_TEST1 // simple example

// compile with: gcc -DREADER_TEST1 -o reader oneread.c ONElib.c merge.c utils.c -lz

int main (int argc, char *argv[])
{
  --argc ; ++argv ; // skip program name
  if (argc != 2) die ("usage: readReport <.1read file>  <i> // reports on i'th read in file") ;
  OneReader *or = oneReaderCreate (*argv, 1) ;
  if (!or) die ("failed to open .1read file %s", *argv) ;
  U64 i = atoi(*++argv) ;
  if (i <= 0) die ("second argument %s must be a positive integer", *argv) ;
  if (!oneReaderGoto (or, i)) die ("failed to locate to %llu", (long long unsigned)i) ;
  printf ("read %llu %s seqLen %d seq %s\n", (long long unsigned)i, or->name, or->seqLen, or->seq) ;
  printf ("maxScore %d dustScore %.1f mLine %s\n", or->maxScore, or->dustScore, or->mLine) ;
  printf   ("nTax %3d taxid   count  bestScore\n", or->nTax) ;
  int j ;
  for (j = 0 ; j < or->nTax ; ++j)
    printf ("        %8d   %5d  %8d\n", or->taxid[j], or->taxCount[j], or->taxBestScore[j]) ;
  oneReaderDestroy (or) ;
}

#endif

#ifdef READER_TEST2 // example using threading

// compile with: gcc -DREADER_TEST2 -o dustbin oneread.c ONElib.c merge.c utils.c -lz

typedef struct {
  OneReader *or ;
  U64        start, end ;
  U64        dustBin[101] ;
} ThreadArg ;

static void *threadProcess (void *arg)
{
  ThreadArg *ta = (ThreadArg*)arg ;
  if (!oneReaderGoto (ta->or, ta->start)) die ("failed oneReaderGoto") ;
  while (ta->start++ < ta->end && oneReaderNext (ta->or))
    ++ta->dustBin[(int)(ta->or->dustScore)] ;
  return 0 ;
}

int main (int argc, char *argv[])
{
  int i, j, n ;
  --argc ; ++argv ; // skip program name
  if (argc != 2)
    die ("usage: dustBin <.1read file> <nThread> // report distribution of dust values") ;
  int nThread = atoi(argv[1]) ;
  if (nThread <= 0) die ("second argument %s must be a positive integer", argv[1]) ;
  OneReader *or = oneReaderCreate (argv[0], nThread) ;
  if (!or) die ("failed to open .1read file %s", argv[0]) ;

  pthread_t *threads = new (nThread, pthread_t) ;
  ThreadArg *ta      = new0 (nThread, ThreadArg) ;
  U64 nRead ; oneReaderStats (or, &nRead, 0, 0, 0, 0) ;
  if (!nRead) die ("no reads found in file %s", *argv) ;
  for (i = 0 ; i < nThread ; ++i)
    { ta[i].or  = or + i ;                     // equivalent to ta[i].or = &or[i]
      ta[i].start = (nRead * i) / nThread ;
      ta[i].end = (nRead * (i+1)) / nThread ;
    }
  for (i = 0 ; i < nThread ; ++i)              // create threads
    pthread_create (&threads[i], 0, threadProcess, &ta[i]) ;
  for (i = 0 ; i < nThread ; ++i)
    pthread_join (threads[i], 0) ;             // wait for threads to complete
  newFree (threads, nThread, pthread_t) ;

  // accumulate scores into the 0 level
  for (i = 1 ; i < nThread ; ++i) for (j = 0 ; j < 101 ; ++j) ta->dustBin[j] += ta[i].dustBin[j] ;

  // find the largest index in the histogram table with data
  for (n = 101 ; n-- ;) if (ta->dustBin[n]) break ;

  printf ("file %s containing %llu reads has dust score distribution:\n", *argv, nRead) ;
  for (j = 0 ; j <= n ; ++j)
    printf ("%d\t%.1f%%\t%llu\n", j, (100.0*ta->dustBin[j])/nRead, ta->dustBin[j]) ;

  newFree (ta, nThread, ThreadArg) ;
  oneReaderDestroy (or) ;
}
#endif
```

---

Richard Durbin 15/10/2025

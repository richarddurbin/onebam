# onebam
ONEcode package to replace SAM/BAM, especially for eDNA mapping to everything.

First, `onebam bam21bam` converts [SAM/BAM format](https://samtools.github.io/hts-specs/SAMv1.pdf) files to a [ONEcode](https://github.com/thegenemyers/ONEcode) equivalent with default file ending `XXX.1bam`. The corresponding schema is given below.  By default the conversion drops the read names, though option `-names` adds them back in.  Also, although .1bam supports arbitrary auxiliary information, you have to explicitly specify the tags that you want to keep using option `-aux` along with their ONEcode characters (by default lower case letters, limiting to 26 auxiliary tags) and their data types.  (The exception is `RG` for readgroups, which is automatically processed.)

The resulting .1bam file has general [ONEcode](https://github.com/thegenemyers/ONEcode) properties of containing its own indices, supporting multithreading, carrying summary statistics, and enabling generic ascii conversion via `oneview`.  However, in initial tests it is larger than the bgzip compressed bam file.

The idea of `onebam makebin` is to convert a BAM file into two reduced representation fixed record length binary files, each of which carries key information for environmental DNA analysis.  These can be efficiently merged by concatenation then sorting with a radix sort, hence supporting integration of the outputs of mapping an eDNA read set against a multi-terabase reference of everything that is split into shards. A synposis is given below.   

###Installation

First you need to install htslib in a directory adjacent to the one in which you will install onebam:

```
github clone https://github.com/samtools/htslib.git
cd htslib
autoreconf -i  # Build the configure script and install files it uses
./configure    # Optional but recommended, for choosing extra functionality
make
make install
cd ..
```

Now onebam:

```
github clone https://github.com/richarddurbin/onebam.git
cd onebam
make
make install
```

###command line arguments/options
The following help information is provided when you run onebam without arguments:

```
Usage: onebam <command> <options>* <args>*
  options and arguments depend on the command
  commands with arguments, each followed by their options:
    bam21bam <XX.bam>             convert BAM/SAM/CRAM file to .1bam
      -o <ZZ.1bam>                     output file - default is XX.1bam for input XX.bam
      -taxid <XX.tsv>                 file with lines <acc>\ttaxid\n sorted on acc
         NB you must give taxids if you plan to use your .1bam to make .1read
      -aux <tag>:<fmt> <char>          record BAM tag with OneCode <char>, e.g. -aux AS:i s
         NB you must use '-aux AS:i s -aux MD:Z m' if you plant to use your .1bam to make .1read
      -names                           keep sequence names [default to drop names]
      -cramref <file_name|URL>         reference for CRAM - needed to read cram
    make1read <XX.1bam>           convert .1bam file to .1read (simplified information per read)
      -T <nthreads>                    number of threads [8]
      -o <ZZ.1read>                    output file - default is XX.1read for input XX.1bam
    makebin <taxid.tsv> <XX.bam>  make fixed-width binary .alb and .txb files from bam file
      -oTxb <ZZ.txb>                   binary taxid output file - default XX.txb from XX.bam
      -oAlb <ZZ.alb>                   binary alignment output file - default XX.alb from XX.bam
      -prefixLen <n>                   ignore first n characters in names [0]
      -maxChars <n>                    write the next n chars of the name, or all if fewer [48]
      -maxEdit <n>                     maximum number of edits for read to be in .alb file [4]
    bin21read <YY.1seq> <XX.bin>  make .1read file from sorted .bin and .1seq
      -T <nthreads>                    number of threads [8]
      -o <ZZ.1read>                    output file - default is XX.1read for input XX.bin

```

###makebin
Synopsis

```
> seqstat -i map_output.bam
bam file, 16947064722 sequences >= 0, 1009029878567 total, 59.54 average, 30 min, 171 max
identifier max_length 40 shared_prefix A00706:796:HTMGLDSX5: shared_prefix_length 21
> onebam makebin -prefixLen 21 -maxChars 19 taxids.tsv 7724-full.bam
```

Typically all the reads in a readset share a common prefix.  When running the makebin option you should specify the length of this using option `-prefixLen`.  You also need an upper bound on the number of additional characters in the longest read name (option `-maxchars`).  These 

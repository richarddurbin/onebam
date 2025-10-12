# makefile for onebam

CFLAGS = -O3
#CFLAGS = -g	# for debugging

LIBS = -lpthread

ALL = onebam ONEview ONEstat seqstat

DESTDIR = ~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	$(RM) -r *.dSYM

### object files

UTILS_OBJS = utils.o array.o hash.o merge.o dict.o
UTILS_HEADERS = utils.h array.h hash.h merge.h dict.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

HTS_DIR = $(PWD)/../htslib/.
HTS_OPTS = -I$(HTS_DIR)/htslib/
HTS_LIBS = -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz 

SEQIO_OPTS = -DONEIO -DBAMIO $(HTS_OPTS)
#SEQIO_LIBS = -lm -lz

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

onebamhts.o: onebamhts.c onebam.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

oneread.o: oneread.c onebam.h merge.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

onebamtax.o: onebamtax.c onebam.h taxonomy.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

seqio.o: seqio.c seqio.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

taxonomy.o: taxonomy.c taxonomy.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

MSDsort.o: MSDsort.c gene_core.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

gene_core.o: gene_core.c gene_core.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

### programs

onebam: onebam.c onebamhts.o oneread.o onebamtax.o taxonomy.o seqio.o ONElib.o MSDsort.o gene_core.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) $(HTS_LIBS) $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

ONEstat: ONEstat.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

seqstat: seqstat.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) $(HTS_LIBS) $(LIBS)

### test

test: onebam
	./onebam

### end of file

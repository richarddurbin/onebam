# makefile for onebam

CFLAGS = -O3
#CFLAGS = -g	# for debugging

LIBS = -lpthread

ALL = onebam ONEview ONEstat seqstat albview txbview

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
ifeq ($(origin HTSSRC), undefined) # Previous behaviour
  HTS_OPTS = -I$(HTS_DIR)/htslib/
  HTS_LIBS = -L$(HTS_DIR) -Wl,-rpath $(HTS_DIR) -lhts -lm -lbz2 -llzma -lcurl -lz 
else # htslib from another location
  HTS_OPTS = -I$(HTSSRC)/include/htslib
  HTS_LIBS = -L$(HTSSRC)/lib -Wl,-rpath $(HTSSRC)/lib -lhts -lm -lbz2 -llzma -lcurl -lz
endif

SEQIO_OPTS = -DONEIO -DBAMIO $(HTS_OPTS)
#SEQIO_LIBS = -lm -lz

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

onebamhts.o: onebamhts.c onebam.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

oneread.o: oneread.c onebam.h merge.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

albcode.o: albcode.c onebam.h
	$(CC) $(CFLAGS) $(HTS_OPTS) -c $^

seqio.o: seqio.c seqio.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

### programs

onebam: onebam.c onebamhts.o oneread.o seqio.o albcode.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) $(HTS_LIBS) $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

ONEstat: ONEstat.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

seqstat: seqstat.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(SEQIO_LIBS) $(HTS_LIBS) $(LIBS)

albview: albview.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lz

txbview: txbview.c $(UTILS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ -lz

### test

test: onebam
	./onebam

### end of file

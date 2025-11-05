# makefile for onebam

CFLAGS = -O3
#CFLAGS = -g	# for debugging

LIBS = -lpthread -lm -lbz2 -llzma -lcurl -lz

ALL = onebam ONEview ONEstat seqstat seqconvert

DESTDIR = ~/bin

# external libraries
HTS_DIR = external/htslib
ZSTD_DIR = external/zstd
LIB_DIR = lib

.PHONY: ${LIB_DIR}

LIBHTS_SOVERSION = 3
LIBZSTD_SOVERSION = 1

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SHARED_EXT = dylib
    LIBHTS_SO_EXT = $(LIBHTS_SOVERSION).dylib
	LIBZSTD_SO_EXT = $(LIBZSTD_SOVERSION).dylib
	R_PATH = -L$(LIB_DIR) -Wl,-rpath,@executable_path/$(LIB_DIR)
else
	SHARED_EXT = so
	LIBHTS_SO_EXT = so.$(LIBHTS_SOVERSION)
	LIBZSTD_SO_EXT = so.$(LIBZSTD_SOVERSION)
    R_PATH = -L$(LIB_DIR) -Wl,-rpath,'$$ORIGIN/$(LIB_DIR)'
endif

HTS_LIB = $(HTS_DIR)/libhts.$(SHARED_EXT)
ZSTD_LIB = $(ZSTD_DIR)/lib/libzstd.$(SHARED_EXT)

all: $(HTS_LIB) $(ZSTD_LIB) $(LIB_DIR) $(ALL)

install:
	cp $(ALL) $(DESTDIR)
	cp -Paf $(LIB_DIR) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	$(RM) -r *.dSYM

### object files
$(HTS_LIB):
	cd $(HTS_DIR) && autoreconf -i && ./configure && $(MAKE)
HTS_OPTS = -I$(HTS_DIR)/htslib/
HTS_LIBS = -lhts

$(ZSTD_LIB):
	cd $(ZSTD_DIR) && $(MAKE) lib-mt
ZSTD_OPTS = -I$(ZSTD_DIR)/lib/
ZSTD_LIBS = -lzstd

SEQIO_OPTS = -DONEIO -DBAMIO $(HTS_OPTS)

UTILS_OBJS = utils.o array.o hash.o merge.o dict.o
UTILS_HEADERS = utils.h array.h hash.h merge.h dict.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

# copy libraries to lib directory
$(LIB_DIR)/libhts.$(SHARED_EXT): $(HTS_LIB)
	@mkdir -p $(LIB_DIR)
	@cp -f $< $@
ifeq ($(UNAME_S),Darwin)
	@install_name_tool -id @rpath/libhts.$(SHARED_EXT) $@
endif

$(LIB_DIR)/libhts.$(LIBHTS_SO_EXT): $(LIB_DIR)/libhts.$(SHARED_EXT)
	@ln -sf libhts.$(SHARED_EXT) $@

$(LIB_DIR)/libzstd.$(SHARED_EXT): $(ZSTD_LIB)
	@mkdir -p $(LIB_DIR)
	@cp -f $< $@
ifeq ($(UNAME_S),Darwin)
	@install_name_tool -id @rpath/libzstd.$(SHARED_EXT) $@
endif

$(LIB_DIR)/libzstd.$(LIBZSTD_SO_EXT): $(LIB_DIR)/libzstd.$(SHARED_EXT)
	@ln -sf libzstd.$(SHARED_EXT) $@

HTS_OBJS = $(LIB_DIR)/libhts.$(SHARED_EXT) $(LIB_DIR)/libhts.$(LIBHTS_SO_EXT)
ZSTD_OBJS = $(LIB_DIR)/libzstd.$(SHARED_EXT) $(LIB_DIR)/libzstd.$(LIBZSTD_SO_EXT)

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

bamsort.o: bamsort.c onebam.h
	$(CC) $(CFLAGS) $(HTS_OPTS) $(ZSTD_OPTS) -c $^

MSDsort.o: MSDsort.c gene_core.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

gene_core.o: gene_core.c gene_core.h
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

### programs

onebam: onebam.c onebamhts.o oneread.o onebamtax.o taxonomy.o bamsort.o seqio.o ONElib.o MSDsort.o gene_core.o $(UTILS_OBJS) | $(HTS_OBJS) $(ZSTD_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(R_PATH) $(HTS_LIBS) $(ZSTD_LIBS) $(LIBS)

ONEview: ONEview.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

ONEstat: ONEstat.c ONElib.o
	$(CC) $(CFLAGS) -o $@ $^

seqstat: seqstat.c seqio.o ONElib.o $(UTILS_OBJS) | $(HTS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(R_PATH) $(HTS_LIBS) $(LIBS)

seqconvert: seqconvert.c seqio.o ONElib.o $(UTILS_OBJS) | $(HTS_OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(R_PATH) $(HTS_LIBS) $(LIBS)

### test

test: onebam
	./onebam

### end of file

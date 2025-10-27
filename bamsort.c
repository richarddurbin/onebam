/*  File: bamsort.c
 *  Author: Chenxi Zhou (cz370@cam.ac.uk)
 *  Copyright (C) University of Cambridge, 2025
 *-------------------------------------------------------------------
 * Description:
 * Exported functions:
 * HISTORY:
 * Last edited: Oct 22 11:20:12 2025 (cz370)
 * Created: Wed Sep  5 22:13:53 2025 (cz370)
 *-------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <zlib.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/resource.h>

#include "zstd.h"
#include "bgzf.h"
#include "sam.h"
#include "hfile.h"
#include "kseq.h"
#include "khash.h"
#include "kstring.h"
#include "thread_pool.h"

#include "onebam.h"
#include "merge.h"

#define PROGRAM "bamsort"

static void infomsg(char *format, ...);
static void errmsg(char *format, ...);
static void warnmsg(char *format, ...);
static void dbgmsg(char *format, ...);

typedef struct timer1_t {
    struct timeval tv;
    struct timer1_t *next;
} timer1_t;
timer1_t *TIME_STACK = NULL;
void push_time(timer1_t **timer);
double pop_time(timer1_t **timer);
double elapsed_time(timer1_t **timer, bool reset);
void free_timer(timer1_t *timer);

#define array_t(type) struct { size_t n, m; type *a; }
#define array_push(type, v, x) do { \
    if ((v).n == (v).m) { \
        (v).m = (v).m? (v).m<<1 : 16; \
        (v).a = (type*)realloc((v).a, sizeof(type) * (v).m); \
    } \
    (v).a[(v).n++] = (x); \
} while (0)

static uint64_t roundup64(uint64_t x) 
{
    if (x == 0) return 1;
    x--; x |= x >> 1; x |= x >> 2; x |= x >> 4; x |= x >> 8; x |= x >> 16; x |= x >> 32; x++;
    return x;
}

static kstring_t *ks_new(size_t size)
{
    kstring_t *ks = (kstring_t *) calloc(1, sizeof(kstring_t));
    if (!ks) return NULL;
    ks->s = (char *) malloc(size);
    if (ks->s) ks->m = size;
    return ks;
}

static inline char *bs_strchrnul(const char *s, const char c) {
    while (*s && *s != c) 
        ++s;
    return (char *) s;
}

KHASH_MAP_INIT_STR(c2i, int)        // ref name -> new tid
KHASH_MAP_INIT_STR(c2c, char*)      // old PG ID -> new PG ID
KHASH_INIT(cset, char*, char, 0, kh_str_hash_func, kh_str_hash_equal) // string set

#define CSTR_SIZE (1ULL<<24) // 16MB
typedef struct {
    char   *string;
    size_t  length;
} cstr_t;
typedef array_t(cstr_t) cstr_array_t;
typedef array_t(int)    cint_array_t;

static void cstr_array_destroy(cstr_array_t *arr);
static char *cstr_array_putsn(const char *s, int l, cstr_array_t *arr);

typedef struct {
    kstring_t    hd_text;           // @HD line
    cstr_array_t sq_text;           // @SQ lines (null-terminated for @SN)
    cstr_array_t rg_text;           // @RG lines
    cstr_array_t pg_text;           // @PG lines
    cstr_array_t co_text;           // @CO lines (verbatim)
    int n_targets;
    char **target_name;             // point to sq_text
    uint32_t *target_len;
    khash_t(c2i)  *sq_map;         // name -> tid point to sq_text
    khash_t(cset) *rg_ids;          // @RG-IDs
    khash_t(cset) *pg_ids;          // @PG-IDs
} sam_hdr_lite_t;

static sam_hdr_lite_t *sam_hdr_lite_init();
static sam_hdr_lite_t *sam_hdr_lite_builder_init();
static sam_hdr_lite_t *sam_hdr_lite_read(samFile *fp);
static void sam_hdr_lite_destroy(sam_hdr_lite_t *h);
static void sam_hdr_lite_free_text(sam_hdr_lite_t *hdr);
static bool sam_hdr_lite_build_tbl(sam_hdr_lite_t *bdr, sam_hdr_lite_t *hdr, int *refmap);
static int sam_parse1_lite(kstring_t *s, sam_hdr_lite_t *h, bam1_t *b);
static int sam_format1_lite(const sam_hdr_lite_t *h, const bam1_t *b, kstring_t *str);
static size_t sam_hdr_lite_min_text_len(sam_hdr_lite_t *hdr);
static bool sam_hdr_lite_build_targets_tbl(sam_hdr_lite_t *hdr);
static bool sam_hdr_lite_build_targets_map(sam_hdr_lite_t *hdr);
static ssize_t sam_hdr_lite_add_lines(sam_hdr_lite_t *hdr, const char *lines, size_t len);
static int sam_hdr_lite_build_HD(sam_hdr_lite_t *hdr, const char *vn, const char *so);
static inline int sam_read1_lite(samFile *fp, sam_hdr_lite_t *header, bam1_t *b);
static inline int bam_read1_lite(samFile *fp, sam_hdr_lite_t *header, bam1_t *b);

typedef struct {
    // this is used for header deduplication
    // if multiple headers share the same reference set
    // they will share the reference map to save memory
    int nref;
    int **refs;                     // reference maps for each (unique) input header
                                    // will not be expanded so enough space should be allocated
    FILE *sfile;                    // file for (unique) original @SQ content - currently unzipped
    char *fbuff;                    // buffer for reading from sfile
    size_t bsize;                   // allocated size of fbuff
    fpos_t *sfpos;                  // file positions of original @SQ lines - size of headers
    khash_t(c2i)  *shash;           // header fingerprint -> leading + trailing + string length
    cint_array_t  *index;           // index to the arrays - more than one when hash collision
} refseq_builder_t;

static refseq_builder_t *refseq_builder_init(int capacity, size_t buff_size, char *fn_template);
static int refseq_builder_add(refseq_builder_t *rdb, sam_hdr_lite_t *hdr, int **_refmap);
static void refseq_builder_destroy(refseq_builder_t *rdb);

#define ZSTD_BAM_MAGIC 0x5A53544442414D01 // "ZSTDBAM\1"
#define ZSTD_COMPRESSION_LEVEL 1
#define ZSTD_MAX_THREADS 64
#define ZSTD_INBUF_SIZE 2*1024*1024ULL
#define ZSTD_OUTBUF_SIZE 16*1024*1024ULL
#define ZSTD_MAXBUF_SIZE 128*1024*1024ULL
typedef struct {
    ZSTD_CCtx *cctx;
    size_t inBufSize, outBufSize;
    char *inBuf, *outBuf;
} zstd_cctx_t;
typedef struct {
    ZSTD_DCtx *dctx;
    ZSTD_inBuffer inBuf; // will be static
    size_t outBufSize;
    char *outBuf;
} zstd_dctx_t;
typedef struct {
    zstd_cctx_t *zstd_cctx;
    ZSTD_CCtx *cctx;
    size_t inBufSize, outBufSize;
    char *inBuf, *outBuf;
    size_t offset;
    FILE *fp;
} zstd_writer_t;

static zstd_cctx_t *zstd_cctx_init(int compressionLevel, int threads);
static void zstd_cctx_destroy(zstd_cctx_t *zstd_cctx);
static zstd_dctx_t *zstd_dctx_init(size_t inBufSize, size_t outBufSize);
static void zstd_dctx_destroy(zstd_dctx_t *zstd_dctx);
static zstd_writer_t *zstd_writer_init(int compressionLevel, int threads, FILE *fp);
static inline int zstd_writer_flush(zstd_writer_t *writer);
static void zstd_writer_destroy(zstd_writer_t *writer);
static inline int zstd_bam_write1(zstd_writer_t *writer, bam1_t *b);
static sam_hdr_lite_t *zstd_sam_hdr_lite_read(FILE *fp);
static int zstd_sam_hdr_lite_write(zstd_writer_t *writer, sam_hdr_lite_t *hdr);

static int IS_BIG_ENDIAN = 0;
const size_t BAM_REC_SIZE = sizeof(bam1_t *) + 300;
const size_t WW_MIN_MEM = BGZF_MAX_BLOCK_SIZE;
const size_t WW_MAX_MEM = 128ULL<<20; // 128 MB
const size_t WW_MIN_REC = WW_MIN_MEM / BAM_REC_SIZE;
const size_t WW_MAX_REC = WW_MAX_MEM / BAM_REC_SIZE;

#define BAM_MAX_HEADER_SIZE UINT32_MAX
#define MAX_WWORKER_THREADS (1<<6) // 64
#define TINY_QUEUE_MASK (MAX_WWORKER_THREADS - 1)

typedef struct { // a simplified version of kdq
	int front, count;
	int a[MAX_WWORKER_THREADS];
} tiny_queue_t;

typedef struct bam_writer_t bam_writer_t;

typedef struct {
    bam_writer_t *writer;
    int tid;
    int status;
    bam1_t **brec;
    size_t   bcnt, bcap;
    // actual data for BAM records 
    uint8_t *ubuf;
    size_t   ulen, ucap;
    // compressed data buffer
    uint8_t *cbuf;
    size_t   clen, ccap;
    // uncompressed data buffer
    // for serialised BAM records
    // size BGZF_BLOCK_SIZE
    uint8_t *bbuf;
    size_t   blen;
    // function args
    void  *args;
    void (*argfree)(void *args);
} bam_wworker_t;

typedef struct bam_writer_t {
    FILE *fp;
    int n_workers;
    bam_wworker_t *workers;
    // current input worker
    bam_wworker_t *input;
    // queues
    tiny_queue_t input_q;  // indices waiting for INPUT
    tiny_queue_t compr_q;  // indices waiting for COMPRESS and WRITE
    // synchronisation
    pthread_mutex_t mtx;
    pthread_cond_t  cond_input;
    pthread_cond_t  cond_compr;
    pthread_cond_t  cond_write;
    int stop;
    // thread handles
    pthread_t *compr_threads;
    pthread_t  main_thread;
    // compressor
    void (*compress)(bam_wworker_t *bw);
    // function args
    void *args;
    void (*argfree)(void *args);
} bam_writer_t;

static bam_writer_t *bam_writer_init(FILE *fp, int n_workers, size_t ubuf_cap);
static void bam_writer_destroy(bam_writer_t *bw);
static int bam_writer_close(bam_writer_t *bw);
static int bw_bam_hdr_lite_write(FILE *fp, sam_hdr_lite_t *hdr);
static void bw_bam_compress(bam_wworker_t *bw);
static int bw_sam_hdr_lite_write(FILE *fp, sam_hdr_lite_t *hdr);
static void bw_sam_compress(bam_wworker_t *bw);
static inline int bw_bam_write1(bam_writer_t *bw, bam1_t *b);

void *sam_compArgs_make(sam_hdr_lite_t *hdr, size_t buf);
void sam_compArgs_free(void *args);

#undef VALIDATE_SORT
#define USE_RADIXSORT

typedef char * cptr_t;

#if defined USE_MERGESORT
#include "ksort.h"
static inline int cptr_lt(const cptr_t a, const cptr_t b)
{
    return strcmp(a, b) < 0;
}
KSORT_INIT(sort, cptr_t, cptr_lt)
#elif defined USE_QUICKSORT
static int cptr_cmp(const void *a, const void *b)
{
    return strcmp(*((cptr_t *)a), *((cptr_t *)b));
}
#else
static int cptr_cmp(const void *a, const void *b)
{
    return strcmp(*((cptr_t *)a), *((cptr_t *)b));
}

static inline int cptr_lt(const cptr_t a, const cptr_t b)
{
    return strcmp(a, b) < 0;
}

static inline void insert_sort(cptr_t *array, int n)
{
    int i, j;
    cptr_t t;
    for (i = 1; i < n; i++) {
        t = array[i];
        for (j = i-1; j >= 0; j--) {
            if (cptr_lt(array[j], t))
                break;
            else
                array[j+1] = array[j];
        }
        array[j+1] = t;
    }
}

#define ISORT_NMAX 10
#define QSORT_NMAX 100

static void radix_subsort(cptr_t *array, size_t nels, size_t index, cptr_t *swaps)
{
    if (nels <= 1) return;

    size_t i, j, k, n, ntop, alive[256] = {0};
    cptr_t *x, *aend, *bucket[256];
    int c, nzero[256];
    char *name;

    // skip common prefix
    aend = array + nels;
    x = array;
    while (1) {
        c = (*x)[index];
        for (x = x + 1; x < aend; x++) {
            if ((*x)[index] != c)
                break;
        }
        if (x < aend)
            break;
        if (!c) // end of string
            return;
        index += 1;
        x = array;
    }

    ntop = 1;
    nzero[0] = c;
    alive[c] = x - array;

    for (; x < aend; x++) {
        c = (*x)[index];
        if (alive[c] == 0)
            nzero[ntop++] = c;
        alive[c]++;
    }

    x = swaps;
    if (ntop <= ISORT_NMAX) {
        for (i = 1; i < ntop; i++) {
            c = nzero[i];
            for (j = i-1; j != (size_t) -1; j--) {
                if (nzero[j] < c)
                    break;
                else
                    nzero[j+1] = nzero[j];
            }
            nzero[j+1] = c;
        }
        for (i = 0; i < ntop; i++) {
            c = nzero[i];
            bucket[c] = x;
            x += alive[c];
        }
    } else {
        ntop = 0;
        for (c = 0; c < 256; c++) {
            if (alive[c] > 0) {
                nzero[ntop++] = c;
                bucket[c] = x;
                x += alive[c];
            }
        }
    }

    int y;
    x = array;
    for (i = 0; i < ntop; i++) {
        c = nzero[i];
        aend = x + alive[c];
        for (; x < aend; x++) {
            y = (*x)[index];
            if (c != y) {
                *bucket[y]++ = *x;
                *x = NULL;
            }
        }
    }
    x = swaps;
    for (i = 0; i < ntop; i++) {
        c = nzero[i];
        bucket[c] = x;
        x += alive[c];
    }
    cptr_t *b;
    x = array;
    for (i = 0; i < ntop; i++) {
        c = nzero[i];
        aend = x + alive[c];
        b = bucket[c];
        for (; x < aend; x++) {
            if (*x == NULL)
                *x = *b++;
        }
    }

    /***
    for (x = array; x < aend; x++) {
        c = (*x)[index];
        *bucket[c]++ = *x;
    }
    ***/

    //memcpy(array, swaps, sizeof(cptr_t) * nels);

    index += 1;
    i = 0;
    if (!nzero[0]) {
        i += 1;
        array += alive[0];
        swaps += alive[0];
    }
    
    // continue sorting in place
    for (; i < ntop; i++) {
        j = alive[nzero[i]];
        if (j > QSORT_NMAX)
            radix_subsort(array, j, index, swaps);
        else {
            x = array;
            aend = x + j;
            for (; x < aend; x++) (*x) += index;
            if (j > ISORT_NMAX)
                qsort(array, j, sizeof(cptr_t), cptr_cmp);
            else
                insert_sort(array, j);
            x = array;
            for (; x < aend; x++) (*x) -= index;
        }
        array += j;
        swaps += j;
    }

    return;
}

static void quick_subsort(cptr_t *array, size_t nels, size_t index)
{
    cptr_t *x, *aend;

    if (index > 0)
        for (x = array, aend = x + nels; x < aend; x++)
            (*x) += index;
    
    if (nels > ISORT_NMAX) 
        qsort(array, nels, sizeof(cptr_t), cptr_cmp);
    else
        insert_sort(array, nels);
    
    if (index > 0)
        for (x = array, aend = array + nels; x < aend; x++)
            (*x) -= index;

    return;
}

#undef TEST_MSDSORT

#ifdef TEST_MSDSORT
void msd_sort(uint8_t *array, size_t nels, int rsize, int ksize, int depth, int mark, int nthreads);

static char *common_prefix(cptr_t *array, size_t nels, int *_len)
#else
static char *common_prefix(cptr_t *array, size_t nels)
#endif
{

    cptr_t *beg = array;
    cptr_t *end = array + nels;
    char *fst, *p, *q;
#ifdef TEST_MSDSORT
    int l, len;
#endif

    if (beg == end) {
#ifdef TEST_MSDSORT
        if (_len) *_len = 0;
#endif
        return NULL;
    }

    fst = strdup(*beg);
#ifdef TEST_MSDSORT
    len = strlen(fst);
#endif
    for (end = end - 1; end > beg; end--) {
        p =  fst;
        q = *end;
        while (*p && *p == *q) p++, q++;
        *p = 0;
#ifdef TEST_MSDSORT
        l = p - fst + strlen(q);
        if (l > len) len = l;
#endif
    }

#ifdef TEST_MSDSORT
    if (_len) *_len = len;
#endif

    return fst;
}

void hybrid_sort(cptr_t *array, size_t nels, int *_lcp) {
    cptr_t *swaps = NULL;
    size_t Mnels;
    int i, lcp;
    char *prefix;
    
    swaps = (cptr_t *) malloc(nels * sizeof(cptr_t));
    if (swaps == NULL) {
        errmsg("swaps memory allocation failed");
        exit (1);
    }

    // find common prefix length
#ifdef TEST_MSDSORT
    int len = 0;
    prefix = common_prefix(array, nels, &len);
#else
    prefix = common_prefix(array, nels);
#endif
    lcp = prefix? strlen(prefix) : 0;
    free(prefix);
    
    if (_lcp) *_lcp = lcp;

#ifdef TEST_MSDSORT
    /*** MSD radix sort - out of place - not used because it is
     * not faster than current implementation
     * data allocation is quite time consuming - could be multithreaded?
     ***/
    int wide = len - lcp;
    if (wide > 0) {
        uint8_t *data, *p;
        char *s;
        int w;
        wide = (wide + 7) & ~7; // make it 8-byte aligned
        data = (uint8_t *) calloc(nels * (wide + 8), sizeof(uint8_t));
        if (data == NULL) {
            errmsg("sorting memory allocation failed");
            exit (1);
        }

        p = data;
        for (i = 0; i < nels; i++) {
            s = array[i] + lcp;
            w = 0;
            while (*s) *p++ = *s++, w++;
            p += wide - w;
            *(cptr_t *)p = array[i];
            p += 8;
        }

        msd_sort(data, nels, wide+8, wide, 0, 0, 1);

        p = data + wide;
        wide += 8;
        for (i = 0; i < nels; i++, p += wide)
            array[i] = *(cptr_t *)(p);

        free(data);
    }
    return;
#endif

    for (i = 0; i < nels; i++)
        array[i] += lcp;

    radix_subsort(array, nels, 0, swaps);
    
    for (i = 0; i < nels; i++)
        array[i] -= lcp;

    free(swaps);
}
#endif

static bool sort_block(bam1_t **data, size_t n, int *lcp)
{
    if (n == 0) return true;

    // save the last pointer because the buffering logic is different
    bam1_t *lastb = *(data + n - 1);
    cptr_t lastd = (cptr_t) lastb->data;

    // change array pointer to query name
    cptr_t *array = (cptr_t *)data;
    cptr_t *x, *aend = array + n - 1;
    for (x = array; x < aend; x++)
        *x += sizeof(bam1_t);
    *aend = lastd;

#if defined USE_RADIXSORT
    hybrid_sort(array, n, lcp);
#elif defined USE_QUICKSORT
    qsort(array, n, sizeof(cptr_t), cptr_cmp);
#elif defined USE_MERGESORT
    ks_mergesort(sort, n, array, 0);
#endif

    aend = array + n;
    for (x = array; x < aend; x++) {
        if (*x == lastd)
            *x = (cptr_t) lastb;
        else
            *x -= sizeof(bam1_t);
    }

#ifdef VALIDATE_SORT
    size_t i;
    for (i = 1; i < n; i++) {
        if (strcmp(bam_get_qname(data[i-1]), bam_get_qname(data[i])) > 0) {
            errmsg("sorting validation failed at index %zu: \"%s\" > \"%s\"",
                   i - 1, bam_get_qname(data[i - 1]), bam_get_qname(data[i]));
            return false;
        }
    }
#endif
    return true;
}

static int N_THREADS = 8; // number of total threads
static int W_THREADS = 4; // number of writing threads

static pthread_mutex_t TMUTEX;
static pthread_cond_t  TCOND;
static int  *Tstack = NULL;
static int   Tavail = 0;

static pthread_mutex_t JMUTEX;
static pthread_cond_t  JCOND;

static pthread_t *Tarray = NULL;
static bool *Tjoins = NULL;

static pthread_mutex_t *DMUTEX;
static pthread_cond_t  *DCOND;

static inline int my_cond_timedwait(pthread_cond_t *cond, pthread_mutex_t *mutex, int seconds, const char *msg)
{
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    ts.tv_sec += seconds;
    int ret = pthread_cond_timedwait(cond, mutex, &ts);
    if (ret == ETIMEDOUT)
        errmsg("timed out waiting for condition: %s", msg);
    return ret;
}

#define N_CHNK 32

typedef struct {
    FILE *fp;
    ZSTD_DCtx *dctx;
    ZSTD_inBuffer inBuf; // will be static
} zst_readerArgs_t;

typedef struct {
    samFile *fp;
    sam_hdr_lite_t *header;
    int (*reader)(samFile *fp, sam_hdr_lite_t *header, bam1_t *b);
    bam1_t *b;
    size_t bbytes; // bytes left in b
    bool at_eof;
} sam_readerArgs_t;

typedef struct {
    size_t (*reader)(void *args, uint8_t *data, size_t size, bool *eof);
    void *readerArgs; // arguments for reader function
    void (*freeArgs)(void *args); // free function for reader arguments
    uint8_t *outBuf; // output buffer for decompressed data
    size_t outBufSize; // size of output buffer
    uint8_t *data[N_CHNK]; // data boundary
    bam1_t **recs[N_CHNK];
    size_t nrec[N_CHNK];
    size_t mrec; // max records per chunk
    size_t pdata; // data position
    size_t bytes; // data in bytes to process
    bool done[N_CHNK]; // track status
    bool wait; // pending data
    int mindex; // current chunk for merging
    int dindex; // current chunk for buffering
    int block; // block number
    int lcp; // shared prefix across the block
    int *lcps[N_CHNK]; // lcp to the previous record
    int *refmap; // reference ID mapping
    int tid; // thread ID
} bam1_stream_t;

typedef struct {
    void *data;
    int n_data;
} data_args_t;

static void zst_readerArgs_free(void *args) {
    zst_readerArgs_t *readerArgs = (zst_readerArgs_t *) args;
    if (readerArgs) {
        if (readerArgs->dctx) {
            ZSTD_freeDCtx(readerArgs->dctx);
            free((void *)readerArgs->inBuf.src);
        }
        free(readerArgs);
    }
}

static zst_readerArgs_t *zst_readerArgs_init(FILE *fp, size_t inBufSize) 
{
    zst_readerArgs_t *readerArgs = NULL;

    readerArgs = (zst_readerArgs_t *) calloc(1, sizeof(zst_readerArgs_t));
    if (!readerArgs) goto err;

    readerArgs->dctx = ZSTD_createDCtx();
    if (!readerArgs->dctx) goto err;
    
    if (!inBufSize) inBufSize = ZSTD_INBUF_SIZE;
    readerArgs->inBuf.src = (char *) malloc(inBufSize);
    readerArgs->inBuf.pos = 0;
    readerArgs->inBuf.size = 0;
    if (!readerArgs->inBuf.src) goto err;

    readerArgs->fp = fp;

    return readerArgs;

err:
    zst_readerArgs_free(readerArgs);
    return NULL;
}

static void sam_readerArgs_free(void *args) {
    sam_readerArgs_t *readerArgs = (sam_readerArgs_t *) args;
    if (readerArgs) {
        if (readerArgs->b)
            bam_destroy1(readerArgs->b);
        free(readerArgs);
    }
}

static sam_readerArgs_t *sam_readerArgs_init(samFile *fp, sam_hdr_lite_t *header, 
    int (*reader)(samFile *fp, sam_hdr_lite_t *header, bam1_t *b)) 
{
    sam_readerArgs_t *readerArgs = NULL;

    readerArgs = (sam_readerArgs_t *) calloc(1, sizeof(sam_readerArgs_t));
    if (!readerArgs) goto err;

    readerArgs->fp = fp;
    readerArgs->header = header;
    readerArgs->reader = reader;
    readerArgs->b = bam_init1();
    if (!readerArgs->b) goto err;
    readerArgs->bbytes = 0;
    readerArgs->at_eof = false;

    return readerArgs;
err:
    sam_readerArgs_free(readerArgs);
    return NULL;
}

static size_t zst_read_data(void *args, uint8_t *beg, size_t size, bool *eof) 
{
    zst_readerArgs_t *readerArgs = (zst_readerArgs_t *) args;
    ZSTD_inBuffer *inBuf = &readerArgs->inBuf;
    ZSTD_outBuffer outBuf = { beg, size, 0 };
    FILE *fp = readerArgs->fp;
    ZSTD_DCtx *dctx = readerArgs->dctx;
    size_t ret, bytes;
    
    while (outBuf.pos < outBuf.size) {
        if (inBuf->pos == inBuf->size) {
            // need to read more compressed data
            bytes = fread((void *)inBuf->src, 1, ZSTD_INBUF_SIZE, fp);

            if (bytes == 0) {
                if (feof(fp))
                    break;
                else if (ferror(fp)) {
                    errmsg("error reading zstd-compressed data");
                    exit(1);
                }
            }

            inBuf->pos = 0;
            inBuf->size = bytes;
        }

        if (inBuf->size > 0) {
            ret = ZSTD_decompressStream(dctx, &outBuf, inBuf);
            if (ZSTD_isError(ret)) {
                errmsg("zstd decompression error: %s", ZSTD_getErrorName(ret));
                exit(1);
            }
        }
    }

    if (eof) *eof = (inBuf->pos == inBuf->size && feof(fp));

    return outBuf.pos;
}

static inline size_t sam_fill_data(bam1_t *b, uint8_t *data, size_t size, size_t bbytes)
{
    size_t offset = sizeof(bam1_t) + b->l_data - bbytes;
    size_t bytes = bbytes < size? bbytes : size;
    size_t fill = bytes; // total bytes to fill
    if (offset < sizeof(bam1_t)) {
        // fill bam1_t structure
        size_t n = sizeof(bam1_t) - offset;
        if (n > bytes) n = bytes;
        memcpy(data, ((uint8_t *) b) + offset, n);
        data  += n;
        bytes -= n;
        offset = 0;
    } else offset -= sizeof(bam1_t);
    if (bytes)
        // fill data
        memcpy(data, b->data + offset, bytes);
    return fill;
}

static size_t sam_read_data(void *args, uint8_t *data, size_t size, bool *eof) 
{
    sam_readerArgs_t *readerArgs = (sam_readerArgs_t *) args;
    bam1_t *b = readerArgs->b;
    samFile *fp = readerArgs->fp;
    sam_hdr_lite_t *header = readerArgs->header;
    int (*reader)(samFile *fp, sam_hdr_lite_t *header, bam1_t *b) = readerArgs->reader;
    size_t fill, bytes, bbytes = readerArgs->bbytes, tbytes = size;
    int ret;

    if (bbytes) {
        fill = sam_fill_data(b, data, size, bbytes);
        size   -= fill;
        data   += fill;
        bbytes -= fill;
    }
    while (size > 0) {
        ret = reader(fp, header, b);
        if (ret < -1) {
            errmsg("error reading SAM/BAM record");
            exit(1);
        } else if (ret == -1) {
            // end of file
            readerArgs->at_eof = true;
            break;
        }
        // got a record
        // bbytes is zero at this point
        bbytes = sizeof(bam1_t) + b->l_data;
        fill = sam_fill_data(b, data, size, bbytes);
        size   -= fill;
        data   += fill;
        bbytes -= fill;
    }

    readerArgs->bbytes = bbytes;
    if (eof) *eof = (bbytes == 0 && readerArgs->at_eof);

    return tbytes - size;
}

static void *bam_stream_data(void *args) 
{
    bam1_stream_t *stream = (bam1_stream_t *) args;
    void *readerArgs = stream->readerArgs;
    size_t (*reader)(void *args, uint8_t *data, size_t size, bool *eof) = stream->reader;
    uint8_t *dst;
    size_t bytes, size, head, tail, tsize, nrec, mrec, sbyte, tbyte;
    bam1_t **recs, *b;
    int *refmap;
    uint8_t *pos, *end, *p, *q;
    int i, lcp, *lcps, block, index, mindex, dindex;
    bool eof = false;
    
    block = stream->block;
    mrec = stream->mrec;
    lcp = stream->lcp;
    dindex = stream->dindex;
    recs = stream->recs[dindex];
    lcps = stream->lcps[dindex];
    refmap = stream->refmap;
    dst = stream->outBuf;
    tsize = stream->outBufSize;
    tbyte = tsize / N_CHNK;
    bytes = stream->bytes;
    pos = dst + stream->pdata;
    end = dst + tsize;

    // find beginning of chunk
    if (pos >= end) pos = dst; // wrap around
    // find size of free space
    mindex = (dindex + 1) % N_CHNK;
    pthread_mutex_lock(DMUTEX + block);
    while (!stream->done[mindex] && mindex != dindex)
        mindex = (mindex + 1) % N_CHNK;
    pthread_mutex_unlock(DMUTEX + block);
    if (mindex == dindex) {
        // all chunks are empty
        tail = end - pos;
        head = pos - dst;
    } else {
        if (stream->data[mindex] >= pos) {
            tail = stream->data[mindex] - pos;
            head = 0;
        } else {
            // wrap around
            tail = end - pos;
            head = stream->data[mindex] - dst;
        }
    }
    size = tail + head;
    if (!size) { // no space available
        pthread_mutex_lock(DMUTEX + block);
        if (stream->done[mindex]) // recheck if the status changed
            stream->wait = true; // pause buffering until further notice
        pthread_mutex_unlock(DMUTEX + block);
        goto release_worker_thread;
    }

    // DO NOT call this routine if nrec > 0
    // otherwise the data could be overwritten
    // buffer more data
    // more data if tail space available
    if (tail > bytes)
        bytes += reader(readerArgs, pos + bytes, tail - bytes, &eof);
    // try memmove and fill if necessary
    while (sizeof(bam1_t) > bytes ||
        sizeof(bam1_t) + (b = (bam1_t *) pos)->l_data > bytes) {
        if (eof) {
            if (bytes) {
                errmsg("incomplete BAM record - probably corrupted file");
                exit(1);
            } else break; // end of file
        } 
        // move to beginning if possible
        if (pos != dst) {
            if (size == tsize) head = size;
            if (head > bytes) {
                memmove(dst, pos, bytes);
                pos = dst;
            }
        }
        // still cannot make a full record
        if (bytes >= tsize) {
            // very unlikely but possible for extremely large records
            // if reaching here, it is guaranteed that pos == dst
            // also all the chunks are empty
            tsize <<= 1;
            uint8_t *newBuf = (uint8_t *) malloc(tsize);
            if (!newBuf) {
                errmsg("failed to allocate memory");
                exit(1);
            }
            memcpy(newBuf, pos, bytes);
            free(dst);
            head = size = tsize;
            tbyte = tsize / N_CHNK;
            pos = dst = newBuf;
            stream->outBuf = newBuf;
            stream->outBufSize = tsize;
        }
        // more data if head space available
        if (head > bytes)
            bytes += reader(readerArgs, pos + bytes, head - bytes, &eof);
        else break; // no more space
    }

    // check end of file
    if (eof && !bytes) {
        stream->nrec[dindex] = 0;
        pthread_mutex_lock(DMUTEX + block);
        stream->done[dindex] = true;
        stream->dindex = -1; // mark end of stream
        pthread_cond_signal(DCOND + block);
        pthread_mutex_unlock(DMUTEX + block);
        goto release_worker_thread;
    }
    // make sure we have enough data for the next record
    if (sizeof(bam1_t) > bytes || 
        sizeof(bam1_t) + (b = (bam1_t *) pos)->l_data > bytes) {
        if (size == tsize) {
            errmsg("incomplete BAM record - probably corrupted file");
            exit(1);
        } else {
            // likely due to full buffer - simply return
            stream->pdata = pos - dst;
            stream->bytes = bytes;
            pthread_mutex_lock(DMUTEX + block);
            if (stream->done[mindex]) // recheck if the status changed
                stream->wait = true; // pause buffering until further notice
            pthread_mutex_unlock(DMUTEX + block);
            goto release_worker_thread;
        }
    }
    // now we have at least one full record
    nrec = 0;
    sbyte = bytes;
    stream->data[dindex] = pos;
    while (nrec < mrec) {
        pos += sizeof(bam1_t);
        b->data = pos;
        pos += b->l_data;
        bytes -= sizeof(bam1_t) + b->l_data;
        recs[nrec++] = b;
        if (sbyte - bytes >= tbyte || // enough data for this chunk
            sizeof(bam1_t) > bytes || 
            sizeof(bam1_t) + (b = (bam1_t *) pos)->l_data > bytes)
            break;
    }
    stream->pdata = pos - dst;
    stream->bytes = bytes;
    stream->nrec[dindex] = nrec;

    // mapping reference
    if (refmap)
        for (i = 0; i < nrec; i++)
            recs[i]->core.tid = refmap[recs[i]->core.tid];

    // calculate lcp array
    lcps[0] = 0;
    for (i = 1; i < nrec; i++) {
        p = recs[i-1]->data + lcp;
        q = recs[i]->data + lcp;
        while (*p && *p == *q) p++, q++;
        if (*p == *q) // reach end of identical strings
            lcps[i] = -1;
        else
            lcps[i] = q - recs[i]->data;
    }

    // mark this chunk as done
    pthread_mutex_lock(DMUTEX + block);
    stream->done[dindex] = true;
    stream->dindex = (dindex + 1) % N_CHNK; // next chunk
    pthread_cond_signal(DCOND + block);
    pthread_mutex_unlock(DMUTEX + block);

release_worker_thread:
    pthread_mutex_lock(&TMUTEX);
    Tstack[Tavail++] = stream->tid;
    pthread_cond_signal(&TCOND);
    pthread_mutex_unlock(&TMUTEX);

    return NULL;
}

static int Jfree = 0;

static int B_THREADS = 1; // threads for data buffering
static int S_THREADS = 1; // threads for sorting
static int Z_THREADS = 0; // this could be changed dynamically

timer1_t *BUFFER_WAIT;
double BUFFER_TIME = 0.0;

static void *data_buffer_thread(void *args) {
    data_args_t *data = (data_args_t *) args;
    int n_streams = data->n_data;
    bam1_stream_t *streams = (bam1_stream_t *) data->data, *stream;
    int i, tid, block, chunk, Javail, Jawait;

    for (i = 0; i < B_THREADS; i++) {
        Tstack[i] = i;
        Tjoins[i] = false;
    }
    Tavail = B_THREADS;
    
    Javail = n_streams;
    while (Javail) {
        pthread_mutex_lock(&JMUTEX);
        Jfree = 0;
        pthread_mutex_unlock(&JMUTEX);
        Javail = n_streams;
        Jawait = 0;
        for (i = 0; i < n_streams; i++) {
            stream = streams + i;
            block = stream->block;
            chunk = stream->dindex;
            // this stream has ended
            if (chunk < 0) {
                Javail--;
                continue;
            }
            // wait for merging to be done
            pthread_mutex_lock(DMUTEX + block);
            if (stream->done[chunk] || stream->wait) {
                Jawait++;
                pthread_mutex_unlock(DMUTEX + block);
                continue;
            }
            pthread_mutex_unlock(DMUTEX + block);
            // buffer data with next available thread
            pthread_mutex_lock(&TMUTEX);
            if (Tavail <= 0)
                pthread_cond_wait(&TCOND,&TMUTEX);
            tid = Tstack[--Tavail];
            pthread_mutex_unlock(&TMUTEX);
            stream->tid = tid;
            if (Tjoins[tid]) pthread_join(Tarray[tid], NULL);
            pthread_create(Tarray + tid, NULL, bam_stream_data, stream);
            Tjoins[tid] = true;
        }
        // wait for all threads to finish
        pthread_mutex_lock(&TMUTEX);
        while (Tavail < B_THREADS)
            pthread_cond_wait(&TCOND,&TMUTEX);
        pthread_mutex_unlock(&TMUTEX);
        // all jobs are waiting for merging
        pthread_mutex_lock(&JMUTEX);
        if (Jawait == Javail && Javail && !Jfree)
            pthread_cond_wait(&JCOND,&JMUTEX);
        pthread_mutex_unlock(&JMUTEX);
    }

    for (i = 0; i < B_THREADS; i++) {
        if (Tjoins[i]) {
            pthread_join(Tarray[i], NULL);
            Tjoins[i] = false;
        }
    }

    return NULL;
}

typedef struct {
    int n;
    int i;
    int *lcps;
    bam1_t **recs;
    bam1_stream_t *data;
    char *dups; // for strdup
} heap_t;

static bool yield(int t, void *data, char **s, int *p)
{
    heap_t *ls = (heap_t *) data + t;
    if (ls->i < ls->n) {
        *s = (char *) ls->recs[ls->i]->data;
        *p = ls->lcps[ls->i++];
        return true;
    }
    if (ls->data) {
        bam1_stream_t *stream = ls->data;
        int chunk = stream->mindex;

        pthread_mutex_lock(DMUTEX+t);
        stream->nrec[chunk] = 0;
        stream->done[chunk] = false;
        stream->wait = false;
        pthread_mutex_unlock(DMUTEX+t);

        pthread_mutex_lock(&JMUTEX);
        Jfree++;
        pthread_cond_signal(&JCOND);
        pthread_mutex_unlock(&JMUTEX);

        chunk = (chunk + 1) % N_CHNK;
        pthread_mutex_lock(DMUTEX+t);
        if (!stream->done[chunk])
            pthread_cond_wait(DCOND+t,DMUTEX+t);
        stream->mindex = chunk;
        pthread_mutex_unlock(DMUTEX+t);
        ls->i = 0;
        ls->n = stream->nrec[chunk];
        ls->recs = stream->recs[chunk];
        ls->lcps = stream->lcps[chunk];
        
        if (ls->i < ls->n) {
            *s = (char *) ls->recs[ls->i]->data;
            *p = ls->lcps[ls->i++];
            return true;
        }
    }
    return false;
}

typedef struct {
    int T;
    heap_t *heaps;
    void *writer;
} merge_args_t;

static void *bam_merge_thread(void *args)
{
    merge_args_t *margs = (merge_args_t *) args;
    int T = margs->T;
    heap_t *heaps = margs->heaps;
    bam_writer_t *writer = (bam_writer_t *) margs->writer;
    int i, j, n, m, t, *list;
    size_t n_recs = 0, n_prev = 0;
    bam1_t **recs;
    int res, *lcps;
    char *dups;
    Merge *merge = mergeCreateString(T, heaps, yield);
    push_time(&TIME_STACK);
    // do merging
    while ((m = mergeNext(merge, &list))) {
        for (i = 0; i < m; i++) {
            t = list[i];
            j = heaps[t].i - 1;
            n = heaps[t].n;
            recs = heaps[t].recs;
            lcps = heaps[t].lcps;
            res = bw_bam_write1(writer, recs[j]);
            while (++j < n && lcps[j] < 0)
                res |= bw_bam_write1(writer, recs[j]);
            if (res < 0) {
                errmsg("failed writing BAM record");
                exit(1);
            }
            n_recs += j - heaps[t].i + 1;
            heaps[t].i = j;
            if (j == n) {
                // end of heap
                // need strdup as the buffer will be refilled
                // and the data pointer could become invalid
                dups = strdup((char *)recs[n-1]->data);
                if (!mergeUpdateString(merge, t, dups)) {
                    errmsg("failed updating merge string");
                    exit(1);
                }
                free(heaps[t].dups);
                heaps[t].dups = dups;
            }
        }
        if (n_recs >= 10000000 + n_prev) {
            infomsg("processed %12zu records in %.3f seconds", 
                n_recs, elapsed_time(&TIME_STACK, false));
            n_prev = n_recs;
        }
    }
    infomsg("processed %12zu records in %.3f seconds", 
        n_recs, pop_time(&TIME_STACK));

    mergeDestroy(merge);
    return NULL;
}

static void *zst_merge_thread(void *args)
{
    merge_args_t *margs = (merge_args_t *) args;
    int T = margs->T;
    heap_t *heaps = margs->heaps;
    zstd_writer_t *fp = (zstd_writer_t *)margs->writer;
    int i, j, n, m, t, *list;
    size_t n_recs = 0, n_prev = 0;
    bam1_t **recs;
    int res, *lcps;
    char *dups;
    Merge *merge = mergeCreateString(T, heaps, yield);
    push_time(&TIME_STACK);
    // do merging
    while ((m = mergeNext(merge, &list))) {
        for (i = 0; i < m; i++) {
            t = list[i];
            j = heaps[t].i - 1;
            n = heaps[t].n;
            recs = heaps[t].recs;
            lcps = heaps[t].lcps;
            res = zstd_bam_write1(fp, recs[j]);
            while (++j < n && lcps[j] < 0)
                res |= zstd_bam_write1(fp, recs[j]);
            if (res < 0) {
                errmsg("failed writing BAM record");
                exit(1);
            }
            n_recs += j - heaps[t].i + 1;
            heaps[t].i = j;
            if (j == n) {
                // end of heap
                // need strdup as the buffer will be refilled
                // and the data pointer could become invalid
                dups = strdup((char *)recs[n-1]->data);
                if (!mergeUpdateString(merge, t, dups)) {
                    errmsg("failed updating merge string");
                    exit(1);
                }
                free(heaps[t].dups);
                heaps[t].dups = dups;
            }
        }
        if (n_recs >= 10000000 + n_prev) {
            infomsg("processed %12zu records in %.3f seconds", 
                n_recs, elapsed_time(&TIME_STACK, false));
            n_prev = n_recs;
        }
    }
    if (zstd_writer_flush(fp) < 0) {
        errmsg("failed to flush ZSTD writer");
        exit(1);
    }
    infomsg("processed %12zu records in %.3f seconds", 
        n_recs, pop_time(&TIME_STACK));
    
    mergeDestroy(merge);
    return NULL;
}

static bool merge_blocks(bam1_stream_t *streams, int n_streams, void *writer, void *(*merge)(void *))
{
    int i;
    bool ret = false;
    int T = n_streams;

    // initialise locks
    DMUTEX = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t) * T);
    DCOND = (pthread_cond_t *) malloc(sizeof(pthread_cond_t) * T);
    // data buffering threads
    B_THREADS = N_THREADS - W_THREADS;
    if (B_THREADS < 1) B_THREADS = 1;
    Tarray = (pthread_t *) malloc(B_THREADS * sizeof(pthread_t));
    Tjoins = (bool *) malloc(B_THREADS * sizeof(bool));
    Tstack  = (int *) malloc(B_THREADS * sizeof(int));
    if (DMUTEX == NULL || DCOND == NULL || Tarray == NULL || 
        Tjoins == NULL || Tstack == NULL) {
        errmsg("threads memory allocation failed");
        goto err;
    }
    for (i = 0; i < T; i++) {
        pthread_mutex_init(DMUTEX+i,NULL);
        pthread_cond_init(DCOND+i,NULL);
    }
    pthread_mutex_init(&JMUTEX,NULL);
    pthread_cond_init(&JCOND,NULL);
    pthread_mutex_init(&TMUTEX,NULL);
    pthread_cond_init(&TCOND,NULL);
    // spawn threads for data buffering
    pthread_t data_thread;
    data_args_t data_args = { .n_data = T, .data = streams };
    pthread_create(&data_thread, NULL, data_buffer_thread, &data_args);
    // wait for first chunk of each stream to be buffered
    for (i = 0; i < T; i++) {
        pthread_mutex_lock(DMUTEX+i);
        if (!streams[i].done[0])
            pthread_cond_wait(DCOND+i,DMUTEX+i);
        pthread_mutex_unlock(DMUTEX+i);
    }
    
    // build data structures for heaps
    heap_t *heaps = (heap_t *) calloc(T, sizeof(heap_t));
    if (!heaps) {
        errmsg("failed to allocate memory for heaps");
        goto err;
    }
    for (i = 0; i < T; i++) {
        heaps[i].i = 0;
        heaps[i].n = streams[i].nrec[0];
        heaps[i].recs = streams[i].recs[0];
        heaps[i].lcps = streams[i].lcps[0];
        heaps[i].data = streams + i;
        heaps[i].dups = NULL;
    }
    // spawn merging thread
    pthread_t merge_thread;
    merge_args_t merge_args = { .T = T, .heaps = heaps, .writer = writer };
    pthread_create(&merge_thread, NULL, merge, &merge_args);
    // wait for merging to finish
    pthread_join(merge_thread, NULL);
    // wait for data buffering to finish
    pthread_join(data_thread, NULL);
    ret = true;

err:
    for (i = 0; i < T; i++) {
        pthread_mutex_destroy(DMUTEX+i);
        pthread_cond_destroy(DCOND+i);
    }
    free(DMUTEX);
    free(DCOND);
    pthread_mutex_destroy(&JMUTEX);
    pthread_cond_destroy(&JCOND);
    pthread_mutex_destroy(&TMUTEX);
    pthread_cond_destroy(&TCOND);

    free(Tarray), Tarray = NULL;
    free(Tjoins), Tjoins = NULL;
    free(Tstack), Tstack = NULL;

    for (i = 0; i < T; i++)
        free(heaps[i].dups);
    free(heaps);

    return ret;
}

#define BAM_BLOCK_SIZE 2*1024*1024
const size_t SORT_DEFAULT_MEGS_PER_THREAD = 768ULL;
const size_t SORT_MIN_MEGS_PER_THREAD = 1ULL;

const size_t MAX_MEM_PER_THREAD = SORT_DEFAULT_MEGS_PER_THREAD<<20;

static void complain_about_memory_setting(size_t max_mem) {
    char  *suffix = "";
    const size_t nine_k = 9<<10;
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "K"; }
    if (max_mem > nine_k) { max_mem >>= 10; suffix = "M"; }

    fprintf(stderr,
"[bamsort] -m setting (%zu%s bytes) is less than the minimum required (%zuM).\n\n"
"Trying to run with -m too small can lead to the creation of a very large number\n"
"of temporary files.  This may make sort fail due to it exceeding limits on the\n"
"number of files it can have open at the same time.\n\n"
"Please check your -m parameter.  It should be an integer followed by one of the\n"
"letters K (for kilobytes), M (megabytes) or G (gigabytes).  You should ensure it\n"
"is at least the minimum above, and much higher if you are sorting a large file.\n",
            max_mem, suffix, SORT_MIN_MEGS_PER_THREAD);
}

enum FILE_TYPE { FT_UNKNOWN, FT_BAM, FT_SAM, FT_CRAM, FT_ZSTD };

static inline int64_t bam_ftell(samFile *fp) 
{
    BGZF *bgzf = fp->fp.bgzf;
    return (bgzf->block_address << 16) | (bgzf->block_offset & 0xFFFF);
}

static inline int64_t sam_ftell(samFile *fp)
{
    struct hFILE *hfile = fp->fp.hfile;
    return (hfile->offset + (hfile->begin - hfile->buffer));
}

typedef struct {
    // consider making reader generic to support ZSTD input
    enum FILE_TYPE ftype;
    samFile *fp;
    int64_t (*ftell)(samFile *fp);
    sam_hdr_lite_t *header;
    int (*reader)(samFile *fp, sam_hdr_lite_t *header, bam1_t *b);
    int64_t boff;
    int64_t eoff;
    uint8_t *data;
    size_t mmem;
    bam1_t **recs;
    size_t mrec;
    bam1_t *b;
    size_t nrec;
    int lcp;
    int *lcps;
    int *refs;
    int jid;
    int tid;
} bio_data_t;

typedef struct {
    enum FILE_TYPE ftype;
    int nfile;
    void **fns;
    int *lcps;
    size_t nrec;
    int *refs;
    sam_hdr_lite_t *header;
} srt_data_t;

static sam_hdr_lite_t *bio_create_partitions(char *infile, int n_part, bio_data_t *bio, int *_part);

static int *SJlist = NULL;
static int  SJn = 0;
static int *WJlist = NULL;
static int  WJn = 0;

static pthread_mutex_t WMUTEX;
static pthread_cond_t  WCOND;

static int   SNdone = 0;
static bool *SJdone = NULL;

static void *bio_sort_thread(void *args)
{
    bio_data_t *bio = (bio_data_t *) args;
    samFile *fp = bio->fp;
    int64_t eoff = bio->eoff;
    int64_t (*ftell)(samFile *fp) = bio->ftell;
    uint8_t *data = bio->data;
    bam1_t **recs = bio->recs;
    size_t mrec = bio->mrec;
    bam1_t *b = bio->b;
    sam_hdr_lite_t *header = bio->header;
    int (*reader)(samFile *fp, sam_hdr_lite_t *header, bam1_t *b) = bio->reader;
    int *refs = bio->refs;
    int *lcps = bio->lcps;
    int tid = bio->tid;
    int jid = bio->jid;
    uint8_t *last = data + bio->mmem;
    size_t nrec = 0;
    bam1_t *rec;
    int i, t, lcp = 0, res = 0;
    bool done = false;

    while ((res = reader(fp, header, b)) >= 0) {
        if (ftell(fp) > eoff) {
            done = true;
            break; // file chunk done
        }
        if (nrec >= mrec)
            break; // record limit reached
        // update reference map
        if (refs)
            __atomic_store_n(refs+b->core.tid, 1, __ATOMIC_RELAXED);
        // copy the record
        if (data + sizeof(*b) + b->l_data < last) {
            rec = recs[nrec++] = (bam1_t *) data;
            *rec = *b;
            rec->data = data + sizeof(bam1_t);
            memcpy(rec->data, b->data, b->l_data);
            data += (sizeof(*b) + b->l_data + 8 - 1) & ~((size_t)(8 - 1));
        } else break; // memory limit reached
    }
    if (res == -1)
        done = true; // end of file
    else if (res < -1) {
        errmsg("error reading BAM record");
        exit(1);
    }
    // include the last record if not done
    // there is extra space preallocated for the last record
    if (!done) {
        recs[nrec++] = b;
        // update reference map
        if (refs)
            __atomic_store_n(refs+b->core.tid, 1, __ATOMIC_RELAXED);
    }

    if (nrec) {
        // sort the records
        if (!sort_block(recs, nrec, &lcp)) {
            errmsg("failed to sort BAM records");
            exit(1);
        }
        // calculate lcp array
        uint8_t *p, *q;
        lcps[0] = 0;
        for (i = 1; i < nrec; i++) {
            p = recs[i-1]->data + lcp;
            q = recs[i]->data + lcp;
            while (*p && *p == *q) p++, q++;
            if (*p == *q) // reach end of identical strings
                lcps[i] = -1;
            else
                lcps[i] = q - recs[i]->data;
        }
    }

    bio->lcp = lcp;
    bio->nrec = nrec;

    // add to writing job list
    pthread_mutex_lock(&WMUTEX);
    WJlist[WJn++] = jid;
    // set job status if done
    if (done) {
        pthread_mutex_lock(&JMUTEX);
        SJdone[jid] = true;
        pthread_mutex_unlock(&JMUTEX);
    }
    pthread_cond_signal(&WCOND);
    pthread_mutex_unlock(&WMUTEX);
    // reclaim worker thread - need broadcast here
    pthread_mutex_lock(&TMUTEX);
    Tstack[Tavail++] = tid;
    pthread_mutex_unlock(&TMUTEX);
    pthread_cond_broadcast(&TCOND);
    
    return NULL;
}

static void *sort_thread(void *args)
{
    data_args_t *data = (data_args_t *) args;
    int n_bio = data->n_data;
    bio_data_t *bio = (bio_data_t *) data->data;
    
    int i, tid, jid;
    
    for (i = 0; i < S_THREADS; i++) {
        Tstack[i] = i;
        Tjoins[i] = false;
    }
    Tavail = S_THREADS;
    for (i = 0; i < n_bio; i++) {
        SJlist[i] = i;
        SJdone[i] = false;
    }
    SJn = n_bio;
    SNdone = 0;

    while (1) {
        pthread_mutex_lock(&JMUTEX);
        if (SNdone >= n_bio) {
            pthread_mutex_unlock(&JMUTEX);
            break;
        }
        if (SJn <= 0)
            pthread_cond_wait(&JCOND,&JMUTEX);
        if (SJn <= 0) {
            pthread_mutex_unlock(&JMUTEX);
            continue;
        }
        jid = SJlist[--SJn];
        pthread_mutex_unlock(&JMUTEX);
        pthread_mutex_lock(&TMUTEX);
        while (Tavail <= Z_THREADS)
            pthread_cond_wait(&TCOND,&TMUTEX);
        tid = Tstack[--Tavail];
        pthread_mutex_unlock(&TMUTEX);
        if (Tjoins[tid]) pthread_join(Tarray[tid], NULL);
        bio[jid].tid = tid;
        pthread_create(Tarray + tid, NULL, bio_sort_thread, bio + jid);
        Tjoins[tid] = true;
    }

    for (i = 0; i < S_THREADS; i++) {
        if (Tjoins[i]) {
            pthread_join(Tarray[i], NULL);
            Tjoins[i] = false;
        }
    }

    return NULL;
}

static FILE *new_tempfile(char *fn_template, bool delete)
{
    char *template = strdup(fn_template);
    int fd = mkstemp(template);
    if (fd == -1) {
        errmsg("mkstemp failed");
        return NULL;
    }
    FILE *fn = fdopen(fd, "w+b");
    if (!fn) {
        errmsg("fdopen failed");
        close(fd);
        return NULL;
    }
    if (delete) unlink(template);
    free(template);
    return fn;
}

static int common_lcp(bio_data_t *bio, int *list, int n)
{
    int i, l, lcp = bio[list[0]].lcp;
    char *f, *p, *q;

    f = NULL;
    for (i = 0; i < n; i++) {
        if (bio[list[i]].nrec) {
            f = (char *) bio[list[i]].recs[0]->data;
            break;
        }
    }
    for (i = i+1; i < n; i++) {
        if (!bio[list[i]].nrec) continue;
        p = f;
        q = (char *) bio[list[i]].recs[0]->data;
        l = bio[list[i]].lcp;
        if (l > lcp) l = lcp;
        while (l-- && *p == *q) p++, q++;
        lcp = p - f;
    }

    return lcp;
}

static bool yield1(int t, void *data, char **s, int *p)
{
    heap_t *ls = (heap_t *) data + t;
    if (ls->i < ls->n) {
        *s = (char *) ls->recs[ls->i]->data;
        *p = ls->lcps[ls->i++];
        return true;
    } else return false;
}

static bool bam_merge(heap_t *heaps, int T, void *writer)
{
    int i, j, n, m, t, res, *list, *lcps;
    bam1_t **recs;
    bool ret = false;
    bam_writer_t *bw = (bam_writer_t *) writer;
    Merge *merge = mergeCreateString(T, heaps, yield1);
    // do merging
    while ((m = mergeNext(merge, &list))) {
        for (i = 0; i < m; i++) {
            t = list[i];
            j = heaps[t].i - 1;
            n = heaps[t].n;
            recs = heaps[t].recs;
            lcps = heaps[t].lcps;
            res = bw_bam_write1(bw, recs[j]);
            while (++j < n && lcps[j] < 0)
                res |= bw_bam_write1(bw, recs[j]);
            if (res < 0) {
                errmsg("failed writing BAM record");
                goto err;
            }
            heaps[t].i = j;
        }
    }
    ret = true;

err:
    mergeDestroy(merge);
    return ret;
}


static bool zst_merge(heap_t *heaps, int T, void *writer)
{
    int i, j, n, m, t, res, *list, *lcps;
    bam1_t **recs;
    bool ret = false;
    zstd_writer_t *fp = (zstd_writer_t *)writer;
    Merge *merge = mergeCreateString(T, heaps, yield1);
    // do merging
    while ((m = mergeNext(merge, &list))) {
        for (i = 0; i < m; i++) {
            t = list[i];
            j = heaps[t].i - 1;
            n = heaps[t].n;
            recs = heaps[t].recs;
            lcps = heaps[t].lcps;
            res = zstd_bam_write1(fp, recs[j]);
            while (++j < n && lcps[j] < 0)
                res |= zstd_bam_write1(fp, recs[j]);
            if (res < 0) {
                errmsg("failed writing BAM record");
                goto err;
            }
            heaps[t].i = j;
        }
    }
    if (zstd_writer_flush(writer) < 0) {
        errmsg("failed to flush ZSTD writer");
        goto err;
    }
    ret = true;

err:
    mergeDestroy(merge);
    return ret;
}

static bool merge_blocks_in_mem(bio_data_t *bio, int n_bio, int *blist, void *writer, bool (*merge)(heap_t *, int, void *))
{
    bool ret = false;
    int i, T = n_bio;

    heap_t *heaps = calloc(T, sizeof(heap_t));
    if (!heaps) {
        errmsg("failed to allocate memory for heaps");
        goto err;
    }
    if (blist) {
        for (i = 0; i < T; i++) {
            heaps[i].i = 0;
            heaps[i].n = bio[blist[i]].nrec;
            heaps[i].recs = bio[blist[i]].recs;
            heaps[i].lcps = bio[blist[i]].lcps;
        }
    } else {
        for (i = 0; i < T; i++) {
            heaps[i].i = 0;
            heaps[i].n = bio[i].nrec;
            heaps[i].recs = bio[i].recs;
            heaps[i].lcps = bio[i].lcps;
        }
    }

    ret = merge(heaps, T, writer);

err:
    free(heaps);
    return ret;
}

static bool bamsort1(char *infile, char *fn_template, bool trim, bool multi_in, bio_data_t *bio, srt_data_t *srt, 
    refseq_builder_t *rdb, bool in_mem_check, void *writer, bool (*merge)(heap_t *, int, void *))
{
    sam_hdr_lite_t *header = NULL;
    size_t n_recs = 0, nrec;
    array_t(void *) fns = {0, 0, NULL};
    array_t(int) lcps = {0, 0, NULL};
    FILE *fn;
    int *refs = NULL;
    int res, i, n_part = 0;
    bool ret = false;

    push_time(&TIME_STACK);

    // open n_threads * 2 files for data streaming
    if ((header = bio_create_partitions(infile, S_THREADS*2, bio, &n_part)) == NULL) {
        errmsg("failed to partition input file: %s", infile);
        goto err;
    }
    srt->header = header;

    if (trim || multi_in) { // can only skip for single non-trimmed BAM input
        push_time(&TIME_STACK);
        // need a reference sequence map
        if (refseq_builder_add(rdb, header, &refs) < 0) {
            errmsg("failed to build reference sequence map");
            goto err;
        }
        infomsg("make reference map in %.3f seconds", pop_time(&TIME_STACK));
        srt->refs = refs;
        refs += 1; refs[-1] = -1; // refs[-1] for unmapped reads
    }

    // set reference sequence map for each partition
    for (i = 0; i < n_part; i++)
        bio[i].refs = refs;

    // spawn threads for sorting
    pthread_t srt_thread;
    data_args_t data = { .n_data = n_part, .data = bio };
    pthread_create(&srt_thread, NULL, sort_thread, &data);

    bool in_mem = true;
    WJn = 0;
    // for a single input file, check if all sorting jobs can be fit in memory
    if (in_mem_check) {
        // wait for the first batch of sorting jobs to be ready
        pthread_mutex_lock(&WMUTEX);
        while (WJn < n_part)
            pthread_cond_wait(&WCOND,&WMUTEX);
        pthread_mutex_unlock(&WMUTEX);
        for (i = 0; i < n_part; i++) {
            if (!SJdone[i]) {
                in_mem = false;
                break;
            }
        }
        if (in_mem) {
            pthread_mutex_lock(&JMUTEX);
            SNdone = n_part;
            pthread_cond_signal(&JCOND);
            pthread_mutex_unlock(&JMUTEX);
            for (i = 0; i < n_part; i++)
                n_recs += bio[i].nrec;
        }
    } else in_mem = false;
    
    if (!in_mem) {
        int wJ = n_part > S_THREADS? S_THREADS : n_part;
        push_time(&TIME_STACK);
        while (wJ) {
            // no writing job yet - wait for sorting to finish
            pthread_mutex_lock(&WMUTEX);
            if (WJn >= wJ) {
                pthread_mutex_lock(&TMUTEX);
                Z_THREADS = W_THREADS;
                pthread_mutex_unlock(&TMUTEX);
            } else {
                pthread_cond_wait(&WCOND,&WMUTEX);
                pthread_mutex_unlock(&WMUTEX);
                continue;
            }
            pthread_mutex_unlock(&WMUTEX);
            
            // writing job available - collect worker threads
            pthread_mutex_lock(&TMUTEX);
            while (Tavail < Z_THREADS)
                pthread_cond_wait(&TCOND,&TMUTEX);
            pthread_mutex_unlock(&TMUTEX);
                
            // do merge writing
            nrec = 0;
            for (i = 0; i < wJ; i++)
                nrec += bio[WJlist[i]].nrec;
            n_recs += nrec;
            if (nrec > 0) {
                fn = new_tempfile(fn_template, true);
                if (!fn) {
                    errmsg("failed to create temporary file");
                    goto err;
                }
                array_push(void *, fns, (void *)fn);
                array_push(int, lcps, common_lcp(bio, WJlist, wJ));
                // set up zstd writer
                ((zstd_writer_t *) writer)->fp = fn;
                if (!merge_blocks_in_mem(bio, wJ, WJlist, writer, merge)) {
                    errmsg("failed to merge write BAM blocks");
                    goto err;
                }
                // rewind for reading later
                if (fseek(fn, 0, SEEK_SET) != 0) {
                    errmsg("failed to rewind temporary file");
                    goto err;
                }
                infomsg("wrote %12zu sorted records to temporary files in %.3f seconds", nrec, elapsed_time(&TIME_STACK, true));
            }
            // reclaim worker threads and signal sorting threads
            pthread_mutex_lock(&TMUTEX);
            Z_THREADS = 0;
            pthread_cond_signal(&TCOND);
            pthread_mutex_unlock(&TMUTEX);
            // release sorting jobs
            pthread_mutex_lock(&JMUTEX);
            for (i = 0; i < wJ; i++)
                if (!SJdone[WJlist[i]])
                    SJlist[SJn++] = WJlist[i];
                else
                    SNdone++;
            pthread_cond_signal(&JCOND);
            pthread_mutex_unlock(&JMUTEX);
            // shift remaining writing jobs to top
            pthread_mutex_lock(&WMUTEX);
            for (i = wJ; i < WJn; i++)
                WJlist[i - wJ] = WJlist[i];
            WJn -= wJ;
            pthread_mutex_unlock(&WMUTEX);
            // adjust number of writing jobs - thread safe
            if (wJ + SNdone > n_part)
                wJ = n_part - SNdone;
        }
        pop_time(&TIME_STACK);
    }
    
    pthread_join(srt_thread, NULL);

    if (in_mem)
        infomsg("wrote %12zu sorted records in memory in %.3f seconds", n_recs, pop_time(&TIME_STACK));
    else
        infomsg("wrote %12zu sorted records to [%d] temporary files in %.3f seconds", n_recs, fns.n, pop_time(&TIME_STACK));

    // if non-trimming, we need to bring back all reference sequences
    if (refs && !trim) // multi-in non-trimming
        for (i = 0; i < header->n_targets; i++)
            refs[i] = 1;

    srt->nrec = n_recs;
    srt->fns = fns.a;
    srt->nfile = fns.n;
    srt->lcps = lcps.a;

    ret = true;

err:
    // close file handles
    for (i = 0; i < n_part; i++) {
        if (bio[i].fp) {
            sam_close((samFile *)bio[i].fp);
            bio[i].fp = NULL;
        }
    }

    return ret;
}

const size_t BAM_STREAM_MIN_MEM = ZSTD_OUTBUF_SIZE;
const size_t BAM_STREAM_MAX_MEM = 1ULL<<32; // 4 GB

static bam1_stream_t *create_bam1_streams(srt_data_t *srt_data, int nIn, size_t max_mem, int *_n_bs)
{
    int i, j, k, n, n_bs;
    uint64_t max_rec;
    enum FILE_TYPE ftype;
    
    if (_n_bs) *_n_bs = 0;

    n_bs = 0;
    for (i = 0; i < nIn; i++)
        n_bs += srt_data[i].nfile;
    if (n_bs == 0) return NULL;

    // estimate memory per file and chunk size
    max_mem /= n_bs;
    if (max_mem < BAM_STREAM_MIN_MEM)
        max_mem = BAM_STREAM_MIN_MEM;
    if (max_mem > BAM_STREAM_MAX_MEM)
        max_mem = BAM_STREAM_MAX_MEM;
    max_rec = max_mem / BAM_REC_SIZE / N_CHNK;
    infomsg("number chunks per stream: %d", N_CHNK);
    infomsg("maximum mems per chunk: %zu Mb", max_mem>>20);
    infomsg("maximum recs per chunk: %d", max_rec);

    // create decompressor for each input file
    bam1_stream_t *bams = (bam1_stream_t *) calloc(n_bs, sizeof(bam1_stream_t));
    if (!bams) {
        errmsg("failed to create BAM stream context");
        return NULL;
    }
    void **fns;
    int nfile, *lcps, *refs;
    sam_hdr_lite_t *hdr;
    zst_readerArgs_t *zstReaderArgs;
    sam_readerArgs_t *samReaderArgs;

    n = 0;
    for (i = 0; i < nIn; i++) {
        hdr = srt_data[i].header;
        ftype = srt_data[i].ftype;
        fns = srt_data[i].fns;
        nfile = srt_data[i].nfile;
        lcps = srt_data[i].lcps;
        refs = srt_data[i].refs;
        if (refs) { refs += 1; refs[-1] = -1; } // for unmapped

        for (j = 0; j < nfile; j++) {
            // for input streams        
            switch (ftype) {
                case FT_ZSTD:
                    bams[n].reader = zst_read_data;
                    bams[n].freeArgs = zst_readerArgs_free;
                    bams[n].readerArgs = zst_readerArgs_init((FILE *) fns[j], 0);
                    if (!bams[n].readerArgs) {
                        errmsg("failed to create ZSTD reader extra arguments");
                        goto err;
                    }
                    break;
                case FT_BAM:
                case FT_SAM:
                    bams[n].reader = sam_read_data;
                    bams[n].freeArgs = sam_readerArgs_free;
                    bams[n].readerArgs = sam_readerArgs_init((samFile *) fns[j], hdr, 
                        (ftype == FT_BAM)? bam_read1_lite : sam_read1_lite);
                    if (!bams[n].readerArgs) {
                        errmsg("failed to create SAM reader extra arguments");
                        goto err;
                    }
                    break;
                default:
                    errmsg("unsupported file type for BAM stream");
                    goto err;
            }
            // other initialisations
            bams[n].outBufSize = max_mem;
            bams[n].outBuf = (uint8_t *) malloc(max_mem);
            bams[n].recs[0] = (bam1_t **) malloc(sizeof(bam1_t *) * max_rec * N_CHNK);
            bams[n].lcps[0] = (int *) malloc(sizeof(int) * max_rec * N_CHNK);
            if (!bams[n].outBuf || !bams[n].recs[0] || !bams[n].lcps[0]) {
                errmsg("failed to allocate memory for BAM records");
                goto err;
            }
            bams[n].data[0] = bams[n].outBuf;
            for (k = 1; k < N_CHNK; k++) {
                bams[n].data[k] = bams[n].outBuf;
                bams[n].recs[k] = bams[n].recs[k-1] + max_rec;
                bams[n].lcps[k] = bams[n].lcps[k-1] + max_rec;
            }
            bams[n].lcp = lcps? lcps[j] : 0;
            bams[n].block = n;
            bams[n].mrec = max_rec;
            bams[n].refmap = refs;

            n++;
        }
    }
    if (_n_bs) *_n_bs = n_bs;

    return bams;

err:
    for (i = 0; i < n_bs; i++) {
        if (bams[i].freeArgs)
            bams[i].freeArgs(bams[i].readerArgs);
        free(bams[i].outBuf);
        free(bams[i].recs[0]);
        free(bams[i].lcps[0]);
    }
    free(bams);

    return NULL;
}

static int sam_hdr_lite_trim(sam_hdr_lite_t *hdr, int *refs);
static void map_refseq(bio_data_t *bio, size_t n_bio, int n_threads, int32_t *refmap);
static void set_s_rlimit();
static bool open_sorted_file(char *infile, srt_data_t *srt, refseq_builder_t *refseqs);
bool mergebam(char **infiles, int nIn, char *outfile, size_t max_mem, int n_threads, bool zstOutput);

bool bamsort(char **infiles, int nIn, char *outfile, size_t max_mem, int n_threads, 
    bool noTrimHeader, bool zstOutput, bool mergeOnly)
{
    if (mergeOnly) return mergebam(infiles, nIn, outfile, max_mem, n_threads, zstOutput);

    bio_data_t *bio_data = NULL;
    srt_data_t *srt_data = NULL;
    zstd_writer_t *zst_writer = NULL;
    sam_hdr_lite_t *header = NULL;
    sam_hdr_lite_t *builder = NULL;
    bam1_stream_t *bam_streams = NULL;
    bam_writer_t *bam_writer = NULL;
    refseq_builder_t *refseqs = NULL;
    struct stat st;
    size_t name_len, max_mem2, max_rec2, tot_mem, n_recs;
    char *fn_template = NULL;
    FILE *fn_out = NULL;
    bool ret = false, derived = false, bam_ok = true;
    bool trim = !noTrimHeader;
    int i, j, *refs, n_part, n_refs, n_streams = 0;

    // set endianness
    long one = 1;
    IS_BIG_ENDIAN = !(*((char *)(&one)));
    
    // sanity checks
    // input files
    if (infiles == NULL || nIn < 1) {
        errmsg("no input files specified");
        goto err;
    }
    for (i = 0; i < nIn; i++) {
        if (infiles[i] == NULL ||
            access(infiles[i], R_OK) ||
            stat(infiles[i], &st) ||
            !S_ISREG(st.st_mode)) {
            errmsg("accessing input file failed: %s", infiles[i]);
            goto err;
        }
    }
    // set resource limits
    set_s_rlimit();

    // total threads
    if (n_threads < 2) n_threads = 1;
    if (n_threads < W_THREADS) W_THREADS = n_threads;
    N_THREADS = n_threads;
    // memory
    if (!max_mem) max_mem = MAX_MEM_PER_THREAD;
    if (max_mem < (SORT_MIN_MEGS_PER_THREAD << 20)) {
        complain_about_memory_setting(max_mem);
        goto err;
    }
    tot_mem = max_mem * N_THREADS;
    // we will open n_part = n_threads*2 partitions
    max_mem2 = max_mem >> 1;
    // estimate number of BAM records that can fit in memory
    // but will be capped by max_mem2
    max_rec2 = max_mem2 / BAM_REC_SIZE;
    infomsg("number parts per file: %d", N_THREADS * 2);
    infomsg("maximum mems per part: %zu Mb", max_mem2>>20);
    infomsg("maximum recs per part: %d", max_rec2);

    // output file
    if (!outfile) {
        outfile = derivedName(infiles[0], "sorted.bam");
        if (!outfile) {
            errmsg("failed to create output file name");
            return false;
        }
        derived = true;
        if (stat(outfile, &st) == 0) {
            errmsg("output file exists: %s", outfile);
            goto err;
        }
    }

    // allocate memory for sorting
    S_THREADS = N_THREADS;
    n_part = S_THREADS * 2;
    bio_data = (bio_data_t *) calloc(n_part, sizeof(bio_data_t));
    if (!bio_data) {
        errmsg("memory allocation failed");
        goto err;
    }
    max_mem2 = max_mem >> 1;
    for (i = 0; i < n_part; i++) {
        bio_data[i].data = (uint8_t *) malloc(max_mem2);
        bio_data[i].recs = (bam1_t **) malloc(sizeof(bam1_t *) * (max_rec2+1));
        bio_data[i].lcps = (int *) malloc(sizeof(int) * (max_rec2+1));
        bio_data[i].b = bam_init1();
        if (!bio_data[i].data || !bio_data[i].recs || 
            !bio_data[i].lcps || !bio_data[i].b) {
            errmsg("memory allocation failed");
            goto err;
        }
        bio_data[i].mrec = max_rec2;
        bio_data[i].mmem = max_mem2;
    }
    srt_data = (srt_data_t *) calloc(nIn, sizeof(srt_data_t));
    if (!srt_data) {
        errmsg("memory allocation failed");
        goto err;
    }
    refseqs = refseq_builder_init(nIn, 0, NULL);
    if (!refseqs) {
        errmsg("failed to initialise reference sequence builder");
        goto err;
    }

    // housekeeping for sorting threads
    Tarray = (pthread_t *) malloc(S_THREADS * sizeof(pthread_t));
    Tjoins = (bool *) malloc(S_THREADS * sizeof(bool));
    Tstack  = (int *) malloc(S_THREADS * sizeof(int));
    SJlist  = (int *) malloc(n_part * sizeof(int));
    WJlist  = (int *) malloc(n_part * sizeof(int));
    SJdone  = (bool *) calloc(n_part, sizeof(bool));
    if (!Tarray || !Tjoins || !Tstack || !SJlist || !WJlist || !SJdone) {
        errmsg("threads memory allocation failed");
        goto err;
    }
    pthread_mutex_init(&JMUTEX,NULL);
    pthread_cond_init(&JCOND,NULL);
    pthread_mutex_init(&TMUTEX,NULL);
    pthread_cond_init(&TCOND,NULL);
    pthread_mutex_init(&WMUTEX,NULL);
    pthread_cond_init(&WCOND,NULL);

    // set output file format
    name_len = strlen(outfile) + 64;
    fn_template = (char *) malloc(name_len);
    if (!fn_template) {
        errmsg("file name buffer allocation failed");
        goto err;
    }

    // initialize zstd writer if needed
    zst_writer = zstd_writer_init(ZSTD_COMPRESSION_LEVEL, W_THREADS, NULL);
    if (!zst_writer) {
        errmsg("failed to initialise zstd writer");
        goto err;
    }

    // sort each input file
    bool in_mem = false, in_mem_check = nIn == 1, multi_in = nIn > 1;
    if (multi_in) {
        // initialise header builder
        builder = sam_hdr_lite_builder_init();
        if (!builder) {
            errmsg("failed to initialize header builder");
            goto err;
        }
    }
    for (i = 0; i < nIn; i++) {
        infomsg("sorting input file [%d/%d]: %s", i+1, nIn, infiles[i]);
        snprintf(fn_template, name_len, "%s.%.6d.XXXXXX", outfile, i+1);
        if (!bamsort1(infiles[i], fn_template, trim, multi_in, bio_data, 
            srt_data+i, refseqs, in_mem_check, zst_writer, zst_merge)) {
            errmsg("failed to sort input file: %s", infiles[i]);
            goto err;
        }
        srt_data[i].ftype = FT_ZSTD; // use zst_writer for merging
        if (in_mem_check && srt_data[i].nfile == 0)
            in_mem = true;
        in_mem_check = false; // only check for the first file
        if (multi_in) {
            // build reference translation table
            push_time(&TIME_STACK);
            header = srt_data[i].header;
            refs = srt_data[i].refs;
            if (refs) { refs += 1; refs[-1] = -1; } // for unmapped 
            if (refs) {
                n_refs = 0;
                for (j = 0; j < header->n_targets; j++)
                    n_refs += refs[j] > 0;
            } else n_refs = header->n_targets;
            if (!sam_hdr_lite_build_tbl(builder, header, refs)) {
                errmsg("failed to build translation table");
                goto err;
            }
            
            infomsg("built translation table in %.3f seconds", pop_time(&TIME_STACK));
            infomsg("BAM file reference sequences: %d", header->n_targets);
            infomsg("retained reference sequences: %d", n_refs);

            // reclaim some memory
            sam_hdr_lite_destroy(header); srt_data[i].header = NULL;
        }
    }

    // clean up sorting thread resources
    pthread_mutex_destroy(&JMUTEX);
    pthread_cond_destroy(&JCOND);
    pthread_mutex_destroy(&TMUTEX);
    pthread_cond_destroy(&TCOND);
    pthread_mutex_destroy(&WMUTEX);
    pthread_cond_destroy(&WCOND);
    free(Tarray), Tarray = NULL;
    free(Tjoins), Tjoins = NULL;
    free(Tstack), Tstack = NULL;
    free(SJlist), SJlist = NULL;
    free(WJlist), WJlist = NULL;
    free(SJdone), SJdone = NULL;

    if (multi_in) {
        header = builder; builder = NULL;
    } else {
        // only one header
        header = srt_data[0].header;
        srt_data[0].header = NULL; // do not destroy
        if (trim) {
            refs = refseqs->refs[0] + 1; // offset for unmapped
            if (sam_hdr_lite_trim(header, refs) < 0) {
                errmsg("failed to trim SAM header");
                goto err;
            }
            if (in_mem) { // in-memory translation
                map_refseq(bio_data, n_part, n_threads, refs);
                free(refseqs->refs[0]); refseqs->refs[0] = NULL;
            }
        }
    }
    if (sam_hdr_lite_build_HD(header, NULL, "queryname") < 0) {
        errmsg("failed to build SAM header HD line");
        goto err;
    }
    
    // reclaim memory for sorting
    if (!in_mem) { // only needed for in-memory merging
        for (i = 0; i < n_part; i++) {
            free(bio_data[i].data); bio_data[i].data = NULL;
            free(bio_data[i].recs); bio_data[i].recs = NULL;
            free(bio_data[i].lcps); bio_data[i].lcps = NULL;
            bam_destroy1(bio_data[i].b); bio_data[i].b = NULL;
        }
        free(bio_data); bio_data = NULL;
    }

    // do not need these anymore
    if (header->sq_map) { kh_destroy(c2i,  header->sq_map); header->sq_map = NULL; }
    if (header->rg_ids) { kh_destroy(cset, header->rg_ids); header->rg_ids = NULL; }
    if (header->pg_ids) { kh_destroy(cset, header->pg_ids); header->pg_ids = NULL; }

    if (zstOutput) {
        // open output file
        fn_out = fopen(outfile, "w+b");
        if (!fn_out) {
            errmsg("failed to open output file: %s", outfile);
            goto err;
        }
        zst_writer->fp = fn_out;
        // write header to output file
        if (zstd_sam_hdr_lite_write(zst_writer, header) != 0) {
            errmsg("failed to write SAM header to ZSTD file \"%s\"", outfile);
            goto err;
        }
    } else {
        // distribute threads accross reading and writing
        int r_threads, w_threads;
        r_threads = w_threads = N_THREADS / 2; // 50% for writing
        if (r_threads > 9) {
            r_threads = 9;
            w_threads = N_THREADS - r_threads;
        }
        if (w_threads < 1) w_threads = 1;
        if (w_threads > MAX_WWORKER_THREADS) w_threads = MAX_WWORKER_THREADS;
        r_threads = N_THREADS - w_threads;
        if (r_threads < 1) r_threads = 1;
        N_THREADS = r_threads;
        W_THREADS = 1; // this is actually the main thread
        infomsg("allocated %d threads for reading", r_threads);
        infomsg("allocated %d threads for writing", w_threads);

        size_t l_text = sam_hdr_lite_min_text_len(header);
        bam_ok = l_text <= BAM_MAX_HEADER_SIZE;
        // open output file
        if (!bam_ok) {
            // change output file name
            char *f = (char *) malloc(name_len);
            int l = strlen(outfile), x = 0;
            pid_t pid = getpid();
            if (l > 4 && strcmp(outfile + l - 4, ".bam") == 0)
                l -= 4;
            snprintf(f, name_len, "%.*s.sam.gz", l, outfile);
            while (stat(f, &st) == 0) {
                // file exists
                // do not give up yet and try to create a new name
                snprintf(f, name_len, "%.*s.%d-%d.sam.gz", l, outfile, pid, ++x);
            }
            if (derived) free(outfile);
            outfile = f;
            derived = true;
            warnmsg("BAM header too large writing as BGZF compressed SAM: \"%s\"", outfile);
        }
        fn_out = fopen(outfile, "w+b");
        if (!fn_out) {
            errmsg("failed to open output file: %s", outfile);
            goto err;
        }
        // write header to output file
        if (bam_ok) {
            if (bw_bam_hdr_lite_write(fn_out, header) != 0) {
                errmsg("failed to write SAM header to BAM file \"%s\"", outfile);
                goto err;
            }
        } else {
            if (bw_sam_hdr_lite_write(fn_out, header) != 0) {
                errmsg("failed to write SAM header to SAM file \"%s\"", outfile);
                goto err;
            }
        }

        bam_writer = bam_writer_init(fn_out, w_threads, WW_MAX_MEM);
        if (!bam_writer) {
            errmsg("failed to initialize multi-threaded BAM writer");
            goto err;
        }
        bam_writer->compress = bam_ok? bw_bam_compress : bw_sam_compress;
    }
    
    // for sam output we need header
    if (!bam_ok) {
        // (re)build target_name array
        // the internal-nulls are removed when writing SAM header
        sam_hdr_lite_build_targets_tbl(header);
        for (i = 0; i < bam_writer->n_workers; i++) {
            bam_writer->workers[i].args = sam_compArgs_make(header, 1024);
            bam_writer->workers[i].argfree = sam_compArgs_free;
        }
    } else cstr_array_destroy(&header->sq_text);
    // header has been written to output file
    // destroy text buffers to reclaim memory    
    ks_free(&header->hd_text);
    cstr_array_destroy(&header->rg_text);
    cstr_array_destroy(&header->pg_text);
    cstr_array_destroy(&header->co_text);

    n_recs = 0;
    for (i = 0; i < nIn; i++)
        n_recs += srt_data[i].nrec;
    
    // merge sorted files
    if (in_mem) {
        infomsg("write in-memory data to %s: %zu records", outfile, n_recs);
        // merge and write to output file
        merge_blocks_in_mem(bio_data, n_part, NULL,
            zstOutput? (void *)zst_writer : (void *)bam_writer,
            zstOutput? zst_merge : bam_merge
        );
    } else {
        infomsg("merging sorted files to %s: %zu records", outfile, n_recs);
        // collect streams from all temporary files
        bam_streams = create_bam1_streams(srt_data, nIn, tot_mem, &n_streams);
        if (!bam_streams) {
            errmsg("failed to create BAM stream contexts");
            goto err;
        }
        // merge and write to output file
        merge_blocks(bam_streams, n_streams, 
            zstOutput? (void *)zst_writer : (void *)bam_writer,
            zstOutput? zst_merge_thread : bam_merge_thread
        );
    }

    if (zstOutput)
        fclose(fn_out), fn_out = NULL;
    else if (bam_writer_close(bam_writer) != 0) {
        errmsg("failed to finalize BAM output file");
        goto err;
    }

    ret = true;

err:
    if (derived) free(outfile);
    if (fn_out) fclose(fn_out);
    if (bio_data) {
        for (i = 0; i < n_part; i++) {
            free(bio_data[i].data);
            free(bio_data[i].recs);
            free(bio_data[i].lcps);
            if (bio_data[i].ftype == FT_ZSTD)
                fclose((FILE *)bio_data[i].fp);
            else
                sam_close((samFile *)bio_data[i].fp);
            bam_destroy1(bio_data[i].b);
        }
        free(bio_data);
    }
    if (srt_data) {
        for (i = 0; i < nIn; i++) {
            for (j = 0; j < srt_data[i].nfile; j++)
                if (srt_data[i].ftype == FT_ZSTD)
                    fclose((FILE *)srt_data[i].fns[j]);
                else
                    sam_close((samFile *)srt_data[i].fns[j]);
            free(srt_data[i].fns);
            free(srt_data[i].lcps);
            sam_hdr_lite_destroy(srt_data[i].header);
        }
        free(srt_data);
    }
    if (bam_streams) {
        for (i = 0; i < n_streams; i++) {
            if (bam_streams[i].freeArgs)
                bam_streams[i].freeArgs(bam_streams[i].readerArgs);
            free(bam_streams[i].outBuf);
            free(bam_streams[i].recs[0]);
            free(bam_streams[i].lcps[0]);
        }
        free(bam_streams);
    }
    free(fn_template);
    refseq_builder_destroy(refseqs);

    sam_hdr_lite_destroy(header);
    zstd_writer_destroy(zst_writer);
    bam_writer_destroy(bam_writer);

    pthread_mutex_destroy(&JMUTEX);
    pthread_cond_destroy(&JCOND);
    pthread_mutex_destroy(&TMUTEX);
    pthread_cond_destroy(&TCOND);
    pthread_mutex_destroy(&WMUTEX);
    pthread_cond_destroy(&WCOND);

    free(Tarray), Tarray = NULL;
    free(Tjoins), Tjoins = NULL;
    free(Tstack), Tstack = NULL;
    free(SJlist), SJlist = NULL;
    free(WJlist), WJlist = NULL;
    free(SJdone), SJdone = NULL;

    return ret;
}

bool mergebam(char **infiles, int nIn, char *outfile, size_t max_mem, int n_threads, bool zstOutput)
{
    srt_data_t *srt_data = NULL;
    sam_hdr_lite_t *header = NULL;
    sam_hdr_lite_t *builder = NULL;
    bam1_stream_t *bam_streams = NULL;
    bam_writer_t *bam_writer = NULL;
    zstd_writer_t *zst_writer = NULL;
    refseq_builder_t *refseqs = NULL;
    struct stat st;
    size_t name_len, tot_mem;
    FILE *fn_out = NULL;
    bool ret = false, derived = false, bam_ok = true;
    int i, j, *refs, n_part, n_refs, n_streams = 0;

    // set endianness
    long one = 1;
    IS_BIG_ENDIAN = !(*((char *)(&one)));
    
    // sanity checks
    // input files
    if (infiles == NULL || nIn < 1) {
        errmsg("no input files specified");
        goto err;
    }
    if (nIn < 2) {
        errmsg("need at least two input files to merge");
        goto err;
    }
    for (i = 0; i < nIn; i++) {
        if (infiles[i] == NULL ||
            access(infiles[i], R_OK) ||
            stat(infiles[i], &st) ||
            !S_ISREG(st.st_mode)) {
            errmsg("accessing input file failed: %s", infiles[i]);
            goto err;
        }
    }
    // set resource limits
    set_s_rlimit();

    // total threads
    if (n_threads < 2) n_threads = 1;
    if (n_threads < W_THREADS) W_THREADS = n_threads;
    N_THREADS = n_threads;
    // memory
    if (!max_mem) max_mem = MAX_MEM_PER_THREAD;
    if (max_mem < (SORT_MIN_MEGS_PER_THREAD << 20)) {
        complain_about_memory_setting(max_mem);
        goto err;
    }
    tot_mem = max_mem * N_THREADS;

    // output file
    if (!outfile) {
        outfile = derivedName(infiles[0], "merged.bam");
        if (!outfile) {
            errmsg("failed to create output file name");
            return false;
        }
        derived = true;
        if (stat(outfile, &st) == 0) {
            errmsg("output file exists: %s", outfile);
            goto err;
        }
    }
    // set output file format
    name_len = strlen(outfile) + 64;

    // allocate memory for inputs
    srt_data = (srt_data_t *) calloc(nIn, sizeof(srt_data_t));
    if (!srt_data) {
        errmsg("memory allocation failed");
        goto err;
    }
    refseqs = refseq_builder_init(nIn, 0, NULL);
    if (!refseqs) {
        errmsg("failed to initialize reference sequence builder");
        goto err;
    }

    // initialise header builder
    builder = sam_hdr_lite_builder_init();
    if (!builder) {
        errmsg("failed to initialize header builder");
        goto err;
    }
    for (i = 0; i < nIn; i++) {
        infomsg("reading header from input file [%d/%d]: %s", i+1, nIn, infiles[i]);
        // build reference translation table
        push_time(&TIME_STACK);
        if (!open_sorted_file(infiles[i], srt_data+i, refseqs)) {
            errmsg("failed to read header from file: %s", infiles[i]);
            goto err;
        }
        header = srt_data[i].header;
        refs = srt_data[i].refs;
        if (refs) { 
            refs += 1;
            // we do not do trimming for merging
            /**refs[-1] = -1;**/ // for unmapped
            // refs[-1] also used as a flag for header duplication
            // if we have refs[-1] > 0 we save time for rebuilding
            if (refs[-1] <= 0) {
                // build full mapping
                for (j = 0; j < header->n_targets; j++)
                    refs[j] = 1;
            }
        }
        if (!sam_hdr_lite_build_tbl(builder, header, refs)) {
            errmsg("failed to build translation table");
            goto err;
        }
        infomsg("built translation table in %.3f seconds", pop_time(&TIME_STACK));
        infomsg("BAM file reference sequences: %d", header->n_targets);
        infomsg("retained reference sequences: %d", builder->n_targets);

        // reclaim some memory unless needed for reading
        if (srt_data[i].ftype != FT_SAM) {
            sam_hdr_lite_destroy(header);
            srt_data[i].header = NULL;
        }
    }
    if (sam_hdr_lite_build_HD(builder, NULL, "queryname") < 0) {
        errmsg("failed to build SAM header HD line");
        goto err;
    }
    header = builder; builder = NULL;

    // do not need these anymore
    if (header->sq_map) { kh_destroy(c2i,  header->sq_map); header->sq_map = NULL; }
    if (header->rg_ids) { kh_destroy(cset, header->rg_ids); header->rg_ids = NULL; }
    if (header->pg_ids) { kh_destroy(cset, header->pg_ids); header->pg_ids = NULL; }

    if (zstOutput) {
        // open output file
        fn_out = fopen(outfile, "w+b");
        if (!fn_out) {
            errmsg("failed to open output file: %s", outfile);
            goto err;
        }
        // initialize zstd writer
        zst_writer = zstd_writer_init(ZSTD_COMPRESSION_LEVEL, W_THREADS, NULL);
        if (!zst_writer) {
            errmsg("failed to initialize zstd writer");
            goto err;
        }
        zst_writer->fp = fn_out;
        // write header to output file
        if (zstd_sam_hdr_lite_write(zst_writer, header) != 0) {
            errmsg("failed to write SAM header to ZSTD file \"%s\"", outfile);
            goto err;
        }
    } else {
        // distribute threads accross reading and writing
        int r_threads, w_threads;
        r_threads = w_threads = N_THREADS / 2; // 50% for writing
        if (r_threads > 9) {
            r_threads = 9;
            w_threads = N_THREADS - r_threads;
        }
        if (w_threads < 1) w_threads = 1;
        if (w_threads > MAX_WWORKER_THREADS) w_threads = MAX_WWORKER_THREADS;
        r_threads = N_THREADS - w_threads;
        if (r_threads < 1) r_threads = 1;
        N_THREADS = r_threads;
        W_THREADS = 1; // this is actually the main thread
        infomsg("allocated %d threads for reading", r_threads);
        infomsg("allocated %d threads for writing", w_threads);

        size_t l_text = sam_hdr_lite_min_text_len(header);
        bam_ok = l_text <= BAM_MAX_HEADER_SIZE;
        // open output file
        if (!bam_ok) {
            // change output file name
            char *f = (char *) malloc(name_len);
            int l = strlen(outfile), x = 0;
            pid_t pid = getpid();
            if (l > 4 && strcmp(outfile + l - 4, ".bam") == 0)
                l -= 4;
            snprintf(f, name_len, "%.*s.sam.gz", l, outfile);
            while (stat(f, &st) == 0) {
                // file exists
                // do not give up yet and try to create a new name
                snprintf(f, name_len, "%.*s.%d-%d.sam.gz", l, outfile, pid, ++x);
            }
            if (derived) free(outfile);
            outfile = f;
            derived = true;
            warnmsg("BAM header too large writing as BGZF compressed SAM: \"%s\"", outfile);
        }
        fn_out = fopen(outfile, "w+b");
        if (!fn_out) {
            errmsg("failed to open output file: %s", outfile);
            goto err;
        }
        // write header to output file
        if (bam_ok) {
            if (bw_bam_hdr_lite_write(fn_out, header) != 0) {
                errmsg("failed to write SAM header to BAM file \"%s\"", outfile);
                goto err;
            }
        } else {
            if (bw_sam_hdr_lite_write(fn_out, header) != 0) {
                errmsg("failed to write SAM header to SAM file \"%s\"", outfile);
                goto err;
            }
        }

        bam_writer = bam_writer_init(fn_out, w_threads, WW_MAX_MEM);
        if (!bam_writer) {
            errmsg("failed to initialize multi-threaded BAM writer");
            goto err;
        }
        bam_writer->compress = bam_ok? bw_bam_compress : bw_sam_compress;
    }
    
    // for sam output we need header
    if (!bam_ok) {
        // (re)build target_name array
        // the internal-nulls are removed when writing SAM header
        sam_hdr_lite_build_targets_tbl(header);
        for (i = 0; i < bam_writer->n_workers; i++) {
            bam_writer->workers[i].args = sam_compArgs_make(header, 1024);
            bam_writer->workers[i].argfree = sam_compArgs_free;
        }
    } else cstr_array_destroy(&header->sq_text);
    // header has been written to output file
    // destroy text buffers to reclaim memory    
    ks_free(&header->hd_text);
    cstr_array_destroy(&header->rg_text);
    cstr_array_destroy(&header->pg_text);
    cstr_array_destroy(&header->co_text);

    infomsg("merging sorted files to %s", outfile);
    // collect streams from all temporary files
    bam_streams = create_bam1_streams(srt_data, nIn, tot_mem, &n_streams);
    if (!bam_streams) {
        errmsg("failed to create BAM stream contexts");
        goto err;
    }
    // merge and write to output file
    merge_blocks(bam_streams, n_streams, 
        zstOutput? (void *)zst_writer : (void *)bam_writer,
        zstOutput? zst_merge_thread : bam_merge_thread
    );

    if (zstOutput)
        fclose(fn_out), fn_out = NULL;
    else if (bam_writer_close(bam_writer) != 0) {
        errmsg("failed to finalize BAM output file");
        goto err;
    }

    ret = true;

err:
    if (derived) free(outfile);
    if (fn_out) fclose(fn_out);
    if (srt_data) {
        for (i = 0; i < nIn; i++) {
            for (j = 0; j < srt_data[i].nfile; j++)
                if (srt_data[i].ftype == FT_ZSTD)
                    fclose((FILE *)srt_data[i].fns[j]);
                else
                    sam_close((samFile *)srt_data[i].fns[j]);
            free(srt_data[i].fns);
            free(srt_data[i].lcps);
            sam_hdr_lite_destroy(srt_data[i].header);
        }
        free(srt_data);
    }
    if (bam_streams) {
        for (i = 0; i < n_streams; i++) {
            if (bam_streams[i].freeArgs)
                bam_streams[i].freeArgs(bam_streams[i].readerArgs);
            free(bam_streams[i].outBuf);
            free(bam_streams[i].recs[0]);
            free(bam_streams[i].lcps[0]);
        }
        free(bam_streams);
    }

    refseq_builder_destroy(refseqs);
    sam_hdr_lite_destroy(header);
    zstd_writer_destroy(zst_writer);
    bam_writer_destroy(bam_writer);

    return ret;
}

static bool open_sorted_file(char *infile, srt_data_t *srt, refseq_builder_t *rdb)
{
    if (!infile || !srt || !rdb) return false;
    // initialise srt_data_t
    memset(srt, 0, sizeof(*srt));
    srt->nfile = 1; // single file
    srt->fns = (void **)calloc(1, sizeof(void *));
    if (!srt->fns) {
        errmsg("open_sorted_file: memory alloc failed");
        return false;
    }

    int r = 0, *refs = NULL;
    
    // Try ZSTD first (magic: 0xFD2FB528 little-endian -> bytes 28 B5 2F FD)
    FILE *fz = fopen(infile, "rb");
    if (!fz) {
        errmsg("open_sorted_file: failed to open as ZSTD %s", infile);
        return false;
    }
    srt->fns[0] = (void *)fz;
    srt->ftype = FT_ZSTD;

    unsigned char head[4];
    size_t bytes = fread(head, 1, 4, fz);
    if (bytes == 4 && head[0] == 0x28 && head[1] == 0xB5 && head[2] == 0x2F && head[3] == 0xFD) {
        // rewind
        if (fseek(fz, 0, SEEK_SET) != 0) {
            errmsg("open_sorted_file: fseek failed %s", infile);
            return false;
        }
        srt->header = zstd_sam_hdr_lite_read(fz);
        if (!srt->header) {
            errmsg("open_sorted_file: failed to read header %s", infile);
            return false;
        }
        if ((r = refseq_builder_add(rdb, srt->header, &refs)) < 0) {
            errmsg("failed to build reference sequence map");
            return false;
        }
        refs[0] = r; // this is a indicator if it is a duplication
        srt->refs = refs;
        return true;
    }
    // Not ZSTD; close and use htslib sam/bam path
    fclose(fz);

    samFile *fs = sam_open(infile, "r");
    if (!fs) {
        errmsg("open_sorted_file: failed to open as SAM/BAM %s", infile);
        return false;
    }
    srt->fns[0] = (void *)fs;
    if (fs->format.format == bam) 
        srt->ftype = FT_BAM;
    else if (fs->format.format == sam)
        srt->ftype = FT_SAM;
    else {
        errmsg("open_sorted_file: unsupported format %s", infile);
        return false;
    }

    srt->header = sam_hdr_lite_read(fs);
    if (!srt->header) {
        errmsg("open_sorted_file: failed to read header %s", infile);
        return false;
    }
    if ((r = refseq_builder_add(rdb, srt->header, &refs)) < 0) {
        errmsg("failed to build reference sequence map");
        return false;
    }
    refs[0] = r; // this is a indicator if it is a duplication
    srt->refs = refs;

    return true;
}

/**********************************************/
/********************** misc ******************/
/**********************************************/

/************************ message *****************/
static void errmsg(char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "[E::%s] ", PROGRAM) ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, ". Aborting.\n") ;
  va_end (args) ;
}

static void infomsg(char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "[M::%s] ", PROGRAM) ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
}

static void warnmsg(char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "[W::%s] ", PROGRAM) ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
}

static void dbgmsg(char *format, ...)
{
  va_list args ;

  va_start (args, format) ;
  fprintf (stderr, "[D::%s] ", PROGRAM) ;
  vfprintf (stderr, format, args) ;
  fprintf (stderr, "\n") ;
  va_end (args) ;
}

static void set_s_rlimit()
{
    struct rlimit rl;
    
    if (getrlimit(RLIMIT_NOFILE, &rl) != 0) {
        warnmsg("getrlimit failed");
        return;
    }

    unsigned long long s_lim = rl.rlim_cur;
    if (s_lim < rl.rlim_max) {
        rl.rlim_cur = rl.rlim_max;
        if(setrlimit(RLIMIT_NOFILE, &rl) != 0)
            warnmsg("setrlimit failed");
    }
    return;
}

/************************ timer stack *****************/

void push_time(timer1_t **timer) 
{
    timer1_t *lap = (timer1_t *) malloc(sizeof(timer1_t));
    if (!lap) return;
    gettimeofday(&lap->tv, NULL);
    lap->next = *timer;
    *timer = lap;
}

// get the elapsed time since the last push and pop the timer stack
double pop_time(timer1_t **timer) 
{
    if (!timer || !*timer) {
        warnmsg("timer stack underflow");
        return -1.0;
    }
    struct timeval end, start;
    gettimeofday(&end, NULL);

    timer1_t *lap = *timer;
    *timer = (*timer)->next;
    start = lap->tv;
    free(lap);

    return (end.tv_sec - start.tv_sec) +
           (end.tv_usec - start.tv_usec) / 1e6;
}

// get the elapsed time since the last push
double elapsed_time(timer1_t **timer, bool reset) 
{
    if (!timer || !*timer) {
        warnmsg("timer stack underflow");
        return -1.0;
    }
    struct timeval end, start;
    gettimeofday(&end, NULL);
    start = (*timer)->tv;

    if (reset) (*timer)->tv = end;

    return (end.tv_sec - start.tv_sec) +
           (end.tv_usec - start.tv_usec) / 1e6;
}

void free_timer(timer1_t *timer) 
{
    while (timer) {
        timer1_t *next = timer->next;
        free(timer);
        timer = next;
    }
}


/************************ ZSTD ************************/
static zstd_cctx_t *zstd_cctx_init(int compressionLevel, int threads) {
    zstd_cctx_t *zstd_cctx = NULL;
    ZSTD_CCtx *cctx = NULL;
    size_t inBufSize, outBufSize;
    char *inBuf = NULL, *outBuf = NULL;

    if (threads < 1) threads = 1;
    if (threads > ZSTD_MAX_THREADS) threads = ZSTD_MAX_THREADS;

    zstd_cctx = (zstd_cctx_t *) malloc(sizeof(zstd_cctx_t));
    if (!zstd_cctx) goto err;
    cctx = ZSTD_createCCtx();
    if (!cctx) goto err;
    
    size_t ret = 0;
    ret = ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, compressionLevel);
    if (ZSTD_isError(ret)) 
        warnmsg("ZSTD failed to set compression level: %s", ZSTD_getErrorName(ret));
    ret = ZSTD_CCtx_setParameter(cctx, ZSTD_c_nbWorkers, threads);
    if (ZSTD_isError(ret)) 
        warnmsg("ZSTD failed to set number of workers: %s", ZSTD_getErrorName(ret));
    ret = ZSTD_CCtx_setParameter(cctx, ZSTD_c_jobSize, ZSTD_INBUF_SIZE);
    if (ZSTD_isError(ret)) 
        warnmsg("ZSTD failed to set job size: %s", ZSTD_getErrorName(ret));

    inBufSize = ZSTD_INBUF_SIZE * threads;
    outBufSize = ZSTD_CStreamOutSize();
    inBuf = (char *) malloc(inBufSize);
    outBuf = (char *) malloc(outBufSize);
    if (!inBuf || !outBuf) goto err;

    zstd_cctx->cctx = cctx;
    zstd_cctx->inBufSize = inBufSize;
    zstd_cctx->outBufSize = outBufSize;
    zstd_cctx->inBuf = inBuf;
    zstd_cctx->outBuf = outBuf;
    return zstd_cctx;

err:
    free(zstd_cctx);
    ZSTD_freeCCtx(cctx);
    free(inBuf);
    free(outBuf);
    return NULL;
}

static void zstd_cctx_destroy(zstd_cctx_t *zstd_cctx) {
    if (!zstd_cctx) return;
    ZSTD_freeCCtx(zstd_cctx->cctx);
    free(zstd_cctx->inBuf);
    free(zstd_cctx->outBuf);
    free(zstd_cctx);
}

static zstd_dctx_t *zstd_dctx_init(size_t inBufSize, size_t outBufSize)
{
    zstd_dctx_t *zstd_dctx = NULL;
    ZSTD_DCtx *dctx = NULL;
    char *inBuf = NULL, *outBuf = NULL;

    if (!inBufSize) inBufSize = ZSTD_INBUF_SIZE;
    if (!outBufSize) outBufSize = ZSTD_OUTBUF_SIZE;

    zstd_dctx = (zstd_dctx_t *) malloc(sizeof(zstd_dctx_t));
    if (!zstd_dctx) goto err;

    dctx = ZSTD_createDCtx();
    if (!dctx) goto err;

    inBuf = (char *) malloc(inBufSize);
    outBuf = (char *) malloc(outBufSize);
    if (!inBuf || !outBuf) goto err;

    zstd_dctx->dctx = dctx;
    zstd_dctx->inBuf = (ZSTD_inBuffer) { inBuf, 0, 0 };
    zstd_dctx->outBufSize = outBufSize;
    zstd_dctx->outBuf = outBuf;

    return zstd_dctx;

err:
    free(zstd_dctx);
    ZSTD_freeDCtx(dctx);
    free(inBuf);
    free(outBuf);
    return NULL;
}

static void zstd_dctx_destroy(zstd_dctx_t *zstd_dctx) 
{
    if (!zstd_dctx) return;
    ZSTD_freeDCtx(zstd_dctx->dctx);
    free((void *)zstd_dctx->inBuf.src);
    free(zstd_dctx->outBuf);
    free(zstd_dctx);
}


static inline bool write_buffered_data(char *inBuf, size_t inSize, char *outBuf, size_t outBufSize, ZSTD_CCtx *cctx, FILE *fp)
{
    ZSTD_inBuffer input = { inBuf, inSize, 0 };
    ZSTD_outBuffer output = { outBuf, outBufSize, 0 };
    size_t ret;
    while (input.pos < input.size) {
        output.pos = 0;
        ret = ZSTD_compressStream2(cctx, &output, &input, ZSTD_e_continue);
        if (ZSTD_isError(ret)) {
            errmsg("zstd compression error: %s\n", ZSTD_getErrorName(ret));
            return false;
        }
        if (output.pos > 0) {
            if (fwrite(outBuf, output.pos, 1, fp) != 1) {
                errmsg("failed to write zstd-compressed data to file");
                return false;
            }
        }
    }

    return true;
}

static inline int zstd_writer_flush(zstd_writer_t *writer)
{
    ZSTD_CCtx *cctx = writer->cctx;
    size_t inBufSize = writer->inBufSize, outBufSize = writer->outBufSize;
    char *inBuf = writer->inBuf, *outBuf = writer->outBuf;
    size_t offset = writer->offset, bytes;
    FILE *fp = writer->fp;

    if (offset > 0) {
        // write remaining data in buffer
        if (!write_buffered_data(inBuf, offset, outBuf, outBufSize, cctx, fp))
            return -1;
    }

    // flush the Zstandard stream
    ZSTD_inBuffer emptyInput = { NULL, 0, 0 };
    ZSTD_outBuffer output = { outBuf, outBufSize, 0 };
    size_t ret;
    do {
        output.pos = 0;
        ret = ZSTD_compressStream2(cctx, &output, &emptyInput, ZSTD_e_end);
        if (ZSTD_isError(ret)) {
            errmsg("zstd compression error: %s\n", ZSTD_getErrorName(ret));
            return -2;
        }
        if (output.pos > 0) {
            if (fwrite(outBuf, output.pos, 1, fp) != 1) {
                errmsg("failed to write zstd-compressed data to file");
                return -3;
            }
        }
    } while (ret != 0);
    fflush(fp);

    writer->offset = 0;

    return offset;
}

static bool write_block(bam1_t **buf, size_t n, FILE *fp, int n_threads)
{
    zstd_cctx_t *zstd_cctx = zstd_cctx_init(ZSTD_COMPRESSION_LEVEL, n_threads);
    if (!zstd_cctx) {
        errmsg("failed to initialise zstd writer");
        return false;
    }
    ZSTD_CCtx *cctx = zstd_cctx->cctx;
    size_t inBufSize = zstd_cctx->inBufSize, outBufSize = zstd_cctx->outBufSize;
    char *inBuf = zstd_cctx->inBuf, *outBuf = zstd_cctx->outBuf;
    size_t i, bytes, offset = 0;
    bam1_t *b;
    for (i = 0; i < n; i++) {
        b = buf[i];
        bytes = sizeof(bam1_t) + b->l_data;
        if (offset + bytes > inBufSize) {
            // buffer is full
            if (!write_buffered_data(inBuf, offset, outBuf, outBufSize, cctx, fp))
                return false;
            offset = 0;
        }
        if (bytes > inBufSize) {
            // large record that cannot fit in the buffer
            if (!write_buffered_data((char *)b, bytes, outBuf, outBufSize, cctx, fp))
                return false;
        } else {
            // copy into input buffer
            memcpy(inBuf + offset, b, sizeof(bam1_t));
            memcpy(inBuf + offset + sizeof(bam1_t), b->data, b->l_data);
            offset += bytes;
        }
    }
    if (offset > 0) {
        // write remaining data in buffer
        if (!write_buffered_data(inBuf, offset, outBuf, outBufSize, cctx, fp))
            return false;
    }

    // flush the Zstandard stream
    ZSTD_inBuffer emptyInput = { NULL, 0, 0 };
    ZSTD_outBuffer output = { outBuf, outBufSize, 0 };
    size_t ret;
    do {
        output.pos = 0;
        ret = ZSTD_compressStream2(cctx, &output, &emptyInput, ZSTD_e_end);
        if (ZSTD_isError(ret)) {
            errmsg("zstd compression error: %s\n", ZSTD_getErrorName(ret));
            return false;
        }
        if (output.pos > 0) {
            if (fwrite(outBuf, output.pos, 1, fp) != 1) {
                errmsg("failed to write zstd-compressed data to file");
                return false;
            }
        }
    } while (ret != 0);
    fflush(fp);
    zstd_cctx_destroy(zstd_cctx);

    return true;
}

static zstd_writer_t *zstd_writer_init(int compressionLevel, int threads, FILE *fp)
{
    zstd_cctx_t *zstd_cctx = zstd_cctx_init(compressionLevel, threads);
    if (!zstd_cctx) {
        errmsg("failed to initialise zstd writer");
        return NULL;
    }

    zstd_writer_t *writer = (zstd_writer_t *) malloc(sizeof(zstd_writer_t));
    if (!writer) {
        errmsg("memory allocation failed");
        zstd_cctx_destroy(zstd_cctx);
        return NULL;
    }

    writer->zstd_cctx = zstd_cctx;
    writer->cctx = zstd_cctx->cctx;
    writer->inBufSize = zstd_cctx->inBufSize;
    writer->outBufSize = zstd_cctx->outBufSize;
    writer->inBuf = zstd_cctx->inBuf;
    writer->outBuf = zstd_cctx->outBuf;
    writer->offset = 0;
    writer->fp = fp;
    
    return writer;
}

static void zstd_writer_destroy(zstd_writer_t *writer)
{
    if (writer) {
        zstd_cctx_destroy(writer->zstd_cctx);
        free(writer);
    }
}

static sam_hdr_lite_t *zstd_sam_hdr_lite_read(FILE *fp)
{
    if (!fp) {
        errmsg("invalid file pointer in zstd_sam_hdr_lite_read");
        return NULL;
    }

    sam_hdr_lite_t *hdr = NULL;
    zstd_dctx_t *zstd_dctx = NULL;

    hdr = sam_hdr_lite_init();
    if (!hdr) {
        errmsg("failed to create SAM header");
        goto err;
    }
    zstd_dctx = zstd_dctx_init(0, 0);
    if (!zstd_dctx) {
        errmsg("failed to initialise ZSTD decompressor");
        goto err;
    }

    ZSTD_inBuffer *inBuf = &zstd_dctx->inBuf;
    ZSTD_outBuffer outBuf;
    char *buf = zstd_dctx->outBuf;
    size_t cap = zstd_dctx->outBufSize, used = 0;
    ssize_t rem, bytes;
    int ret = 1;
    while (ret) {
        // ensure space
        if (used == cap) {
            // no newline so far in an overlong line: grow buffer
            size_t new_cap = cap << 1;
            char *tmp = (char *) realloc(buf, new_cap + 1);
            if (!tmp) {
                errmsg("header buffer realloc failed (from %zu to %zu)", cap, new_cap);
                goto err;
            }
            buf = zstd_dctx->outBuf = tmp;
            cap = zstd_dctx->outBufSize = new_cap;
        }
        rem = cap - used; // remaining space
        
        // fill input buffer
        while (rem && ret) {
            outBuf.dst = buf + used;
            outBuf.size = rem;
            outBuf.pos = 0;
            if (inBuf->pos == inBuf->size) {
                // need to read more compressed data
                bytes = fread((void *)inBuf->src, 1, ZSTD_INBUF_SIZE, fp);

                if (bytes == 0) {
                    if (feof(fp)) {
                        ret = 0; // end of file
                        break;
                    } else if (ferror(fp)) {
                        errmsg("error reading zstd-compressed data");
                        goto err;
                    }
                }

                inBuf->pos = 0;
                inBuf->size = bytes;
            }

            if (inBuf->size > 0) {
                ret = ZSTD_decompressStream(zstd_dctx->dctx, &outBuf, inBuf);
                if (ZSTD_isError(ret)) {
                    errmsg("zstd decompression error: %s", ZSTD_getErrorName(ret));
                    goto err;
                }
            }
            used += outBuf.pos;
            rem -= outBuf.pos;
        }
        // consume lines        
        bytes = sam_hdr_lite_add_lines(hdr, buf, used);
        if (bytes < 0) {
            errmsg("failed to add SAM header lines");
            goto err;
        }
        // shift remainder
        rem = used - bytes;
        if (rem)
            memmove(buf, buf + bytes, rem);
        used = rem;
    }
    
    // rewind file pointer to the end of header
    fseek(fp, -(long)inBuf->size + (long)inBuf->pos, SEEK_CUR);

    zstd_dctx_destroy(zstd_dctx);
    return hdr;

err:
    zstd_dctx_destroy(zstd_dctx);
    sam_hdr_lite_destroy(hdr);
    return NULL;
}

static inline void fix_internal_nulls(cstr_t *cs)
{
    if (!cs || !cs->string || cs->length <= 0) return;
    char *s = cs->string;
    char *e = s + cs->length;
    *e = '@'; // fake terminator
    while (s++ < e)
        if (s[-1] == '\0')
            s[-1] = *s == '@'? '\n' : '\t';
    *e = '\0'; // real terminator
}

static int zstd_sam_hdr_lite_write(zstd_writer_t *writer, sam_hdr_lite_t *hdr)
{
    if (!writer || !writer->fp || !hdr) {
        errmsg("invalid arguments to zstd_sam_hdr_lite_write");
        return -1;
    }
    
    size_t i, n, l_text = 0;
    cstr_array_t *text;
    if (hdr->hd_text.l && /*@HD content*/
        !write_buffered_data(hdr->hd_text.s, hdr->hd_text.l,
        writer->outBuf, writer->outBufSize, writer->cctx, writer->fp)) {
        errmsg("failed to write zstd BAM header @HD line");
        return -2;
    }
    // write @SQ lines
    text = &hdr->sq_text;
    for (i = 0, n = text->n; i < n; i++) {
        fix_internal_nulls(&text->a[i]);
        if (text->a[i].length && /*@SQ content*/
            !write_buffered_data(text->a[i].string, text->a[i].length,
            writer->outBuf, writer->outBufSize, writer->cctx, writer->fp)) {
            errmsg("failed to write zstd BAM header @SQ line");
            return -2;
        }
        l_text += text->a[i].length;
    }
    // write @RG lines
    text = &hdr->rg_text;
    for (i = 0, n = text->n; i < n; i++) {
        fix_internal_nulls(&text->a[i]);
        if (text->a[i].length && /*@RG content*/
            !write_buffered_data(text->a[i].string, text->a[i].length,
            writer->outBuf, writer->outBufSize, writer->cctx, writer->fp)) {
            errmsg("failed to write zstd BAM header @RG line");
            return -2;
        }
        l_text += text->a[i].length;
    }
    // write @PG lines
    text = &hdr->pg_text;
    for (i = 0, n = text->n; i < n; i++) {
        fix_internal_nulls(&text->a[i]);
        if (text->a[i].length && /*@PG content*/
            !write_buffered_data(text->a[i].string, text->a[i].length,
            writer->outBuf, writer->outBufSize, writer->cctx, writer->fp)) {
            errmsg("failed to write zstd BAM header @PG line");
            return -2;
        }
        l_text += text->a[i].length;
    }
    // write @CO lines
    text = &hdr->co_text;
    for (i = 0, n = text->n; i < n; i++) {
        fix_internal_nulls(&text->a[i]);
        if (text->a[i].length && /*@CO content*/
            !write_buffered_data(text->a[i].string, text->a[i].length,
            writer->outBuf, writer->outBufSize, writer->cctx, writer->fp)) {
            errmsg("failed to write zstd BAM header @CO line");
            return -2;
        }
        l_text += text->a[i].length;
    }

    if ((zstd_writer_flush(writer) /*end the frame*/) != 0) {
        errmsg("failed to flush zstd BAM header");
        return -3;
    }

    infomsg("written SAM header: %d sequences; %zu characters", hdr->n_targets, l_text);

    return 0;
}

static inline int zstd_bam_write1(zstd_writer_t *writer, bam1_t *b)
{
    ZSTD_CCtx *cctx = writer->cctx;
    size_t inBufSize = writer->inBufSize, outBufSize = writer->outBufSize;
    char *inBuf = writer->inBuf, *outBuf = writer->outBuf;
    size_t offset = writer->offset, bytes;
    FILE *fp = writer->fp;
    
    bytes = sizeof(bam1_t) + b->l_data;
    if (offset + bytes > inBufSize) {
        // buffer is full
        if (!write_buffered_data(inBuf, offset, outBuf, outBufSize, cctx, fp))
            return -1;
        offset = 0;
    }
    if (bytes > inBufSize) {
        // large record that cannot fit in the buffer
        if (!write_buffered_data((char *)b, bytes, outBuf, outBufSize, cctx, fp))
            return -2;
    } else {
        // copy into input buffer
        memcpy(inBuf + offset, b, sizeof(bam1_t));
        memcpy(inBuf + offset + sizeof(bam1_t), b->data, b->l_data);
        offset += bytes;
    }

    writer->offset = offset;

    return bytes;
}
/************************ ZSTD END ************************/

/******************** File Partitioning *******************/

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8

// Heuristic checks for BAM record validity
static inline void parse_bam1_core(const uint8_t *x, bam1_core_t *c) {
    c->tid        = le_to_u32(x);
    c->pos        = le_to_i32(x+4);
    uint32_t x2   = le_to_u32(x+8);
    c->bin        = x2>>16;
    c->qual       = x2>>8&0xff;
    c->l_qname    = x2&0xff;
    c->l_extranul = (c->l_qname%4 != 0)? (4 - c->l_qname%4) : 0;
    uint32_t x3   = le_to_u32(x+12);
    c->flag       = x3>>16;
    c->n_cigar    = x3&0xffff;
    c->l_qseq     = le_to_u32(x+16);
    c->mtid       = le_to_u32(x+20);
    c->mpos       = le_to_i32(x+24);
    c->isize      = le_to_i32(x+28);
}

static inline int bam_read1_lite(samFile *fp, sam_hdr_lite_t *header, bam1_t *b)
{
    return bam_read1(fp->fp.bgzf, b);
}

static int bgzf_next_sam(BGZF *fp, uint8_t *buffer, sam_hdr_lite_t *header) {
    const int min_bytes = 32 + 4;
    int64_t i, offset, coffset, bytes;
    int32_t block_len;
    uint32_t new_l_data;
    bam1_t *b = bam_init1();
    bam1_core_t *c = &b->core;
    int ret = 0;

    while (1) {
        // Read from current block-aligned position
        offset  = bgzf_tell(fp);
        bytes   = bgzf_read(fp, buffer, BGZF_MAX_BLOCK_SIZE);
        coffset = bgzf_tell(fp);
        if (bytes == 0 && fp->fp->at_eof)
            goto found; // EOF
        // Scan through the buffer looking for valid record starts
        for (i = 0; i + min_bytes < bytes; i++) {
            block_len = le_to_u32(buffer+i);
            if (fp->is_be)
                ed_swap_4p(&block_len);
            if (block_len < 32) continue;
            parse_bam1_core(buffer+i+4, c);
            if (header) {
                if (c->tid  >= header->n_targets || c->tid  < -1 ||
                    c->mtid >= header->n_targets || c->mtid < -1) {
                    continue;
                }
            }
            new_l_data = block_len - 32 + c->l_extranul;
            if (new_l_data > INT_MAX || c->l_qseq < 0 || c->l_qname < 1)
                continue;
            if (((uint64_t) c->n_cigar << 2) + c->l_qname + c->l_extranul
                + (((uint64_t) c->l_qseq + 1) >> 1) + c->l_qseq > (uint64_t) new_l_data)
                continue;
            // try read the full record
            // so it accutally positions to the next record
            if (bgzf_seek(fp, offset + i, SEEK_SET) < 0)
                break;
            if (bam_read1(fp, b) >= -1)
                goto found;
        }

        if (bgzf_seek(fp, coffset, SEEK_SET) < 0)
            break;
    }
    ret = -1;

found:
    bam_destroy1(b);
    return ret;
}

static inline int unpackInt16(const uint8_t *buffer)
{
    return buffer[0] | buffer[1] << 8;
}

// Returns: 0 on success (BGZF header); -1 on non-BGZF GZIP header; -2 on error
static inline int check_header(const uint8_t *header)
{
    if ( header[0] != 31 || header[1] != 139 || header[2] != 8 ) return -2;
    return ((header[3] & 4) != 0
            && unpackInt16((uint8_t*)&header[10]) == 6
            && header[12] == 'B' && header[13] == 'C'
            && unpackInt16((uint8_t*)&header[14]) == 2) ? 0 : -1;
}

// Try to find the start of a BGZF block by scanning forward
// Returns 0 on success, -1 on failure
static int bgzf_resync(BGZF *fp, uint8_t *buffer) {
    int64_t i, offset, voffset, bsize, bytes = 0;
    uint32_t isize, cksum;

    offset = bgzf_tell(fp) >> 16; // current compressed position
    while (1) {
        if (fp->fp->at_eof) return 0; // reached EOF
        // read next chunk of raw compressed data
        bytes += hread(fp->fp, buffer+bytes, BGZF_MAX_BLOCK_SIZE-bytes);
        for (i = 0; i + BLOCK_HEADER_LENGTH < bytes; i++) {
            // verify bgzf header
            if (check_header(buffer+i)) continue;
            // verify bgzf block - by reading it
            voffset = (offset + i) << 16;
            if (bgzf_seek(fp, voffset, SEEK_SET) < 0)
                return -1;
            if (bgzf_read_block(fp) < 0) continue;
            if (bgzf_seek(fp, voffset, SEEK_SET) < 0)
                return -1;
            return 0;
        }
        offset += i;
        bytes -= i;
        memmove(buffer, buffer+i, bytes);
    }

    return -1;
}

static int64_t bam_next_sam(samFile *fp, int64_t offset, sam_hdr_lite_t *header)
{
    uint8_t buffer[BGZF_MAX_BLOCK_SIZE];
    BGZF *bgzf = fp->fp.bgzf;
    if (bgzf->fp->at_eof) {
        return bgzf_tell(bgzf);
    } else if (bgzf_seek(bgzf, offset, SEEK_SET) < 0) {
        errmsg("failed to seek to the file point");
        return -1;
    }
    if (bgzf->fp->at_eof) {
        return bgzf_tell(bgzf);
    } else if (bgzf_resync(bgzf, buffer) < 0) {
        errmsg("failed to resync to BGZF block");
        return -1;
    }
    if (bgzf->fp->at_eof) {
        return bgzf_tell(bgzf);
    } else if (bgzf_next_sam(bgzf, buffer, header) < 0) {
        errmsg("failed to find next BAM record");
        return -1;
    }
    return bgzf_tell(bgzf);
}

static inline int sam_read1_lite(samFile *fp, sam_hdr_lite_t *header, bam1_t *b)
{
    // parse SAM line in fp->line
    int ret;
    ret = hts_getline(fp, KS_SEP_LINE, &fp->line);
    if (ret < 0) return ret; // -1: EOF; -2: error
    ret = sam_parse1_lite(&fp->line, header, b);
    fp->line.l = 0;
    if (ret < 0) ret = -3; // indicate parse error
    return ret;
}

static int sam_resync(samFile *fp, sam_hdr_lite_t *header)
{
    if (hts_getline(fp, KS_SEP_LINE, &fp->line) < 0 && 
        !fp->fp.hfile->at_eof) {
        errmsg("failed to read next line");
        return -1;
    }
    int ret = -1;
    bam1_t *b = bam_init1();
    if (sam_read1_lite(fp, header, b) < -1) {
        errmsg("failed to read next SAM record");
        goto err;
    }
    ret = 0;
err:
    bam_destroy1(b);
    return ret;
}

static int64_t sam_next_sam(samFile *fp, int64_t offset, sam_hdr_lite_t *header)
{
    if (fp->fp.hfile->at_eof) {
        return htell(fp->fp.hfile);
    } else if (hseek(fp->fp.hfile, offset, SEEK_SET) != offset) {
        errmsg("failed to seek to the file point");
        return -1;
    }
    if (fp->fp.hfile->at_eof) {
        return htell(fp->fp.hfile);
    } else if (sam_resync(fp, header) < 0) {
        errmsg("failed to resync to next line");
        return -1;
    }
    return htell(fp->fp.hfile);
}

static sam_hdr_lite_t *bio_create_partitions(char *infile, int n_part, bio_data_t *bio, int *_part)
{
    sam_hdr_lite_t *header = NULL, *ret = NULL;
    samFile *fp = NULL;
    size_t file_size;
    struct stat st;
    int64_t boff, eoff, coff, chunk_size = 0;
    int i, part;

    if (n_part < 1) n_part = 1;
    part = n_part; // default to n_part

    if (stat(infile, &st) != 0) {
        errmsg("failed to stat input file: %s", infile);
        goto err;
    }
    file_size = st.st_size;
    
    fp = sam_open(infile, "r");
    if (!fp) {
        errmsg("failed to open input file: %s", infile);
        goto err;
    }

    header = sam_hdr_lite_read(fp);
    if (!header) {
        errmsg("failed to read header: %s", infile);
        goto err;
    }

    switch(fp->format.format) {
        case sam:
            if (fp->format.compression != no_compression) {
                errmsg("compressed SAM files are not supported");
                goto err;
            }
            boff = htell(fp->fp.hfile);
            chunk_size = (file_size - boff) / n_part;
            for (i = 0; i < n_part; i++) {
                eoff = sam_next_sam(fp, boff + chunk_size, header);
                if (eoff < 0) {
                    errmsg("failed to find partition end point");
                    goto err;
                }
                bio[i].jid = i;
                bio[i].fp = sam_open(infile, "r");
                bio[i].ftype = FT_SAM;
                bio[i].ftell = sam_ftell;
                bio[i].header = header;
                bio[i].reader = sam_read1_lite;
                bio[i].refs = NULL;
                bio[i].boff = boff;
                bio[i].eoff = eoff;
                if (boff != hseek(((samFile *)bio[i].fp)->fp.hfile, boff, SEEK_SET)) {
                    errmsg("failed to seek to the partition point");
                    goto err;
                }
                boff = eoff;
                if (fp->fp.hfile->at_eof) {
                    part = i + 1; // actual number of partitions
                    break;
                }
            }
            bio[part-1].eoff = INT64_MAX; // ensure last partition goes to the end of file
            break;
        case bam:
            boff = bgzf_tell(fp->fp.bgzf);
            assert((boff & 0xFFFF) == 0);
            file_size <<= 16;
            chunk_size = (file_size - boff) / n_part;
            for (i = 0; i < n_part; i++) {
                eoff = bam_next_sam(fp, boff + chunk_size, header);
                if (eoff < 0) {
                    errmsg("failed to find partition end point");
                    goto err;
                }
                bio[i].jid = i;
                bio[i].fp = sam_open(infile, "r");
                bio[i].ftype = FT_BAM;
                bio[i].ftell = bam_ftell;
                bio[i].header = header;
                bio[i].reader = bam_read1_lite;
                bio[i].refs = NULL;
                bio[i].boff = boff;
                bio[i].eoff = eoff;
                if (bgzf_seek(((samFile *)bio[i].fp)->fp.bgzf, boff, SEEK_SET) ||
                    boff != bgzf_tell(((samFile *)bio[i].fp)->fp.bgzf)) {
                    errmsg("failed to seek to the partition point");
                    goto err;
                }
                boff = eoff;
                if (fp->fp.bgzf->fp->at_eof) {
                    part = i + 1; // actual number of partitions
                    break;
                }
            }
            bio[part-1].eoff = INT64_MAX; // ensure last partition goes to the end of file
            break;
        case cram:
            errmsg("CRAM files are not supported");
            goto err;
        default:
            errmsg("unknown file format");
            goto err;
    }

    ret = header; header = NULL; // release header so we don't destroy it

err:
    if (_part) *_part = part;
    sam_hdr_lite_destroy(header);
    sam_close(fp);

    return ret;
}
/****************** END File Partitioning *****************/

/************** Mapping Reference Sequences ***************/

typedef struct {
    bam1_t **beg;
    bam1_t **end;
    int32_t *refmap;
} refmap_data_t;

static void *map_refseq_core(void *args) {
    refmap_data_t *data = (refmap_data_t *)args;
    bam1_t **beg = data->beg, **end = data->end, *rec;
    int32_t *refmap = data->refmap;
    while (beg < end) {
        rec = *beg++;
        rec->core.tid = refmap[rec->core.tid];
    }
    return NULL;
}

static void map_refseq(bio_data_t *bio, size_t n_bio, int n_threads, int32_t *refmap)
{
    size_t i, b, nrec;
    bam1_t **recs;

    if (n_threads <= 1) {
        for (b = 0; b < n_bio; b++) {
            recs = bio[b].recs;
            nrec = bio[b].nrec;
            for (i = 0; i < nrec; i++)
                recs[i]->core.tid = refmap[recs[i]->core.tid];
        }
        return;
    }

    pthread_t *threads = malloc(n_threads * sizeof(pthread_t));
    refmap_data_t *data = malloc(n_threads * sizeof(refmap_data_t));
    if (!threads || !data) {
        errmsg("failed to allocate memory for threads");
        exit(1);
    }
    for (b = 0; b < n_bio; b++) {
        recs = bio[b].recs;
        nrec = bio[b].nrec;

        if (!nrec) continue;

        for (i = 0; i < n_threads; i++) {
            data[i].beg = recs + (nrec * i) / n_threads;
            data[i].refmap = refmap;
            if (i > 0) data[i-1].end = data[i].beg;
        }
        data[n_threads-1].end = recs + nrec;
        for (i = 1; i < n_threads; i++)
            pthread_create(&threads[i], NULL, map_refseq_core, data + i);
        map_refseq_core(data);
        for (i = 1; i < n_threads; i++)
            pthread_join(threads[i], NULL);
    }
    free(threads);
    free(data);

    return;
}

/************ END Mapping Reference Sequences *************/

/*************** Multithreaded BAM Writing ****************/

// BGZF constants
#define BGZF_MAX_BLOCK_SIZE 0x10000
#define BGZF_BLOCK_SIZE 0xff00
#define BGZF_HEADER_SIZE 18
#define BGZF_FOOTER_SIZE 8
static const unsigned char BGZF_BLOCK_HEADER[18] = {
    0x1f, 0x8b, 0x08, 0x04, // gzip magic, CM=deflate, FLG=FEXTRA
    0x00, 0x00, 0x00, 0x00, // MTIME
    0x00,                   // XFL
    0xff,                   // OS = unknown
    0x06, 0x00,             // XLEN = 6
    0x42, 0x43, 0x02, 0x00, // BGZF subfield 'BC', SLEN=2
    0x00, 0x00              // BSIZE (patched below to 27 -> total 28 bytes)
};

static inline void tq_push(tiny_queue_t *q, int x)
{
	q->a[((q->count++) + q->front) & TINY_QUEUE_MASK] = x;
}

static inline void tq_shift(tiny_queue_t *q)
{
	if (q->count == 0) return;
	q->front++;
	q->front &= TINY_QUEUE_MASK;
	q->count--;
}

static inline int tq_front(tiny_queue_t *q)
{
	if (q->count == 0) return -1;
	return q->a[q->front];
}

static inline int tq_pop(tiny_queue_t *q)
{
    int x = tq_front(q);
    tq_shift(q);
    return x;
}

enum WWORKER_STATUS { WW_INPUT, WW_COMPR, WW_WRITE };

static void *bw_compr_worker(void *arg);
static void *bw_main_worker(void *arg);

static bam_writer_t *bam_writer_init(FILE *fp, int n_workers, size_t ubuf_cap)
{
    // sanity checks
    if (!fp) return NULL;
    if (n_workers < 1) n_workers = 1;
    if (n_workers > MAX_WWORKER_THREADS) n_workers = MAX_WWORKER_THREADS;
    if (ubuf_cap < WW_MIN_MEM) ubuf_cap = WW_MIN_MEM;
    if (ubuf_cap > WW_MAX_MEM) ubuf_cap = WW_MAX_MEM;
    // estimate buffer capacities
    size_t brec_cap = ubuf_cap / BAM_REC_SIZE;
    if (brec_cap < WW_MIN_REC) brec_cap = WW_MIN_REC;
    if (brec_cap > WW_MAX_REC) brec_cap = WW_MAX_REC;

    bam_writer_t *w = (bam_writer_t*)calloc(1, sizeof(bam_writer_t));
    if (!w) return NULL;
    w->fp = fp;
    w->n_workers = n_workers;
    w->workers = (bam_wworker_t*)calloc(n_workers, sizeof(bam_wworker_t));
    if (!w->workers) return NULL;
    w->compr_threads = (pthread_t*)calloc(n_workers, sizeof(pthread_t));
    if (!w->compr_threads) return NULL;

    w->input_q.count = w->input_q.front = 0;
    w->compr_q.count = w->compr_q.front = 0;

    pthread_mutex_init(&w->mtx, NULL);
    pthread_cond_init(&w->cond_input, NULL);
    pthread_cond_init(&w->cond_compr, NULL);
    pthread_cond_init(&w->cond_write, NULL);
    w->stop = 0;

    int i;
    for (i = 0; i < n_workers; i++) {
        bam_wworker_t *bw = w->workers + i;
        bw->writer = w;
        bw->tid = i;
        bw->status = WW_INPUT;
        tq_push(&w->input_q, i); // add to input queue
        bw->brec = (bam1_t **)malloc(brec_cap * sizeof(bam1_t *));
        if (!bw->brec) return NULL;
        bw->bcap = brec_cap;
        bw->bcnt = 0;
        bw->ubuf = (uint8_t *)malloc(ubuf_cap);
        if (!bw->ubuf) return NULL;
        bw->ucap = ubuf_cap;
        bw->ulen = 0;
        bw->cbuf = (uint8_t *)malloc(ubuf_cap);
        if (!bw->cbuf) return NULL;
        bw->ccap = ubuf_cap;
        bw->clen = 0;
        bw->bbuf = (uint8_t *)malloc(BGZF_BLOCK_SIZE);
        bw->blen = 0;
        if (!bw->bbuf) return NULL;
    }

    // set input worker to the first worker
    int t = tq_pop(&w->input_q);
    tq_push(&w->compr_q, t);
    w->input = w->workers + t;
    
    // spawn worker threads
    for (i = 0; i < n_workers; i++)
        pthread_create(&w->compr_threads[i], NULL, &bw_compr_worker, &w->workers[i]);
    pthread_create(&w->main_thread, NULL, &bw_main_worker, w);
    return w;
}

static void bam_writer_destroy(bam_writer_t *bw)
{
    if (!bw) return;
    int i;
    bam_wworker_t *worker;
    for (i = 0; i < bw->n_workers; i++) {
        worker = bw->workers + i;
        free(worker->brec);
        free(worker->cbuf);
        free(worker->bbuf);
        free(worker->ubuf);
        if (worker->args && worker->argfree)
            worker->argfree(worker->args);
    }
    free(bw->workers);
    free(bw->compr_threads);
    pthread_mutex_destroy(&bw->mtx);
    pthread_cond_destroy(&bw->cond_input);
    pthread_cond_destroy(&bw->cond_compr);
    pthread_cond_destroy(&bw->cond_write);
    if (bw->args && bw->argfree)
        bw->argfree(bw->args);
    free(bw);
}

static void *bw_main_worker(void *arg)
{
    bam_writer_t *worker = (bam_writer_t*)arg;
    for (;;) {
        int w = -1, t;
        // wait until compr/write job is available OR stop
        pthread_mutex_lock(&worker->mtx);
        while (worker->compr_q.count > 0 || !worker->stop) {
            t = tq_front(&worker->compr_q);
            if (t >= 0 && worker->workers[t].status == WW_WRITE) {
                w = t;
                break;
            }
            pthread_cond_wait(&worker->cond_write, &worker->mtx);
        }
        if (w < 0 && worker->stop) { // drained
            pthread_mutex_unlock(&worker->mtx);
            break;
        }
        pthread_mutex_unlock(&worker->mtx);
        // write compressed data
        bam_wworker_t *bw = &worker->workers[w];
        if (bw->clen &&
            fwrite(bw->cbuf, 1, bw->clen, worker->fp) != bw->clen) {
            errmsg("Error writing to output BAM file");
            exit(1);
        }
        // reset to accepting new input
        // safe to be unprotected by mutex
        bw->bcnt = 0;
        bw->ulen = 0;
        pthread_mutex_lock(&worker->mtx);
        bw->status = WW_INPUT;
        tq_shift(&worker->compr_q); // remove from compr_q
        tq_push(&worker->input_q, w); // add to input_q
        pthread_cond_signal(&worker->cond_input);
        pthread_mutex_unlock(&worker->mtx);
    }
    return NULL;
}

static void *bw_compr_worker(void *arg)
{
    bam_wworker_t *bw = (bam_wworker_t*)arg;
    bam_writer_t *worker = bw->writer;
    for (;;) {
        pthread_mutex_lock(&worker->mtx);
        while (worker->compr_q.count > 0 || !worker->stop) {
            if (bw->status == WW_COMPR)
                break;
            else
                pthread_cond_wait(&worker->cond_compr, &worker->mtx);
        }
        if (bw->status != WW_COMPR && worker->stop) { // drained
            pthread_mutex_unlock(&worker->mtx);
            break;
        }
        pthread_mutex_unlock(&worker->mtx);

        // perform compression
        worker->compress(bw);

        // signal main thread for compression completion
        pthread_mutex_lock(&worker->mtx);
        bw->status = WW_WRITE;
        pthread_cond_signal(&worker->cond_write);
        pthread_mutex_unlock(&worker->mtx);
    }
    return NULL;
}

static inline int bw_bam_write1(bam_writer_t *bw, bam1_t *b)
{
    bam_wworker_t *worker = bw->input;
    bam1_t **brec = worker->brec, *rec;
    size_t bcnt = worker->bcnt;
    size_t bcap = worker->bcap;
    size_t ulen = worker->ulen;
    size_t ucap = worker->ucap;
    uint8_t *data = worker->ubuf + worker->ulen;

    size_t bytes = (sizeof(*b) + b->l_data + 8 - 1) & ~((size_t)(8 - 1));
    
    if (bcnt >= bcap || ulen + bytes > ucap) {
        // switch to next worker
        pthread_mutex_lock(&bw->mtx);
        bw->input->status = WW_COMPR;
        pthread_cond_broadcast(&bw->cond_compr);
        while (bw->input_q.count == 0)
            pthread_cond_wait(&bw->cond_input, &bw->mtx);
        int t = tq_pop(&bw->input_q);
        tq_push(&bw->compr_q, t); // add to compr_q
        bw->input = bw->workers + t;
        pthread_mutex_unlock(&bw->mtx);
        worker = bw->input;
        ulen = worker->ulen;
        ucap = worker->ucap;
        brec = worker->brec;
        bcnt = worker->bcnt;
        data = worker->ubuf; // ulen and bcnt are both zeros here
    }

    if (ulen + bytes > ucap) {
        // try to increase buffer size for this extremely large record
        size_t ubuf_cap = ucap ? ucap : WW_MIN_MEM;
        if (ulen + bytes > ubuf_cap)
            ubuf_cap = (ulen + bytes) * 1.1 + BGZF_MAX_BLOCK_SIZE;
        if (ubuf_cap > WW_MAX_MEM)
            warnmsg("large BAM writer buffer size: %zu Mb", ubuf_cap >> 20);
        // ulen and bcnt is zero here
        worker->ubuf = (uint8_t *)realloc(worker->ubuf, ubuf_cap);
        if (!worker->ubuf) {
            errmsg("failed to allocate %zu Mb memory for BAM writer buffer", ubuf_cap >> 20);
            return -1;
        }
        // estimate buffer capacities
        size_t brec_cap = ubuf_cap / BAM_REC_SIZE;
        if (brec_cap > WW_MAX_REC)
            warnmsg("large BAM writer record count: %zu", brec_cap);
        worker->brec = (bam1_t **)realloc(worker->brec, brec_cap * sizeof(bam1_t *));
        if (!worker->brec) {
            errmsg("failed to allocate memory for %zu BAM records", brec_cap);
            return -1;
        }
        worker->ucap = ubuf_cap;
        worker->bcap = brec_cap;
        brec = worker->brec;
        bcnt = worker->bcnt;
        data = worker->ubuf;
    }
    
    rec = brec[bcnt] = (bam1_t *) data;
    *rec = *b;
    rec->data = data + sizeof(bam1_t);
    memcpy(rec->data, b->data, b->l_data);

    worker->bcnt += 1;
    worker->ulen += bytes;

    return 0;
}

static void bw_writer_join(bam_writer_t *writer)
{
    if (!writer) return;
    // wait for all input workers to finish
    pthread_mutex_lock(&writer->mtx);
    // add the last job to the compression queue
    writer->input->status = WW_COMPR;
    pthread_cond_broadcast(&writer->cond_compr);
    // there will be a cond_input signal when writing done
    while (writer->compr_q.count > 0)
        pthread_cond_wait(&writer->cond_input, &writer->mtx);
    writer->stop = 1;
    pthread_cond_broadcast(&writer->cond_compr);
    pthread_cond_broadcast(&writer->cond_write);
    pthread_mutex_unlock(&writer->mtx);

    // join all threads
    int i;
    for (i = 0; i < writer->n_workers; i++)
        pthread_join(writer->compr_threads[i], NULL);
    pthread_join(writer->main_thread, NULL);
}

static void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host)
{
    uint32_t *cigar = (uint32_t*)(data + c->l_qname);
    uint32_t i;
    for (i = 0; i < c->n_cigar; ++i) ed_swap_4p(&cigar[i]);
}

static inline int bw_flush(bam_wworker_t *bw, FILE *fp)
{
    if (bw->blen == 0) return 0;

    // Ensure space for at least one full BGZF block
    if (bw->clen + BGZF_MAX_BLOCK_SIZE > bw->ccap) {
        size_t ccap = bw->ccap ? bw->ccap : BGZF_MAX_BLOCK_SIZE;
        while (bw->clen + BGZF_MAX_BLOCK_SIZE > ccap)
            ccap = ccap * 1.1 + BGZF_MAX_BLOCK_SIZE;
        uint8_t *cbuf = (uint8_t*)realloc(bw->cbuf, ccap);
        if (!cbuf) return -1;
        bw->cbuf = cbuf;
        bw->ccap = ccap;
    }

    if (bw->blen > BGZF_BLOCK_SIZE) return -1; // caller's responsibility

    uint8_t *dst = bw->cbuf + bw->clen;
    memcpy(dst, BGZF_BLOCK_HEADER, BGZF_HEADER_SIZE); // header template

    // Set up raw DEFLATE (no zlib/gzip wrapper, windowBits = -15)
    z_stream zs; memset(&zs, 0, sizeof(zs));
    if (deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
        return -1;
    zs.next_in   = bw->bbuf;
    zs.avail_in  = (uInt)bw->blen;
    zs.next_out  = dst + BGZF_HEADER_SIZE;
    zs.avail_out = (uInt)(BGZF_MAX_BLOCK_SIZE - BGZF_HEADER_SIZE - BGZF_FOOTER_SIZE);

    int zret = deflate(&zs, Z_FINISH);
    if (zret != Z_STREAM_END) {
        deflateEnd(&zs);
        return -1;
    }
    size_t comp_size = zs.total_out;
    deflateEnd(&zs);

    size_t total = BGZF_HEADER_SIZE + comp_size + BGZF_FOOTER_SIZE;
    if (total > BGZF_MAX_BLOCK_SIZE) return -1; // should not happen

    // Patch BSIZE (total - 1) at bytes 16..17
    uint16_t bsize = (uint16_t)(total - 1);
    if (IS_BIG_ENDIAN) {
        uint16_t tmp = bsize;
        ed_swap_2p(&tmp);
        memcpy(dst + 16, &tmp, 2);
    } else {
        dst[16] = (uint8_t)(bsize & 0xff);
        dst[17] = (uint8_t)(bsize >> 8);
    }

    // Footer: CRC32 and ISIZE over uncompressed payload
    uint32_t crc = crc32(crc32(0L, Z_NULL, 0), bw->bbuf, (uInt)bw->blen);
    uint32_t isize = (uint32_t)bw->blen; // use actual uncompressed size
    if (IS_BIG_ENDIAN) { // file format is little-endian; swap if host big-endian
        uint32_t tcrc = crc, tisize = isize;
        ed_swap_4p(&tcrc); ed_swap_4p(&tisize);
        memcpy(dst + BGZF_HEADER_SIZE + comp_size,     &tcrc,   4);
        memcpy(dst + BGZF_HEADER_SIZE + comp_size + 4, &tisize, 4);
    } else {
        memcpy(dst + BGZF_HEADER_SIZE + comp_size,     &crc,   4);
        memcpy(dst + BGZF_HEADER_SIZE + comp_size + 4, &isize, 4);
    }

    bw->blen = 0; // clear the block buffer

    if (fp) {
        // write newly flushed block out
        // bw->clen won't change here
        if (fwrite(bw->cbuf + bw->clen, 1, total, fp) != total)
            return -1;
    } else bw->clen += total;

    return 0;
}

static inline int bw_flush_try(bam_wworker_t *bw, ssize_t size)
{
    if (bw->blen + size > BGZF_BLOCK_SIZE)
        return bw_flush(bw, NULL);
    return 0;
}

static inline ssize_t bw_write_small(bam_wworker_t *bw, const void *data, size_t size, FILE *fp)
{
    const uint8_t *src = (const uint8_t *) data;
    size_t copied = 0;
    while (copied < size) {
        size_t avail = BGZF_BLOCK_SIZE - bw->blen;
        if (avail == 0) {
            if (bw_flush(bw, fp) < 0) return -1;
            avail = BGZF_BLOCK_SIZE - bw->blen;
        }
        size_t chunk = size - copied;
        if (chunk > avail) chunk = avail;
        memcpy(bw->bbuf + bw->blen, src + copied, chunk);
        bw->blen += chunk;
        copied += chunk;
    }
    return (ssize_t) size;
}

static inline int bw_bam_write(bam1_t *b, bam_wworker_t *bw) 
{
    const bam1_core_t *c = &b->core;
    uint32_t x[8], block_len = b->l_data - c->l_extranul + 32, y;
    int i, ok;
    if (c->l_qname - c->l_extranul > 255) {
        errmsg("QNAME \"%s\" is longer than 254 characters", bam_get_qname(b));
        return -1;
    }
    if (c->n_cigar > 0xffff) {
        errmsg("Long CIGAR string is not supported yet: %u", c->n_cigar);
        return -1;
    }
    if (c->pos > INT_MAX ||
        c->mpos > INT_MAX ||
        c->isize < INT_MIN || c->isize > INT_MAX) {
        errmsg("Positional data is too large for BAM format");
        return -1;
    }
    x[0] = c->tid;
    x[1] = c->pos;
    x[2] = (uint32_t)c->bin<<16 | c->qual<<8 | (c->l_qname - c->l_extranul);
    if (c->n_cigar > 0xffff) x[3] = (uint32_t)c->flag << 16 | 2;
    else x[3] = (uint32_t)c->flag << 16 | (c->n_cigar & 0xffff);
    x[4] = c->l_qseq;
    x[5] = c->mtid;
    x[6] = c->mpos;
    x[7] = c->isize;
    ok = (bw_flush_try(bw, 4 + block_len) >= 0);
    if (IS_BIG_ENDIAN) {
        for (i = 0; i < 8; ++i) ed_swap_4p(x + i);
        y = block_len;
        if (ok) ok = (bw_write_small(bw, ed_swap_4p(&y), 4, NULL) >= 0);
        swap_data(c, b->l_data, b->data, 1);
    } else {
        if (ok) ok = (bw_write_small(bw, &block_len, 4, NULL) >= 0);
    }
    if (ok) ok = (bw_write_small(bw, x, 32, NULL) >= 0);
    if (ok) ok = (bw_write_small(bw, b->data, c->l_qname - c->l_extranul, NULL) >= 0);
    if (ok) ok = (bw_write_small(bw, b->data + c->l_qname, b->l_data - c->l_qname, NULL) >= 0);
    if (IS_BIG_ENDIAN) swap_data(c, b->l_data, b->data, 0);
    return ok? 4 + block_len : -1;
}

static void bw_bam_compress(bam_wworker_t *bw)
{
    size_t n = bw->bcnt;
    if (n == 0) return;
    
    bw->blen = 0;
    bw->clen = 0;

    size_t i;
    for (i = 0; i < n; i++) {
        if (bw_bam_write(bw->brec[i], bw) < 0) {
            errmsg("Error writing BAM record");
            exit(1);
        }
    }
    if (bw_flush(bw, NULL) < 0) {
        errmsg("Error flushing BGZF block");
        exit(1);
    }
}

typedef struct {
    kstring_t *str;        // string buffer for formatting SAM lines
    sam_hdr_lite_t *hdr;   // SAM header
} sam_compArgs_t;

void *sam_compArgs_make(sam_hdr_lite_t *hdr, size_t buf)
{
    sam_compArgs_t *args = NULL;
    kstring_t *str = NULL;
    
    args = (sam_compArgs_t *)calloc(1, sizeof(sam_compArgs_t));
    if (!args) goto fail;
    
    str = (kstring_t *)calloc(1, sizeof(kstring_t));
    if (!str) return NULL;
    str->s = (char *)malloc(buf);
    if (str->s) str->m = buf;

    args->str = str;
    args->hdr = hdr;

    return (void *)args;

fail:
    ks_free(str);
    free(str);
    free(args);
    return NULL;
}

void sam_compArgs_free(void *args)
{
    if (!args) return;
    sam_compArgs_t *cargs = (sam_compArgs_t *)args;
    ks_free(cargs->str);
    free(cargs->str);
    free(cargs);
}

static void bw_sam_compress(bam_wworker_t *bw)
{
    size_t n = bw->bcnt;
    if (n == 0) return;

    sam_compArgs_t *args = (sam_compArgs_t *)bw->args;
    sam_hdr_lite_t *hdr = args->hdr;
    kstring_t *str = args->str;
    bw->blen = 0;
    bw->clen = 0;

    size_t i;
    for (i = 0; i < n; i++) {
        if (sam_format1_lite(hdr, bw->brec[i], str) < 0 ||
            kputc_('\n', str) < 0 ||
            bw_flush_try(bw, str->l) < 0 || 
            bw_write_small(bw, str->s, str->l, NULL) < 0) {
            errmsg("Error writing SAM record");
            exit(1);
        }
    }
    if (bw_flush(bw, NULL) < 0) {
        errmsg("Error flushing BGZF block");
        exit(1);
    }
}

static size_t bw_hdr_lite_write_text(bam_wworker_t *bw, sam_hdr_lite_t *hdr, FILE *fp, int more)
{
    int32_t i, j, n;
    cstr_array_t *cstrs;
    size_t l_text = 0;

    // write @SQ lines
    cstrs = &hdr->sq_text;
    for (i = 0, n = cstrs->n; i < n; i++) {
        fix_internal_nulls(&cstrs->a[i]);
        if (bw_write_small(bw, cstrs->a[i].string, cstrs->a[i].length, fp) < 0)
            return -1;
        l_text += cstrs->a[i].length;
    }
    // write @RG lines
    cstrs = &hdr->rg_text;
    for (i = 0, n = cstrs->n; i < n; i++) {
        fix_internal_nulls(&cstrs->a[i]);
        if (bw_write_small(bw, cstrs->a[i].string, cstrs->a[i].length, fp) < 0)
            return -2;
        l_text += cstrs->a[i].length;
    }
    if (more < 1) return l_text; // no more lines to write
    // write @PG lines
    cstrs = &hdr->pg_text;
    for (i = 0, n = cstrs->n; i < n; i++) {
        fix_internal_nulls(&cstrs->a[i]);
        if (bw_write_small(bw, cstrs->a[i].string, cstrs->a[i].length, fp) < 0)
            return -3;
        l_text += cstrs->a[i].length;
    }
    if (more < 2) return l_text; // no more lines to write
    // write @CO lines
    cstrs = &hdr->co_text;
    for (i = 0, n = cstrs->n; i < n; i++) {
        fix_internal_nulls(&cstrs->a[i]);
        if (bw_write_small(bw, cstrs->a[i].string, cstrs->a[i].length, fp) < 0)
            return -4;
        l_text += cstrs->a[i].length;
    }

    return l_text;
}

static int bw_bam_hdr_lite_write(FILE *fp, sam_hdr_lite_t *hdr)
{
    if (!fp || !hdr) {
        errmsg("invalid arguments to bw_bam_hdr_lite_write");
        return -1;
    }

    int32_t i, x, ret = -1;
    uint32_t l_text = 0, l;
    bam_wworker_t _bw, *bw = &_bw;
    int more = 0;

    memset(bw, 0, sizeof(*bw));

    bw->bbuf = (uint8_t*)malloc(BGZF_BLOCK_SIZE);
    if (!bw->bbuf) goto err;
    bw->blen = 0;
    bw->cbuf = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    if (!bw->cbuf) goto err;
    bw->ccap = BGZF_MAX_BLOCK_SIZE;
    bw->clen = 0;

    // @HD line
    l_text += hdr->hd_text.l;
    for (i = 0; i < hdr->sq_text.n; i++)
        l_text += hdr->sq_text.a[i].length;
    for (i = 0; i < hdr->rg_text.n; i++)
        l_text += hdr->rg_text.a[i].length;
    if (l_text > BAM_MAX_HEADER_SIZE) {
        errmsg("BAM header text too long: %u characters > %u maximum", l_text, BAM_MAX_HEADER_SIZE);
        goto err;
    }
    // check if we can include @PG and @CO lines
    l = 0;
    for (i = 0; i < hdr->pg_text.n; i++)
        l += hdr->pg_text.a[i].length;
    if (l_text + l <= BAM_MAX_HEADER_SIZE) {
        l_text += l;
        more += 1;
        l = 0;
        for (i = 0; i < hdr->co_text.n; i++)
            l += hdr->co_text.a[i].length;
        if (l_text + l <= BAM_MAX_HEADER_SIZE) {
            l_text += l;
            more += 1;
        } else warnmsg("skipping @CO lines due to header size limit");
    } else warnmsg("skipping @PG lines due to header size limit");

    // write "BAM1"
    if (bw_write_small(bw, "BAM\1", 4, fp) < 0) goto err;
    // write plain text and the number of reference sequences
    if (IS_BIG_ENDIAN) {
        x = ed_swap_4(l_text);
        if (bw_write_small(bw, &x, 4, fp) < 0) goto err;
        if (bw_write_small(bw, hdr->hd_text.s, hdr->hd_text.l, fp) < 0) goto err;
        if (bw_hdr_lite_write_text(bw, hdr, fp, more) < 0) goto err;
        x = ed_swap_4(hdr->n_targets);
        if (bw_write_small(bw, &x, 4, fp) < 0) goto err;
    } else {
        if (bw_write_small(bw, &l_text, 4, fp) < 0) goto err;
        if (bw_write_small(bw, hdr->hd_text.s, hdr->hd_text.l, fp) < 0) goto err;
        if (bw_hdr_lite_write_text(bw, hdr, fp, more) < 0) goto err;
        if (bw_write_small(bw, &hdr->n_targets, 4, fp) < 0) goto err;
    }

    // write sequence names and lengths
    // @SQ lines are formalised already
    cstr_array_t *cstrs = &hdr->sq_text;
    char *line, *last, *endl, *tbeg, *tend, *name, *tlen;
    uint32_t t, name_len, target_len;
    for (i = 0; i < cstrs->n; i++) {
        line = cstrs->a[i].string;
        last = line + cstrs->a[i].length;
        *last = '\0'; // terminate cstring
        while (line < last) {
            endl = bs_strchrnul(line, '\n');
            *endl = '\0'; // terminate line
            name = tlen = NULL;
            tbeg = line;
            t = 0;
            while (tbeg < endl && t < 2) {
                tend = bs_strchrnul(tbeg, '\t');
                if (tbeg[0] == 'S' && tbeg[1] == 'N') {
                    name = tbeg + 3;
                    name_len = tend - name;
                    t++;
                } else if (tbeg[0] == 'L' && tbeg[1] == 'N') {
                    tlen = tbeg + 3;
                    t++;
                }
                tbeg = tend + 1;
            }
            if (!name || !tlen) {
                errmsg("malformed @SQ line (missing SN or LN): \"%s\"", line);
                goto err;
            }
            target_len = (uint32_t)strtoul(tlen, NULL, 10);
            // make sure name is null-terminated
            name[name_len] = '\0';
            name_len += 1; // include null character
            if (IS_BIG_ENDIAN) {
                x = ed_swap_4(name_len);
                if (bw_write_small(bw, &x, 4, fp) < 0) goto err;
                if (bw_write_small(bw, name, name_len, fp) < 0) goto err;
                x = ed_swap_4(target_len);
                if (bw_write_small(bw, &x, 4, fp) < 0) goto err;
            } else {
                if (bw_write_small(bw, &name_len, 4, fp) < 0) goto err;
                if (bw_write_small(bw, name, name_len, fp) < 0) goto err;
                if (bw_write_small(bw, &target_len, 4, fp) < 0) goto err;
            }
            line = endl + 1;
        }
    }
    if (bw->blen && bw_flush(bw, fp) < 0) goto err;

    infomsg("written SAM header: %d sequences; %zu characters", hdr->n_targets, l_text);
    
    ret = 0;

err:
    free(bw->bbuf);
    free(bw->cbuf);
    return ret;
}

static int bw_sam_hdr_lite_write(FILE *fp, sam_hdr_lite_t *hdr)
{
    if (!fp || !hdr) {
        errmsg("invalid arguments to bw_sam_hdr_lite_write");
        return -1;
    }

    int32_t i, x, ret = -1;
    bam_wworker_t _bw, *bw = &_bw;
    size_t l_text = 0;
    int more = 2; // write all lines

    memset(bw, 0, sizeof(*bw));

    bw->bbuf = (uint8_t*)malloc(BGZF_BLOCK_SIZE);
    if (!bw->bbuf) goto err;
    bw->blen = 0;
    bw->cbuf = (uint8_t*)malloc(BGZF_MAX_BLOCK_SIZE);
    if (!bw->cbuf) goto err;
    bw->ccap = BGZF_MAX_BLOCK_SIZE;
    bw->clen = 0;

    // @HD line
    if (bw_write_small(bw, hdr->hd_text.s, hdr->hd_text.l, fp) < 0) goto err;
    if ((l_text = bw_hdr_lite_write_text(bw, hdr, fp, more)) < 0) goto err;
    if (bw->blen && bw_flush(bw, fp) < 0) goto err;

    l_text += hdr->hd_text.l;

    infomsg("written SAM header: %d sequences; %zu characters", hdr->n_targets, l_text);

    ret = 0;

err:
    free(bw->bbuf);
    free(bw->cbuf);
    return ret;
}

static size_t sam_hdr_lite_min_text_len(sam_hdr_lite_t *hdr)
{
    if (!hdr) return 0;

    // minimum is @HD + @SQ + @RG lines
    // @HD line
    int i;
    size_t l_text = hdr->hd_text.l;
    for (i = 0; i < hdr->sq_text.n; i++)
        l_text += hdr->sq_text.a[i].length;
    for (i = 0; i < hdr->rg_text.n; i++)
        l_text += hdr->rg_text.a[i].length;
    return l_text;
}

// Write standard 28-byte BGZF EOF marker block by constructing its fields.
static int bw_eof_write(FILE *fp) 
{
    unsigned char hdr[18];
    memcpy(hdr, BGZF_BLOCK_HEADER, sizeof(BGZF_BLOCK_HEADER));
    unsigned char deflate_empty[2] = { 0x03, 0x00 }; // minimal empty deflate block
    unsigned char footer[8] = { 0,0,0,0, 0,0,0,0 };  // CRC32(0) + ISIZE(0)
    uint16_t bsize = 28 - 1;
    hdr[16] = (unsigned char)(bsize & 0xff);
    hdr[17] = (unsigned char)(bsize >> 8);
    // Resulting 28 bytes are the canonical BGZF EOF marker.
    if (fwrite(hdr, 1, sizeof(hdr), fp) != sizeof(hdr)) return -1;
    if (fwrite(deflate_empty, 1, sizeof(deflate_empty), fp) != sizeof(deflate_empty)) return -1;
    if (fwrite(footer, 1, sizeof(footer), fp) != sizeof(footer)) return -1;
    return 0;
}

static int bam_writer_close(bam_writer_t *bw)
{
    if (!bw) return -1;
    // wait for all workers to finish
    bw_writer_join(bw);
    // write EOF block
    if (bw_eof_write(bw->fp) != 0) {
        errmsg("Failed to write BGZF EOF block");
        return -1;
    }
    return 0;
}

/************* END Multithreaded BAM Writing **************/


/******************* SAM Header Lite *********************/
#include "hts_os.h"

sam_hdr_lite_t *sam_hdr_lite_init()
{
    sam_hdr_lite_t *hdr = (sam_hdr_lite_t *) calloc(1, sizeof(sam_hdr_lite_t));
    if (!hdr) return NULL;
    return hdr;
}

void sam_hdr_lite_destroy(sam_hdr_lite_t *hdr)
{
    if (!hdr) return;
    ks_free(&hdr->hd_text);
    cstr_array_destroy(&hdr->sq_text);
    cstr_array_destroy(&hdr->rg_text);
    cstr_array_destroy(&hdr->pg_text);
    cstr_array_destroy(&hdr->co_text);
    kh_destroy(c2i,  hdr->sq_map);
    kh_destroy(cset, hdr->rg_ids);
    kh_destroy(cset, hdr->pg_ids);
    free(hdr->target_name);
    free(hdr->target_len);
    free(hdr);
}

static sam_hdr_lite_t *sam_hdr_lite_builder_init() 
{
    sam_hdr_lite_t *builder = sam_hdr_lite_init();
    if (!builder) return NULL;
    
    builder->sq_map = kh_init(c2i);
    builder->rg_ids = kh_init(cset);
    builder->pg_ids = kh_init(cset);
    if (!builder->sq_map || 
        !builder->rg_ids || 
        !builder->pg_ids) {
        sam_hdr_lite_destroy(builder);
        return NULL;
    }
    hts_srand48((long)time(NULL));
    
    return builder;
}


static inline char *get_eol(char *ptr, char *eptr)
{
    char *p = ptr + 1, prev = *ptr;
    while (p < eptr) {
        if (*p == '@' && (prev == '\n' || prev == '\0'))
            return p - 1; // EOL found
        prev = *p++;
    }
    return eptr - 1; // no further header line start found
}

static inline char *get_tag(char *ptr, char *eptr, const char *tag)
{
    char *p, *q;
    eptr -= 4; // ensure space for tag+':'+delimiter
    for (p = ptr; p < eptr; p++) {
        if ((p[-1] == '\t' || p[-1] == '\0') &&
            p[0] == tag[0] && p[1] == tag[1] && p[2] == ':') {
            q = p += 3;
            while (q < eptr && *q && *q != '\t' && *q != '\n')
                q++;
            *q = '\0'; // terminate at delimiter
            return p;
        }
    }
    return NULL;
}

static const char *HD_VN_NUMBER = "1.6"; // default version number
static const char *HD_SO_STRING = "queryname"; // default sort order
static const char *HD_LINE_FORMAT = "@HD\tVN:%s\tSO:%s\n";

static int sam_hdr_lite_build_HD(sam_hdr_lite_t *hdr, const char *_vn, const char *_so)
{
    char *vn = NULL, *so = NULL;
    bool new_vn = false;
    if (!_vn || !*_vn) {
        if (hdr->hd_text.l > 3) {
            // try extracting version from existing @HD line
            char *ptr = hdr->hd_text.s, *eptr = ptr + hdr->hd_text.l;
            char *vtag = get_tag(ptr+4, eptr, "VN");
            if (vtag) {
                vn = strdup(vtag);
                new_vn = true;
            }
        }
        if (!vn) vn = (char *)HD_VN_NUMBER;
    } else vn = (char *)_vn;
    if (!_so || !*_so) 
        so = (char *)HD_SO_STRING;
    else so = (char *)_so;
    hdr->hd_text.l = 0;
    ksprintf(&hdr->hd_text, HD_LINE_FORMAT, vn, so);

    if (new_vn) free(vn);
    return 0;
}

static bool sam_hdr_lite_build_targets_tbl(sam_hdr_lite_t *hdr) 
{
    ssize_t i, n_targets;
    int ret;
    khiter_t k;
    char *line, *endl, *last, *tag, **target_name;

    if (hdr->target_name) { 
        free(hdr->target_name); 
        hdr->target_name = NULL;
    }

    hdr->target_name = (char **)malloc(hdr->n_targets * sizeof(char *));
    if (!hdr->target_name) goto err;
    target_name = hdr->target_name;

    // process @SQ lines
    n_targets = 0;
    for (i = 0; i < hdr->sq_text.n; i++) {
        line = hdr->sq_text.a[i].string;
        last = line + hdr->sq_text.a[i].length;
        while (line < last) {
            // get the end of line
            endl = get_eol(line, last);
            // this is built inplace
            tag = get_tag(line+4, endl, "SN");
            if (!tag) goto err;
            // put this reference in map
            target_name[n_targets++] = tag;
            line = endl + 1; // next line
        }
    }

    if (n_targets != hdr->n_targets) goto err;
    return true;

err:
    errmsg("failed to build sequence name to ID map");
    return false;
}

static bool sam_hdr_lite_build_targets_map(sam_hdr_lite_t *hdr) 
{
    ssize_t i, n_targets;
    int ret;
    khiter_t k;
    char *line, *endl, *last, *tag;

    if (hdr->sq_map) 
        kh_clear(c2i, hdr->sq_map);
    else 
        hdr->sq_map = kh_init(c2i);
    khash_t(c2i) *sq_map = hdr->sq_map;

    if (!sq_map) {
        errmsg("failed to create sequence name to ID map");
        return false;
    }

    // process @SQ lines
    n_targets = 0;
    for (i = 0; i < hdr->sq_text.n; i++) {
        line = hdr->sq_text.a[i].string;
        last = line + hdr->sq_text.a[i].length;
        while (line < last) {
            // get the end of line
            endl = get_eol(line, last);
            // this is built inplace
            tag = get_tag(line+4, endl, "SN");
            if (!tag) goto err;
            // put this reference in map
            k = kh_put(c2i, sq_map, tag, &ret);
            if (ret < 0) goto err;
            kh_val(sq_map, k) = n_targets++;
            line = endl + 1; // next line
        }
    }

    if (n_targets != hdr->n_targets) goto err;
    return true;

err:
    errmsg("failed to build sequence name to ID map");
    return false;
}

enum SAM_HDR_TAG { 
    SAM_HDR_HD = 0x4844, 
    SAM_HDR_SQ = 0x5351,
    SAM_HDR_RG = 0x5247,
    SAM_HDR_PG = 0x5047,
    SAM_HDR_CO = 0x434F
};

static inline ssize_t sam_hdr_lite_add_line(sam_hdr_lite_t *hdr, const char *line, size_t len)
{
    if (len < 4 || 
        line[0] != '@' || 
        line[len-1] != '\n')
        return -1;
    uint16_t htag = (uint16_t)line[1]<<8 | (uint16_t)line[2];
    switch (htag) {
        case SAM_HDR_HD: // @HD
            if (kputsn(line, len, &hdr->hd_text) < 0)
                goto err;
            break;
        case SAM_HDR_SQ: // @SQ
            if (cstr_array_putsn(line, len, &hdr->sq_text) < 0)
                goto err;
            hdr->n_targets++; // count targets
            break;
        case SAM_HDR_RG: // @RG
            if (cstr_array_putsn(line, len, &hdr->rg_text) < 0)
                goto err;
            break;
        case SAM_HDR_PG: // @PG
            if (cstr_array_putsn(line, len, &hdr->pg_text) < 0)
                goto err;
            break;
        case SAM_HDR_CO: // @CO
            if (cstr_array_putsn(line, len, &hdr->co_text) < 0)
                goto err;
            break;
        default:
            errmsg("unrecognized SAM header line type: \"%.*s\"", len, line);
            return -1;
    }
    return len;

err:
    errmsg("failed to add SAM header line");
    return -1;
}

static inline ssize_t sam_hdr_lite_add_lines(sam_hdr_lite_t *hdr, const char *lines, size_t len)
{
    char *line, *endl, *last;
    ssize_t llen;
    line = (char *) lines;
    last = line + len;
    while (line < last) {
        endl = line;
        while (++endl < last && *endl != '\n');
        if (endl == last || *endl != '\n')
            return (ssize_t) (line - lines);
        llen = (size_t) (endl - line) + 1; // +1 for '\n'
        if (sam_hdr_lite_add_line(hdr, line, llen) < 0) {
            infomsg("invalid SAM header line: \"%.*s\"", (int)llen, line);
            goto err;
        }
        line = endl + 1;
    }

    return (ssize_t) (line - lines);

err:
    errmsg("failed to add SAM header lines");
    return -1;
}

sam_hdr_lite_t *sam_hdr_lite_create(samFile *fp)
{
    if (!fp) return NULL;
    
    sam_hdr_lite_t *hdr = sam_hdr_lite_init();
    if (!hdr) return NULL;

    int ret;
    int next_c = '@';
    while (next_c == '@' && (ret = hts_getline(fp, KS_SEP_LINE, &fp->line)) >= 0) {
        if (fp->line.s[0] != '@')
            break;

        kputc('\n', &fp->line); // hts_getline does not retain newline

        if (sam_hdr_lite_add_line(hdr, fp->line.s, fp->line.l) < 0) {
            errmsg("Invalid header line: %s", fp->line.s);
            goto error;
        }

        if (fp->is_bgzf) {
            next_c = bgzf_peek(fp->fp.bgzf);
        } else {
            unsigned char nc;
            ssize_t pret = hpeek(fp->fp.hfile, &nc, 1);
            next_c = pret > 0 ? nc : pret - 1;
        }
        if (next_c < -1)
            goto error;
    }
    if (next_c != '@')
        fp->line.l = 0;
    
    // build target name2id map
    if (!sam_hdr_lite_build_targets_map(hdr)) {
        errmsg("Failed to build target name to ID map");
        goto error;
    }

    return hdr;

error:
    sam_hdr_lite_destroy(hdr);
    return NULL;
}

sam_hdr_lite_t *bam_hdr_lite_read(BGZF *fp)
{
    sam_hdr_lite_t *h;
    uint8_t *buf = NULL;
    int magic_len, has_EOF;
    int32_t i, name_len;
    size_t bufsize = CSTR_SIZE << 2;
    ssize_t bytes, chunk, remain, consume, l_text;
    // check EOF
    has_EOF = bgzf_check_EOF(fp);
    if (has_EOF < 0) {
        perror("[W::bam_hdr_read] bgzf_check_EOF");
    } else if (has_EOF == 0) {
        hts_log_warning("EOF marker is absent. The input is probably truncated");
    }
    // allocate buffer
    buf = (uint8_t *)malloc(bufsize + 1);
    if (!buf) goto nomem;
    // read "BAM1"
    magic_len = bgzf_read(fp, buf, 4);
    if (magic_len != 4 || memcmp(buf, "BAM\1", 4)) {
        errmsg("Invalid BAM binary header");
        return 0;
    }
    h = sam_hdr_lite_init();
    if (!h) goto nomem;
    // read plain text and parse lines
    bytes = bgzf_read(fp, buf, 4);
    if (bytes != 4) goto read_err;
    l_text = le_to_u32(buf);

    remain = 0;
    while (l_text > 0) {
        chunk = bufsize - remain;
        if (chunk == 0) {
            // a very long header line
            // expand buffer
            bufsize <<= 1;
            uint8_t *tmp = (uint8_t *) realloc(buf, bufsize + 1);
            if (!tmp) goto nomem;
            buf = tmp;
            chunk = bufsize - remain;
        }
        chunk = l_text > chunk ? chunk : l_text;
        bytes = bgzf_read(fp, buf + remain, chunk);
        if (bytes != chunk) goto read_err;
        consume = sam_hdr_lite_add_lines(h, (char *)buf, remain + bytes);
        if (consume < 0) goto parse_err;
        remain += bytes - consume;
        // move remaining data to front
        if (remain > 0)
            memmove(buf, buf + consume, remain);
        l_text -= bytes;
    }
    if (remain > 0) {
        // process remaining data
        consume = sam_hdr_lite_add_lines(h, (char *)buf, remain);
        if (consume < remain) goto parse_err;
    }

    bytes = bgzf_read(fp, &h->n_targets, 4);
    if (bytes != 4) goto read_err;
    if (fp->is_be) ed_swap_4p(&h->n_targets);

    if (h->n_targets < 0) goto invalid;

    // read reference sequence names and lengths
    // no data stored
    for (i = 0; i != h->n_targets; ++i) {
        // name length
        bytes = bgzf_read(fp, &name_len, 4);
        if (bytes != 4) goto read_err;
        if (fp->is_be) ed_swap_4p(&name_len);
        if (name_len <= 0) goto invalid;

        // target name
        while (name_len > 0) {
            chunk = name_len > bufsize ? bufsize : name_len;
            bytes = bgzf_read(fp, buf, chunk);
            if (bytes != chunk) goto read_err;
            name_len -= bytes;
        }

        // target length
        bytes = bgzf_read(fp, buf, 4);
        if (bytes != 4) goto read_err;
    }
    free(buf);
    return h;

 nomem:
    errmsg("Out of memory");
    goto clean;

 read_err:
    if (bytes < 0) {
        errmsg("Error reading BGZF stream");
    } else {
        errmsg("Truncated BAM header");
    }
    goto clean;

parse_err:
    errmsg("Failed to parse BAM header");
    goto clean;

 invalid:
    errmsg("Invalid BAM binary header");

 clean:
    free(buf);
    if (h != NULL)
        sam_hdr_lite_destroy(h);
    return NULL;
}

#ifndef EFTYPE
#define EFTYPE ENOEXEC
#endif
#ifndef EOVERFLOW
#define EOVERFLOW ERANGE
#endif

static sam_hdr_lite_t *sam_hdr_lite_read(samFile *fp)
{
    if (!fp) {
        errno = EINVAL;
        return NULL;
    }

    switch (fp->format.format) {
    case bam:
        return bam_hdr_lite_read(fp->fp.bgzf);

    case sam:
        return sam_hdr_lite_create(fp);

    case cram:

    case fastq_format:
    case fasta_format:
        errno = ENOTSUP;
        return NULL;

    case empty_format:
        errno = EPIPE;
        return NULL;

    default:
        errno = EFTYPE;
        return NULL;
    }
}

// unique @PG-ID generator - copied from htslib/bam_sort.c
static int gen_unique_id(char *prefix, khash_t(cset) *existing_ids,
                         bool always_add_suffix, kstring_t *dest) {
    khiter_t iter;

    if (!always_add_suffix) {
        // Try prefix on its own first
        iter = kh_get(cset, existing_ids, prefix);
        if (iter == kh_end(existing_ids)) { // prefix isn't used yet
            dest->l = 0;
            if (kputs(prefix, dest) == EOF) return -1;
            return 0;
        }
    }

    do {
        dest->l = 0;
        ksprintf(dest, "%s-%08lX", prefix, lrand48());
        iter = kh_get(cset, existing_ids, ks_str(dest));
    } while (iter != kh_end(existing_ids));

    return 0;
}

static void cstr_free(cstr_t *s) {
    if (s && s->string) {
        free(s->string);
        s->string = NULL;
        s->length = 0;
    }
}

static cstr_t *cstr_array_expand(cstr_array_t *arr, size_t size)
{
    if (!arr) return NULL;
    if (arr->n == arr->m) {
        size_t m = arr->m? arr->m << 1 : 16;
        void *a = realloc(arr->a, m * sizeof(cstr_t));
        if (!a) return NULL;
        arr->a = (cstr_t *) a;
        arr->m = m;
    }
    size_t cap = (size > CSTR_SIZE) ? size : CSTR_SIZE;
    char *mem = (char *)malloc(cap + 1); // +1 for '\0'
    if (!mem) return NULL;
    arr->a[arr->n].string = mem;
    arr->a[arr->n].length = 0;
    arr->a[arr->n].string[0] = '\0';
    return &arr->a[arr->n++];
}

static char *cstr_array_putsn(const char *s, int l, cstr_array_t *arr)
{
    if (!arr || !s || l <= 0) return NULL;

    // oversized string: its own chunk.
    if (l > CSTR_SIZE) {
        cstr_t *chunk = cstr_array_expand(arr, (size_t) l);
        if (!chunk) return NULL;
        memcpy(chunk->string, s, (size_t) l);
        chunk->length = l;
        chunk->string[chunk->length] = '\0';
        return chunk->string;
    }

    // ensure at least one chunk exists.
    if (arr->n == 0)
        if (!cstr_array_expand(arr, 0)) return NULL;

    cstr_t *cur = &arr->a[arr->n - 1];

    // if current chunk cannot hold the whole string, start a new one.
    if (cur->length + l > CSTR_SIZE) {
        cur = cstr_array_expand(arr, CSTR_SIZE);
        if (!cur) return NULL;
    }

    memcpy(cur->string + cur->length, s, (size_t) l);
    cur->length += l;
    cur->string[cur->length] = '\0';
    return cur->string + cur->length - l;
}

static void cstr_array_destroy(cstr_array_t *arr)
{
    if (!arr) return;
    for (size_t i = 0; i < arr->n; i++)
        free(arr->a[i].string);
    free(arr->a);
    arr->n = arr->m = 0; arr->a = NULL;
}

static void sam_hdr_lite_free_text(sam_hdr_lite_t *hdr)
{
    if (!hdr) return;

    ks_free(&hdr->hd_text);

    cstr_array_destroy(&hdr->sq_text);
    cstr_array_destroy(&hdr->rg_text);
    cstr_array_destroy(&hdr->pg_text);
    cstr_array_destroy(&hdr->co_text);
}

// Build translation table for one input header; update @SQ, @RG, @PG text
// The 'text' field of the original header will be modified (in-place)
static bool sam_hdr_lite_build_tbl(sam_hdr_lite_t *bdr, sam_hdr_lite_t *hdr, int *refmap)
{
    if (!bdr || !hdr) return false;

    // temporary hash for building @PG translation table
    khash_t(c2c) *pmap = kh_init(c2c);
    if (!pmap) goto err;
    typedef struct { char *line; char *pg; char *pp; char *dup; int llen; } pg_info_t;
    array_t(pg_info_t) pinfo = {0,0,NULL};
    
    int ret;
    khiter_t k;
    char *line, *endl, *last, *tag, *dup, *pg, *pp;
    ssize_t i, new, old = 0;
    kstring_t tstr = {0,0,NULL}, *vn = &bdr->hd_text;
    pg_info_t *pg_line;
    bool res = false;

    // process @HD line
    if (!bdr->hd_text.l) {
        // no @HD line present; add the new one
        if (kputsn(hdr->hd_text.s, hdr->hd_text.l, &bdr->hd_text) < 0) goto err;
    }

    // process @SQ lines
    // we check refmap[-1] to see if we need to filter unused references
    // this might only be effective under '-m' - merge mode
    // which means no trimming of references and this header is a duplicate
    if (!refmap || (refmap && refmap[-1] <= 0)) {
        for (i = 0; i < hdr->sq_text.n; i++) {
            line = hdr->sq_text.a[i].string;
            last = line + hdr->sq_text.a[i].length;
            while (line < last) {
                // get the end of line
                endl = get_eol(line, last);
                // this is built inplace
                tag = get_tag(line+4, endl, "SN");
                if (!tag) goto err;
                // check if this reference is in use
                if (refmap && !refmap[old]) {
                    old++;
                    line = endl + 1; // next line
                    continue; // skip unused reference
                }
                k = kh_get(c2i, bdr->sq_map, tag);
                if (k == kh_end(bdr->sq_map)) { // new reference
                    dup = cstr_array_putsn(line, endl-line+1, &bdr->sq_text);
                    if (!dup) goto err;
                    dup += (tag - line); // point to SN value
                    new = bdr->n_targets++;
                    k = kh_put(c2i, bdr->sq_map, dup, &ret);
                    if (ret < 0) goto err;
                    kh_val(bdr->sq_map, k) = new;
                } else new = kh_val(bdr->sq_map, k); // existing reference
                if (refmap) refmap[old] = new;
                old++;
                line = endl + 1; // next line
            }
        }
    }

    // process @RG lines
    for (i = 0; i < hdr->rg_text.n; i++) {
        line = hdr->rg_text.a[i].string;
        last = line + hdr->rg_text.a[i].length;
        while (line < last) {
            // get the end of line
            endl = get_eol(line, last);
            // this is built inplace
            tag = get_tag(line+4, endl, "ID");
            if (!tag) goto err;
            if (kh_get(cset, bdr->rg_ids, tag) == kh_end(bdr->rg_ids)) {
                // append @RG line
                dup = cstr_array_putsn(line, endl-line+1, &bdr->rg_text);
                if (!dup) goto err;
                dup += (tag - line); // point to ID value
                k = kh_put(cset, bdr->rg_ids, dup, &ret);
                if (ret < 0) goto err;
            }
            line = endl + 1; // next line
        }
    }

    // process @PG lines
    for (i = 0; i < hdr->pg_text.n; i++) {
        line = hdr->pg_text.a[i].string;
        last = line + hdr->pg_text.a[i].length;
        while (line < last) {
            // get the end of line
            endl = get_eol(line, last);
            *endl = '\0'; // terminate line to avoid bs_strchrnul bounds issues
            // this is built inplace
            pg = get_tag(line+4, endl, "ID");
            if (!pg) goto err;
            pp = get_tag(line+4, endl, "PP");
            if (gen_unique_id(pg, bdr->pg_ids, true, &tstr) < 0) goto err;
            dup = strdup(tstr.s);
            if (!dup) goto err;
            kh_put(cset, bdr->pg_ids, dup, &ret); // store in global set
            if (ret < 0) goto err;
            k = kh_put(c2c, pmap, pg, &ret); // map old ID -> new ID
            if (ret < 0) goto err;
            kh_val(pmap, k) = dup;
            // save line info for second pass
            array_push(pg_info_t, pinfo, ((pg_info_t){line, pg, pp, dup, endl-line}));
            line = endl + 1; // next line
        }
    }
    // second pass: update @PG lines with new IDs & fixed PP
    for (i = 0; i < pinfo.n; i++) {
        pg_line = pinfo.a + i;
        // now build a new @PG line
        new = -1;
        tstr.l = 0;
        line = pg_line->line;
        last = line + pg_line->llen;
        while (line < last) {
            endl = bs_strchrnul(line, '\t');
            if (line+3 == pg_line->pg) { // ID: -> new ID
                new = tstr.l + 3; // remember where ID value starts
                if (kputsn("ID:", 3, &tstr) == EOF) goto err;
                if (kputs(kh_val(pmap, kh_get(c2c, pmap, pg_line->pg)), &tstr) == EOF) goto err;
            } else if (line+3 == pg_line->pp) { // PP: -> new PP
                if (kputsn("PP:", 3, &tstr) == EOF) goto err;
                k = kh_get(c2c, pmap, pg_line->pp);
                tag = k < kh_end(pmap)? kh_val(pmap, k) : NULL;
                if (!tag) {
                    warnmsg("PG line with ID:%s has a PP link to missing program '%s'", pg_line->pg, pg_line->pp);
                    tag = pg_line->pp; // fall back to old PP
                }
                if (kputs(tag, &tstr) == EOF) goto err;
            } else { // copy other tags verbatim
                if (kputsn(line, (int)(endl-line), &tstr) == EOF) goto err;
            }
            if (kputc('\0', &tstr) == EOF) goto err; // embedded null
            line = endl + 1;
        }
        dup = cstr_array_putsn(tstr.s, tstr.l, &bdr->pg_text);
        if (!dup) goto err;
        dup += new; // point to ID value
        k = kh_get(cset, bdr->pg_ids, dup);
        kh_key(bdr->pg_ids, k) = dup; // replace with stable copy
    }

    // process @CO lines
    for (i = 0; i < hdr->co_text.n; i++)
        if (!cstr_array_putsn(hdr->co_text.a[i].string, 
            hdr->co_text.a[i].length, &bdr->co_text)) 
            goto err;

    res = true;

err:
    if (pmap) kh_destroy(c2c, pmap);
    for (i = 0; i < pinfo.n; i++)
        free(pinfo.a[i].dup);
    free(pinfo.a);
    free(tstr.s);
    return ret;
}

/****************** END SAM Header Lite *******************/

/***************** Header Deduplication *******************/


static int sam_hdr_lite_trim(sam_hdr_lite_t *hdr, int *refs)
{
    if (!hdr) return -1;

    if (hdr->sq_map)      { free(hdr->sq_map); hdr->sq_map = NULL; }
    if (hdr->target_name) { free(hdr->target_name); hdr->target_name = NULL; }
    if (hdr->target_len)  { free(hdr->target_len); hdr->target_len = NULL; }
    
    // rebuild header @SQ lines
    int32_t i, s_idx = 0, n_refs = 0, n_targets = hdr->n_targets;
    char *line, *last, *endl, *fill, *tag;
    ssize_t llen;
    refs[-1] = -1; // unmapped
    // process @SQ lines
    for (i = 0; i < hdr->sq_text.n; i++) {
        line = hdr->sq_text.a[i].string;
        last = line + hdr->sq_text.a[i].length;
        fill = line;
        while (line < last) {
            // get the end of line
            endl = get_eol(line, last);
            // this is built inplace
            tag = get_tag(line+4, endl, "SN");
            if (!tag) return -1;
            // check if this reference is in use
            if (!refs[s_idx]) {
                s_idx++;
                line = endl + 1; // next line
                continue; // keep this reference
            }
            // move this @SQ line to the right place
            llen = endl - line + 1;
            memmove(fill, line, llen);
            fill += llen;
            refs[s_idx++] = n_refs++;
            line = endl + 1; // next line
        }
        hdr->sq_text.a[i].length = fill - hdr->sq_text.a[i].string;
        hdr->sq_text.a[i].string[hdr->sq_text.a[i].length] = '\0';
    }

    hdr->n_targets = n_refs;

    infomsg("BAM file reference sequences: %d", n_targets);
    infomsg("retained reference sequences: %d", n_refs);

    return n_refs;
}

static refseq_builder_t *refseq_builder_init(int capacity, size_t buff_size, char *fn_template)
{
    if (capacity < 1) capacity = 1;
    if (!buff_size) buff_size = 4ULL << 20; // 4MB

    char *tmpl = NULL;
    refseq_builder_t *rdb = NULL;
    
    rdb = (refseq_builder_t *)calloc(1, sizeof(refseq_builder_t));
    if (!rdb) goto fail;

    rdb->nref  = 0;
    rdb->refs  = (int **) calloc(capacity, sizeof(int*));
    rdb->sfpos = (fpos_t *) calloc(capacity, sizeof(fpos_t));
    rdb->index = (cint_array_t *) calloc(capacity, sizeof(cint_array_t));
    rdb->bsize = buff_size;
    rdb->fbuff = (char *) malloc(buff_size);
    rdb->shash = kh_init(c2i);
    if (!rdb->refs || !rdb->sfpos || !rdb->index || !rdb->fbuff || !rdb->shash)
        goto fail;

    if (fn_template) {
        tmpl = malloc(strlen(fn_template) + 1);
        if (!tmpl) goto fail;
        strcpy(tmpl, fn_template);
    } else {
        tmpl = malloc(24);
        if (!tmpl) goto fail;
        strcpy(tmpl, ".__bamsort.rdb.XXXXXX");
    }
    int fd = mkstemp(tmpl);
    if (fd == -1) goto fail;
    rdb->sfile = fdopen(fd, "w+b");
    if (!rdb->sfile) goto fail;
    unlink(tmpl); // anonymous
    
    free(tmpl);
    return rdb;

fail:
    free(tmpl);
    refseq_builder_destroy(rdb);

    return NULL;
}

#define FP_SLEN 100

static char *build_sq_fingerprint(cstr_array_t *sq_text, size_t sq_len) {
    if (!sq_text || sq_len == 0) return strdup("<empty>"); // special case for empty @SQ
    int t = sq_text->n - 1; // last @SQ line
    size_t head = sq_text->a[0].length < FP_SLEN ? sq_text->a[0].length : FP_SLEN;
    size_t tail = sq_text->a[t].length < FP_SLEN ? sq_text->a[t].length : FP_SLEN;
    char *hseq = sq_text->a[0].string;
    char *tseq = sq_text->a[t].length > FP_SLEN ? 
        (sq_text->a[t].string + (sq_text->a[t].length - tail)) : sq_text->a[t].string;
    // max digits for size_t (assume 20) + null
    size_t cap = head + tail + 20 + 1;
    char *fp = (char *) malloc(cap);
    if (!fp) return NULL;
    size_t off = 0;
    if (head) { memcpy(fp + off, hseq, head); off += head; }
    if (tail) { memcpy(fp + off, tseq, tail); off += tail; }
    off += (size_t) snprintf(fp + off, cap - off, "%zu", sq_len);
    fp[off] = 0;
    return fp;
}

static int refseq_builder_add(refseq_builder_t *rdb, sam_hdr_lite_t *hdr, int **_refmap)
{
    int ret = 0;
    int *refmap = NULL;
    char *fprint = NULL;
    
    if (!rdb || !hdr) goto fail;

    // build fingerprint
    cstr_array_t *sq_text = &hdr->sq_text;
    size_t i, sq_len = 0;

    for (i = 0; i < sq_text->n; i++) sq_len += sq_text->a[i].length;

    fprint = build_sq_fingerprint(&hdr->sq_text, sq_len);
    if (!fprint) goto fail;

    // search in hash
    cint_array_t *bucket = NULL;
    khiter_t k = kh_get(c2i, rdb->shash, fprint);
    if (k != kh_end(rdb->shash)) {
        int idx, base = kh_val(rdb->shash, k);
        bucket = &rdb->index[base];
        if (sq_len == 0) {
            // special case for empty @SQ
            refmap = rdb->refs[bucket->a[0]];
            ret = bucket->a[0] + 1;
            goto done;
        }
        size_t b, remain, chunk, bytes;
        char *string;
        bool equal;    
        for (b = 0; b < bucket->n; b++) {
            idx = bucket->a[b];
            // compare content
            if (fsetpos(rdb->sfile, &rdb->sfpos[idx]) != 0) goto fail;
            equal = true;
            for (i = 0; i < sq_text->n && equal; i++) {
                remain = sq_text->a[i].length;
                string = sq_text->a[i].string;
                while (remain) {
                    chunk = remain < rdb->bsize ? remain : rdb->bsize;
                    bytes = fread(rdb->fbuff, 1, chunk, rdb->sfile);
                    if (bytes != chunk || memcmp(rdb->fbuff, string, chunk) != 0) {
                        equal = false;
                        break;
                    }
                    string += chunk;
                    remain -= chunk;
                }
            }
            if (equal) {
                infomsg("found identical reference entry to findex %d", idx);
                refmap = rdb->refs[idx];
                ret = idx + 1;
                goto done;
            }
        }
    } else {
        int ret;
        k = kh_put(c2i, rdb->shash, fprint, &ret);
        if (ret < 0) goto fail;
        kh_val(rdb->shash, k) = rdb->nref;
        fprint = NULL; // hide it from free
        bucket = &rdb->index[rdb->nref];
    }
    // add new sequence content
    // either new fingerprint
    // or in rare cases hash collision with identical fingerprint
    // but different @SQ content: create new entry, add to bucket
    int newidx = rdb->nref;
    refmap = (int *) calloc(hdr->n_targets+1, sizeof(int)); // +1 for unmapped
    if (!refmap) goto fail;
    rdb->refs[newidx] = refmap;
    rdb->nref++;
    if (fseek(rdb->sfile, 0, SEEK_END) != 0) goto fail;
    if (fgetpos(rdb->sfile, &rdb->sfpos[newidx]) != 0) goto fail;
    for (i = 0; i < sq_text->n; i++) {
        size_t bytes = sq_text->a[i].length;
        if (fwrite(sq_text->a[i].string, 1, bytes, rdb->sfile) != bytes) goto fail;
    }
    array_push(int, *bucket, newidx);
    
    infomsg("found new reference entry findex %d", newidx);

done:
    free(fprint);
    if (_refmap) *_refmap = refmap;
    return ret; // 0 for not found, >0 for found (index+1)

fail:
    free(fprint);
    free(refmap);
    return -1; // error
}

static void refseq_builder_destroy(refseq_builder_t *rdb)
{
    if (!rdb) return;
    
    int i;
    if (rdb->sfile) fclose(rdb->sfile);
    for (i = 0; i < rdb->nref; i++)
        free(rdb->refs[i]);
    free(rdb->refs);
    free(rdb->sfpos);
    if (rdb->index) {
        for (i = 0; i < rdb->nref; i++)
            free(rdb->index[i].a);
        free(rdb->index);
    }
    free(rdb->fbuff);
    if (rdb->shash) {
        khiter_t k;
        for (k = kh_begin(rdb->shash); k != kh_end(rdb->shash); ++k)
            if (kh_exist(rdb->shash, k))
                free((char*)kh_key(rdb->shash, k));
        kh_destroy(c2i, rdb->shash);
    }
    free(rdb);
}

/***************** END Header Deduplication ***************/


/***************** SAM IOs Copied from HTSLIB *************/
#include "hts_endian.h"

KHASH_SET_INIT_INT(tag)

static inline int isprint_c(char c) { return isprint((unsigned char) c); }
const char *
hts_strprint(char *buf, size_t buflen, char quote, const char *s, size_t len)
{
    const char *slim = (len < SIZE_MAX)? &s[len] : NULL;
    char *t = buf, *bufend = buf + buflen;

    size_t qlen = quote? 1 : 0;
    if (quote) *t++ = quote;

    for (; slim? (s < slim) : (*s); s++) {
        char c;
        size_t clen;
        switch (*s) {
        case '\n': c = 'n'; clen = 2; break;
        case '\r': c = 'r'; clen = 2; break;
        case '\t': c = 't'; clen = 2; break;
        case '\0': c = '0'; clen = 2; break;
        case '\\': c = '\\'; clen = 2; break;
        default:
            c = *s;
            if (c == quote) clen = 2;
            else clen = isprint_c(c)? 1 : 4;
            break;
        }

        if (t-buf + clen + qlen >= buflen) {
            while (t-buf + 3 + qlen >= buflen) t--;
            if (quote) *t++ = quote;
            strcpy(t, "...");
            return buf;
        }

        if (clen == 4) {
            snprintf(t, bufend - t, "\\x%02X", (unsigned char) c);
            t += clen;
        }
        else {
            if (clen == 2) *t++ = '\\';
            *t++ = c;
        }
    }

    if (quote) *t++ = quote;
    *t = '\0';
    return buf;
}

static int sam_realloc_bam_data(bam1_t *b, size_t desired)
{
    uint32_t new_m_data;
    uint8_t *new_data;
    new_m_data = desired;
    kroundup32(new_m_data); // next power of 2
    new_m_data += 32; // reduces malloc arena migrations?
    if (new_m_data < desired) {
        errno = ENOMEM; // Not strictly true but we can't store the size
        return -1;
    }
#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (new_m_data > FUZZ_ALLOC_LIMIT) {
        errno = ENOMEM;
        return -1;
    }
#endif
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = realloc(b->data, new_m_data);
    } else {
        if ((new_data = malloc(new_m_data)) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) return -1;
    b->data = new_data;
    b->m_data = new_m_data;
    return 0;
}

static inline int possibly_expand_bam_data(bam1_t *b, size_t bytes) {
    size_t new_len = (size_t) b->l_data + bytes;

    if (new_len > INT32_MAX || new_len < bytes) { // Too big or overflow
        errno = ENOMEM;
        return -1;
    }
    if (new_len <= b->m_data) return 0;
    return sam_realloc_bam_data(b, new_len);
}

static inline uint64_t hts_str2uint(const char *in, char **end, int bits,
                                    int *failed) {
    uint64_t n = 0, limit = (bits < 64 ? (1ULL << bits) : 0) - 1;
    const unsigned char *v = (const unsigned char *) in;
    const unsigned int ascii_zero = '0'; // Prevents conversion to signed
    uint32_t fast = bits * 1000 / 3322 + 1; // log(10)/log(2) ~= 3.322
    unsigned int d;

    if (*v == '+')
        v++;

    while (--fast && *v>='0' && *v<='9')
        n = n*10 + *v++ - ascii_zero;

    if ((unsigned)(*v - ascii_zero) < 10 && !fast) {
        uint64_t limit_d_10 = limit / 10;
        uint64_t limit_m_10 = limit - 10 * limit_d_10;
        while ((d = *v - ascii_zero) < 10) {
            if (n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10)) {
                n = n*10 + d;
                v++;
            } else {
                do { v++; } while (*v - ascii_zero < 10);
                n = limit;
                *failed = 1;
                break;
            }
        }
    }

    *end = (char *)v;
    return n;
}

static inline int64_t hts_str2int(const char *in, char **end, int bits,
                                    int *failed) {
    uint64_t n = 0, limit = (1ULL << (bits - 1)) - 1;
    uint32_t fast = (bits - 1) * 1000 / 3322 + 1; // log(10)/log(2) ~= 3.322
    const unsigned char *v = (const unsigned char *) in;
    const unsigned int ascii_zero = '0'; // Prevents conversion to signed
    unsigned int d;

    int neg;
    switch(*v) {
    case '-':
        limit++;
        neg=1;
        v++;
        // See "dup" comment below
        while (--fast && *v>='0' && *v<='9')
            n = n*10 + *v++ - ascii_zero;
        break;

    case '+':
        v++;
        // fall through

    default:
        neg = 0;
        // dup of above.  This is somewhat unstable and mainly for code
        // size cheats to prevent instruction cache lines spanning 32-byte
        // blocks in the sam_parse_B_vals calling code.  It's been tested
        // on gcc7, gcc13, clang10 and clang16 with -O2 and -O3.  While
        // not exhaustive, this code duplication gives stable fast results
        // while a single copy does not.
        // (NB: system was "seq4d", so quite old)
        while (--fast && *v>='0' && *v<='9')
            n = n*10 + *v++ - ascii_zero;
        break;
    }

    // NB gcc7 is slow with (unsigned)(*v - ascii_zero) < 10,
    // while gcc13 prefers it.
    if (*v>='0' && !fast) { // rejects ',' and tab
        uint64_t limit_d_10 = limit / 10;
        uint64_t limit_m_10 = limit - 10 * limit_d_10;
        while ((d = *v - ascii_zero) < 10) {
            if (n < limit_d_10 || (n == limit_d_10 && d <= limit_m_10)) {
                n = n*10 + d;
                v++;
            } else {
                do { v++; } while (*v - ascii_zero < 10);
                n = limit;
                *failed = 1;
                break;
            }
        }
    }

    *end = (char *)v;

    return neg ? (int64_t)-n : (int64_t)n;
}

static inline int aux_type2size(uint8_t type)
{
    switch (type) {
    case 'A': case 'c': case 'C':
        return 1;
    case 's': case 'S':
        return 2;
    case 'i': case 'I': case 'f':
        return 4;
    case 'd':
        return 8;
    case 'Z': case 'H': case 'B':
        return type;
    default:
        return 0;
    }
}

static int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning) // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
{
    bam1_core_t *c = &b->core;

    // Bail out as fast as possible for the easy case
    uint32_t test_CG = BAM_CSOFT_CLIP | (c->l_qseq << BAM_CIGAR_SHIFT);
    if (c->n_cigar == 0 || test_CG != *bam_get_cigar(b))
        return 0;

    // The above isn't fool proof - we may have old CIGAR tags that aren't used,
    // but this is much less likely so do as a secondary check.
    if (c->tid < 0 || c->pos < 0)
        return 0;

    // Do we have a CG tag?
    uint8_t *CG = bam_aux_get(b, "CG");
    int saved_errno = errno;
    if (!CG) {
        if (errno != ENOENT) return -1;  // Bad aux data
        errno = saved_errno; // restore errno on expected no-CG-tag case
        return 0;
    }

    // Now we start with the serious work migrating CG to CIGAR
    uint32_t cigar_st, n_cigar4, CG_st, CG_en, ori_len = b->l_data,
        *cigar0, CG_len, fake_bytes;
    cigar0 = bam_get_cigar(b);
    fake_bytes = c->n_cigar * 4;
    if (CG[0] != 'B' || !(CG[1] == 'I' || CG[1] == 'i'))
        return 0; // not of type B,I
    CG_len = le_to_u32(CG + 2);
    // don't move if the real CIGAR length is shorter than the fake cigar length
    if (CG_len < c->n_cigar || CG_len >= 1U<<29) return 0;

    // move from the CG tag to the right position
    cigar_st = (uint8_t*)cigar0 - b->data;
    c->n_cigar = CG_len;
    n_cigar4 = c->n_cigar * 4;
    CG_st = CG - b->data - 2;
    CG_en = CG_st + 8 + n_cigar4;
    if (possibly_expand_bam_data(b, n_cigar4 - fake_bytes) < 0) return -1;
    // we need c->n_cigar-fake_bytes bytes to swap CIGAR to the right place
    b->l_data = b->l_data - fake_bytes + n_cigar4;
    // insert c->n_cigar-fake_bytes empty space to make room
    memmove(b->data + cigar_st + n_cigar4, b->data + cigar_st + fake_bytes, ori_len - (cigar_st + fake_bytes));
    // copy the real CIGAR to the right place; -fake_bytes for the fake CIGAR
    memcpy(b->data + cigar_st, b->data + (n_cigar4 - fake_bytes) + CG_st + 8, n_cigar4);
    if (ori_len > CG_en) // move data after the CG tag
        memmove(b->data + CG_st + n_cigar4 - fake_bytes, b->data + CG_en + n_cigar4 - fake_bytes, ori_len - CG_en);
    b->l_data -= n_cigar4 + 8; // 8: CGBI (4 bytes) and CGBI length (4)
    if (recal_bin)
        b->core.bin = hts_reg2bin(b->core.pos, bam_endpos(b), 14, 5);
    if (give_warning)
        hts_log_warning("%s encodes a CIGAR with %d operators at the CG tag", bam_get_qname(b), c->n_cigar);
    return 1;
}

/**********************
 *** SAM record I/O ***
 **********************/

// The speed of this code can vary considerably depending on minor code
// changes elsewhere as some of the tight loops are particularly prone to
// speed changes when the instruction blocks are split over a 32-byte
// boundary.  To protect against this, we explicitly specify an alignment
// for this function.  If this is insufficient, we may also wish to
// consider alignment of blocks within this function via
// __attribute__((optimize("align-loops=5"))) (gcc) or clang equivalents.
// However it's not very portable.
// Instead we break into separate functions so we can explicitly specify
// use __attribute__((aligned(32))) instead and force consistent loop
// alignment.
static inline int64_t grow_B_array(bam1_t *b, uint32_t *n, size_t size) {
    // Avoid overflow on 32-bit platforms, but it breaks BAM anyway
    if (*n > INT32_MAX*0.666) {
        errno = ENOMEM;
        return -1;
    }

    size_t bytes = (size_t)size * (size_t)(*n>>1);
    if (possibly_expand_bam_data(b, bytes) < 0) {
        errmsg("Out of memory");
        return -1;
    }

    (*n)+=*n>>1;
    return 0;
}

// This ensures that q always ends up at the next comma after
// reading a number even if it's followed by junk.  It
// prevents the possibility of trying to read more than n items.
#define skip_to_comma_(q) do { while (*(q) > '\t' && *(q) != ',') (q)++; } while (0)

HTS_ALIGN32
static char *sam_parse_Bc_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        *(b->data + b->l_data) = hts_str2int(q + 1, &q, 8, overflow);
        b->l_data++;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BC_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 1) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            *(b->data + b->l_data) = hts_str2uint(q + 1, &q, 8, overflow);
            b->l_data++;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bs_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        i16_to_le(hts_str2int(q + 1, &q, 16, overflow),
                  b->data + b->l_data);
        b->l_data += 2;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BS_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 2) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u16_to_le(hts_str2uint(q + 1, &q, 16, overflow),
                      b->data + b->l_data);
            b->l_data += 2;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bi_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        i32_to_le(hts_str2int(q + 1, &q, 32, overflow),
                  b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_BI_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        if (q[1] != '-') {
            u32_to_le(hts_str2uint(q + 1, &q, 32, overflow),
                      b->data + b->l_data);
            b->l_data += 4;
        } else {
            *overflow = 1;
            q++;
            skip_to_comma_(q);
        }
    }
    return q;
}

HTS_ALIGN32
static char *sam_parse_Bf_vals(bam1_t *b, char *q, uint32_t *nused,
                               uint32_t *nalloc, int *overflow) {
    while (*q == ',') {
        if ((*nused)++ >= (*nalloc)) {
            if (grow_B_array(b, nalloc, 4) < 0)
                return NULL;
        }
        float_to_le(strtod(q + 1, &q), b->data + b->l_data);
        b->l_data += 4;
    }
    return q;
}

HTS_ALIGN32
static int sam_parse_B_vals_r(char type, uint32_t nalloc, char *in,
                              char **end, bam1_t *b,
                              int *ctr) {
    // Protect against infinite recursion when dealing with invalid input.
    // An example string is "XX:B:C,-".  The lack of a number means min=0,
    // but it overflowed due to "-" and so we repeat ad-infinitum.
    //
    // Loop detection is the safest solution incase there are other
    // strange corner cases with malformed inputs.
    if (++(*ctr) > 2) {
        errmsg("Malformed data in B:%c array", type);
        return -1;
    }

    int orig_l = b->l_data;
    char *q = in;
    int32_t size;
    size_t bytes;
    int overflow = 0;

    size = aux_type2size(type);
    if (size <= 0 || size > 4) {
        errmsg("Unrecognized type B:%c", type);
        return -1;
    }

    // Ensure space for type + values.
    // The first pass through here we don't know the number of entries and
    // nalloc == 0.  We start with a small working set and then parse the
    // data, growing as needed.
    //
    // If we have a second pass through we do know the number of entries
    // and nalloc is already known.  We have no need to expand the bam data.
    if (!nalloc)
         nalloc=7;

    // Ensure allocated memory is big enough (for current nalloc estimate)
    bytes = (size_t) nalloc * (size_t) size;
    if (bytes / size != nalloc
        || possibly_expand_bam_data(b, bytes + 2 + sizeof(uint32_t))) {
        errmsg("Out of memory");
        return -1;
    }

    uint32_t nused = 0;

    b->data[b->l_data++] = 'B';
    b->data[b->l_data++] = type;
    // 32-bit B-array length is inserted later once we know it.
    int b_len_idx = b->l_data;
    b->l_data += sizeof(uint32_t);

    if (type == 'c') {
        if (!(q = sam_parse_Bc_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'C') {
        if (!(q = sam_parse_BC_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 's') {
        if (!(q = sam_parse_Bs_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'S') {
        if (!(q = sam_parse_BS_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'i') {
        if (!(q = sam_parse_Bi_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'I') {
        if (!(q = sam_parse_BI_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    } else if (type == 'f') {
        if (!(q = sam_parse_Bf_vals(b, q, &nused, &nalloc, &overflow)))
            return -1;
    }
    if (*q != '\t' && *q != '\0') {
        // Unknown B array type or junk in the numbers
        errmsg("Malformed B:%c", type);
        return -1;
    }
    i32_to_le(nused, b->data + b_len_idx);

    if (!overflow) {
        *end = q;
        return 0;
    } else {
        int64_t max = 0, min = 0, val;
        // Given type was incorrect.  Try to rescue the situation.
        char *r = q;
        q = in;
        overflow = 0;
        b->l_data = orig_l;
        // Find out what range of values is present
        while (q < r) {
            val = hts_str2int(q + 1, &q, 64, &overflow);
            if (max < val) max = val;
            if (min > val) min = val;
            skip_to_comma_(q);
        }
        // Retry with appropriate type
        if (!overflow) {
            if (min < 0) {
                if (min >= INT8_MIN && max <= INT8_MAX) {
                    return sam_parse_B_vals_r('c', nalloc, in, end, b, ctr);
                } else if (min >= INT16_MIN && max <= INT16_MAX) {
                    return sam_parse_B_vals_r('s', nalloc, in, end, b, ctr);
                } else if (min >= INT32_MIN && max <= INT32_MAX) {
                    return sam_parse_B_vals_r('i', nalloc, in, end, b, ctr);
                }
            } else {
                if (max < UINT8_MAX) {
                    return sam_parse_B_vals_r('C', nalloc, in, end, b, ctr);
                } else if (max <= UINT16_MAX) {
                    return sam_parse_B_vals_r('S', nalloc, in, end, b, ctr);
                } else if (max <= UINT32_MAX) {
                    return sam_parse_B_vals_r('I', nalloc, in, end, b, ctr);
                }
            }
        }
        // If here then at least one of the values is too big to store
        errmsg("Numeric value in B array out of allowed range");
        return -1;
    }
#undef skip_to_comma_
}

HTS_ALIGN32
static int sam_parse_B_vals(char type, char *in, char **end, bam1_t *b)
{
    int ctr = 0;
    uint32_t nalloc = 0;
    return sam_parse_B_vals_r(type, nalloc, in, end, b, &ctr);
}

static inline unsigned int parse_sam_flag(char *v, char **rv, int *overflow) {
    if (*v >= '1' && *v <= '9') {
        return hts_str2uint(v, rv, 16, overflow);
    }
    else if (*v == '0') {
        // handle single-digit "0" directly; otherwise it's hex or octal
        if (v[1] == '\t') { *rv = v+1; return 0; }
        else {
            unsigned long val = strtoul(v, rv, 0);
            if (val > 65535) { *overflow = 1; return 65535; }
            return val;
        }
    }
    else {
        // TODO implement symbolic flag letters
        *rv = v;
        return 0;
    }
}

// Parse tag line and append to bam object b.
// Shared by both SAM and FASTQ parsers.
//
// The difference between the two is how lenient we are to recognising
// non-compliant strings.  The FASTQ parser glosses over arbitrary
// non-SAM looking strings.
static inline int isspace_c(char c) { return isspace((unsigned char) c); }
static inline int aux_parse(char *start, char *end, bam1_t *b, int lenient,
                            khash_t(tag) *tag_whitelist) {
    int overflow = 0;
    int checkpoint;
    char logbuf[40];
    char *q = start, *p = end;

#define _parse_err(cond, ...)                   \
    do {                                        \
        if (cond) {                             \
            if (lenient) {                      \
                while (q < p && !isspace_c(*q))   \
                    q++;                        \
                while (q < p && isspace_c(*q))    \
                    q++;                        \
                b->l_data = checkpoint;         \
                goto loop;                      \
            } else {                            \
                errmsg(__VA_ARGS__);            \
                goto err_ret;                   \
            }                                   \
        }                                       \
    } while (0)

    while (q < p) loop: {
        char type;
        checkpoint = b->l_data;
        if (p - q < 5) {
            if (lenient) {
                break;
            } else {
                errmsg("Incomplete aux field");
                goto err_ret;
            }
        }
        _parse_err(q[0] < '!' || q[1] < '!', "invalid aux tag id");

        if (lenient && (q[2] | q[4]) != ':') {
            while (q < p && !isspace_c(*q))
                q++;
            while (q < p && isspace_c(*q))
                q++;
            continue;
        }

        if (tag_whitelist) {
            int tt = q[0]*256 + q[1];
            if (kh_get(tag, tag_whitelist, tt) == kh_end(tag_whitelist)) {
                while (q < p && *q != '\t')
                    q++;
                continue;
            }
        }

        // Copy over id
        if (possibly_expand_bam_data(b, 2) < 0) goto err_ret;
        memcpy(b->data + b->l_data, q, 2); b->l_data += 2;
        q += 3; type = *q++; ++q; // q points to value
        if (type != 'Z' && type != 'H') // the only zero length acceptable fields
            _parse_err(*q <= '\t', "incomplete aux field");

        // Ensure enough space for a double + type allocated.
        if (possibly_expand_bam_data(b, 16) < 0) goto err_ret;

        if (type == 'A' || type == 'a' || type == 'c' || type == 'C') {
            b->data[b->l_data++] = 'A';
            b->data[b->l_data++] = *q++;
        } else if (type == 'i' || type == 'I') {
            if (*q == '-') {
                int32_t x = hts_str2int(q, &q, 32, &overflow);
                if (x >= INT8_MIN) {
                    b->data[b->l_data++] = 'c';
                    b->data[b->l_data++] = x;
                } else if (x >= INT16_MIN) {
                    b->data[b->l_data++] = 's';
                    i16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'i';
                    i32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            } else {
                uint32_t x = hts_str2uint(q, &q, 32, &overflow);
                if (x <= UINT8_MAX) {
                    b->data[b->l_data++] = 'C';
                    b->data[b->l_data++] = x;
                } else if (x <= UINT16_MAX) {
                    b->data[b->l_data++] = 'S';
                    u16_to_le(x, b->data + b->l_data);
                    b->l_data += 2;
                } else {
                    b->data[b->l_data++] = 'I';
                    u32_to_le(x, b->data + b->l_data);
                    b->l_data += 4;
                }
            }
        } else if (type == 'f') {
            b->data[b->l_data++] = 'f';
            float_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(float);
        } else if (type == 'd') {
            b->data[b->l_data++] = 'd';
            double_to_le(strtod(q, &q), b->data + b->l_data);
            b->l_data += sizeof(double);
        } else if (type == 'Z' || type == 'H') {
            char *end = strchr(q, '\t');
            if (!end) end = q + strlen(q);
            _parse_err(type == 'H' && ((end-q)&1) != 0,
                       "hex field does not have an even number of digits");
            b->data[b->l_data++] = type;
            if (possibly_expand_bam_data(b, end - q + 1) < 0) goto err_ret;
            memcpy(b->data + b->l_data, q, end - q);
            b->l_data += end - q;
            b->data[b->l_data++] = '\0';
            q = end;
        } else if (type == 'B') {
            type = *q++; // q points to the first ',' following the typing byte
            _parse_err(*q && *q != ',' && *q != '\t',
                       "B aux field type not followed by ','");

            if (sam_parse_B_vals(type, q, &q, b) < 0)
                goto err_ret;
        } else _parse_err(1, "unrecognized type %s", hts_strprint(logbuf, sizeof logbuf, '\'', &type, 1));

        while (*q > '\t') { q++; } // Skip any junk to next tab
        q++;
    }

    _parse_err(!lenient && overflow != 0, "numeric value out of allowed range");
#undef _parse_err

    return 0;

err_ret:
    return -2;
}

static inline int bam_name2id_lite(sam_hdr_lite_t *h, const char *name)
{
    kh_c2i_t *sq_map = h->sq_map;
    khint_t k = kh_get(c2i, sq_map, name);
    return k == kh_end(sq_map) ? -1 : kh_val(sq_map, k);
}

int sam_parse1_lite(kstring_t *s, sam_hdr_lite_t *h, bam1_t *b)
{
#define _read_token(_p) (_p); do { char *tab = strchr((_p), '\t'); if (!tab) goto err_ret; *tab = '\0'; (_p) = tab + 1; } while (0)

#if HTS_ALLOW_UNALIGNED != 0 && ULONG_MAX == 0xffffffffffffffff

// Macro that operates on 64-bits at a time.
#define COPY_MINUS_N(to,from,n,l,failed)                        \
    do {                                                        \
        uint64_u *from8 = (uint64_u *)(from);                   \
        uint64_u *to8 = (uint64_u *)(to);                       \
        uint64_t uflow = 0;                                     \
        size_t l8 = (l)>>3, i;                                  \
        for (i = 0; i < l8; i++) {                              \
            to8[i] = from8[i] - (n)*0x0101010101010101UL;       \
            uflow |= to8[i];                                    \
        }                                                       \
        for (i<<=3; i < (l); ++i) {                             \
            to[i] = from[i] - (n);                              \
            uflow |= to[i];                                     \
        }                                                       \
        failed = (uflow & 0x8080808080808080UL) > 0;            \
    } while (0)

#else

// Basic version which operates a byte at a time
#define COPY_MINUS_N(to,from,n,l,failed) do {                \
        uint8_t uflow = 0;                                   \
        for (i = 0; i < (l); ++i) {                          \
            (to)[i] = (from)[i] - (n);                       \
            uflow |= (uint8_t) (to)[i];                      \
        }                                                    \
        failed = (uflow & 0x80) > 0;                         \
    } while (0)

#endif

#define _get_mem(type_t, x, b, l) if (possibly_expand_bam_data((b), (l)) < 0) goto err_ret; *(x) = (type_t*)((b)->data + (b)->l_data); (b)->l_data += (l)
#define _parse_err(cond, ...) do { if (cond) { errmsg(__VA_ARGS__); goto err_ret; } } while (0)
#define _parse_warn(cond, ...) do { if (cond) { warnmsg(__VA_ARGS__); } } while (0)

    uint8_t *t;

    char *p = s->s, *q;
    int i, overflow = 0;
    char logbuf[40];
    hts_pos_t cigreflen;
    bam1_core_t *c = &b->core;

    b->l_data = 0;
    memset(c, 0, 32);

    // qname
    q = _read_token(p);

    _parse_warn(p - q <= 1, "empty query name");
    _parse_err(p - q > 255, "query name too long");
    // resize large enough for name + extranul
    if (possibly_expand_bam_data(b, (p - q) + 4) < 0) goto err_ret;
    memcpy(b->data + b->l_data, q, p-q); b->l_data += p-q;

    c->l_extranul = (4 - (b->l_data & 3)) & 3;
    memcpy(b->data + b->l_data, "\0\0\0\0", c->l_extranul);
    b->l_data += c->l_extranul;

    c->l_qname = p - q + c->l_extranul;

    // flag
    c->flag = parse_sam_flag(p, &p, &overflow);
    if (*p++ != '\t') goto err_ret; // malformated flag

    // chr
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(h->n_targets == 0, "no SQ lines present in the header");
        c->tid = bam_name2id_lite(h, q);
        _parse_err(c->tid < -1, "failed to parse header");
        _parse_warn(c->tid < 0, "unrecognized reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    } else c->tid = -1;

    // pos
    c->pos = hts_str2uint(p, &p, 62, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->pos < 0 && c->tid >= 0) {
        _parse_warn(1, "mapped query cannot have zero coordinate; treated as unmapped");
        c->tid = -1;
    }
    if (c->tid < 0) c->flag |= BAM_FUNMAP;

    // mapq
    c->qual = hts_str2uint(p, &p, 8, &overflow);
    if (*p++ != '\t') goto err_ret;
    // cigar
    if (*p != '*') {
        uint32_t *cigar = NULL;
        int old_l_data = b->l_data;
        int n_cigar = bam_parse_cigar(p, &p, b);
        if (n_cigar < 1 || *p++ != '\t') goto err_ret;
        cigar = (uint32_t *)(b->data + old_l_data);

        // can't use bam_endpos() directly as some fields not yet set up
        cigreflen = (!(c->flag&BAM_FUNMAP))? bam_cigar2rlen(c->n_cigar, cigar) : 1;
        if (cigreflen == 0) cigreflen = 1;
    } else {
        _parse_warn(!(c->flag&BAM_FUNMAP), "mapped query must have a CIGAR; treated as unmapped");
        c->flag |= BAM_FUNMAP;
        q = _read_token(p);
        cigreflen = 1;
    }
    _parse_err(HTS_POS_MAX - cigreflen <= c->pos,
               "read ends beyond highest supported position");
    c->bin = hts_reg2bin(c->pos, c->pos + cigreflen, 14, 5);
    // mate chr
    q = _read_token(p);
    if (strcmp(q, "=") == 0) {
        c->mtid = c->tid;
    } else if (strcmp(q, "*") == 0) {
        c->mtid = -1;
    } else {
        c->mtid = bam_name2id_lite(h, q);
        _parse_err(c->mtid < -1, "failed to parse header");
        _parse_warn(c->mtid < 0, "unrecognized mate reference name %s; treated as unmapped", hts_strprint(logbuf, sizeof logbuf, '"', q, SIZE_MAX));
    }
    // mpos
    c->mpos = hts_str2uint(p, &p, 62, &overflow) - 1;
    if (*p++ != '\t') goto err_ret;
    if (c->mpos < 0 && c->mtid >= 0) {
        _parse_warn(1, "mapped mate cannot have zero coordinate; treated as unmapped");
        c->mtid = -1;
    }
    // tlen
    c->isize = hts_str2int(p, &p, 63, &overflow);
    if (*p++ != '\t') goto err_ret;
    _parse_err(overflow, "number outside allowed range");
    // seq
    q = _read_token(p);
    if (strcmp(q, "*")) {
        _parse_err(p - q - 1 > INT32_MAX, "read sequence is too long");
        c->l_qseq = p - q - 1;
        hts_pos_t ql = bam_cigar2qlen(c->n_cigar, (uint32_t*)(b->data + c->l_qname));
        _parse_err(c->n_cigar && ql != c->l_qseq, "CIGAR and query sequence are of different length");
        i = (c->l_qseq + 1) >> 1;
        _get_mem(uint8_t, &t, b, i);

        unsigned int lqs2 = c->l_qseq&~1, i;
        for (i = 0; i < lqs2; i+=2)
            t[i>>1] = (seq_nt16_table[(unsigned char)q[i]] << 4) | seq_nt16_table[(unsigned char)q[i+1]];
        for (; i < c->l_qseq; ++i)
            t[i>>1] = seq_nt16_table[(unsigned char)q[i]] << ((~i&1)<<2);
    } else c->l_qseq = 0;
    // qual
    _get_mem(uint8_t, &t, b, c->l_qseq);
    if (p[0] == '*' && (p[1] == '\t' || p[1] == '\0')) {
        memset(t, 0xff, c->l_qseq);
        p += 2;
    } else {
        int failed = 0;
        _parse_err(s->l - (p - s->s) < c->l_qseq
                   || (p[c->l_qseq] != '\t' && p[c->l_qseq] != '\0'),
                   "SEQ and QUAL are of different length");
        COPY_MINUS_N(t, p, 33, c->l_qseq, failed);
        _parse_err(failed, "invalid QUAL character");
        p += c->l_qseq + 1;
    }

    // aux
    if (aux_parse(p, s->s + s->l, b, 0, NULL) < 0)
        goto err_ret;

    if (bam_tag2cigar(b, 1, 1) < 0)
        return -2;
    return 0;

#undef _parse_warn
#undef _parse_err
#undef _get_mem
#undef _read_token
err_ret:
    return -2;
}

/*
 * Convert a nibble encoded BAM sequence to a string of bases.
 *
 * We do this 2 bp at a time for speed. Equiv to:
 *
 * for (i = 0; i < len; i++)
 *    seq[i] = seq_nt16_str[bam_seqi(nib, i)];
 */
static inline void nibble2base_default(uint8_t *nib, char *seq, int len) {
    static const char code2base[512] =
        "===A=C=M=G=R=S=V=T=W=Y=H=K=D=B=N"
        "A=AAACAMAGARASAVATAWAYAHAKADABAN"
        "C=CACCCMCGCRCSCVCTCWCYCHCKCDCBCN"
        "M=MAMCMMMGMRMSMVMTMWMYMHMKMDMBMN"
        "G=GAGCGMGGGRGSGVGTGWGYGHGKGDGBGN"
        "R=RARCRMRGRRRSRVRTRWRYRHRKRDRBRN"
        "S=SASCSMSGSRSSSVSTSWSYSHSKSDSBSN"
        "V=VAVCVMVGVRVSVVVTVWVYVHVKVDVBVN"
        "T=TATCTMTGTRTSTVTTTWTYTHTKTDTBTN"
        "W=WAWCWMWGWRWSWVWTWWWYWHWKWDWBWN"
        "Y=YAYCYMYGYRYSYVYTYWYYYHYKYDYBYN"
        "H=HAHCHMHGHRHSHVHTHWHYHHHKHDHBHN"
        "K=KAKCKMKGKRKSKVKTKWKYKHKKKDKBKN"
        "D=DADCDMDGDRDSDVDTDWDYDHDKDDDBDN"
        "B=BABCBMBGBRBSBVBTBWBYBHBKBDBBBN"
        "N=NANCNMNGNRNSNVNTNWNYNHNKNDNBNN";

    int i, len2 = len/2;
    seq[0] = 0;

    for (i = 0; i < len2; i++)
        // Note size_t cast helps gcc optimiser.
        memcpy(&seq[i*2], &code2base[(size_t)nib[i]*2], 2);

    if ((i *= 2) < len)
        seq[i] = seq_nt16_str[bam_seqi(nib, i)];
}

#if defined HAVE_ATTRIBUTE_CONSTRUCTOR && \
    ((defined __x86_64__ && defined HAVE_ATTRIBUTE_TARGET_SSSE3 && defined HAVE_BUILTIN_CPU_SUPPORT_SSSE3) || \
     (defined __ARM_NEON))
#define BUILDING_SIMD_NIBBLE2BASE
#endif

static inline void nibble2base(uint8_t *nib, char *seq, int len) {
#ifdef BUILDING_SIMD_NIBBLE2BASE
    extern void (*htslib_nibble2base)(uint8_t *nib, char *seq, int len);
    htslib_nibble2base(nib, seq, len);
#else
    nibble2base_default(nib, seq, len);
#endif
}

// With gcc, -O3 or -ftree-loop-vectorize is really key here as otherwise
// this code isn't vectorised and runs far slower than is necessary (even
// with the restrict keyword being used).
static inline void HTS_OPT3
add33(uint8_t *a, const uint8_t * b, int32_t len) {
    uint32_t i;
    for (i = 0; i < len; i++)
        a[i] = b[i]+33;
}

static int sam_format1_lite(const sam_hdr_lite_t *h, const bam1_t *b, kstring_t *str)
{
    int i, r = 0;
    uint8_t *s, *end;
    const bam1_core_t *c = &b->core;

    str->l = 0;

    if (c->l_qname == 0)
        return -1;
    r |= kputsn_(bam_get_qname(b), c->l_qname-1-c->l_extranul, str);
    r |= kputc_('\t', str); // query name
    r |= kputw(c->flag, str); r |= kputc_('\t', str); // flag
    if (c->tid >= 0) { // chr
        r |= kputs(h->target_name[c->tid] , str);
        r |= kputc_('\t', str);
    } else r |= kputsn_("*\t", 2, str);
    r |= kputll(c->pos + 1, str); r |= kputc_('\t', str); // pos
    r |= kputw(c->qual, str); r |= kputc_('\t', str); // qual
    if (c->n_cigar) { // cigar
        uint32_t *cigar = bam_get_cigar(b);
        for (i = 0; i < c->n_cigar; ++i) {
            r |= kputw(bam_cigar_oplen(cigar[i]), str);
            r |= kputc_(bam_cigar_opchr(cigar[i]), str);
        }
    } else r |= kputc_('*', str);
    r |= kputc_('\t', str);
    if (c->mtid < 0) r |= kputsn_("*\t", 2, str); // mate chr
    else if (c->mtid == c->tid) r |= kputsn_("=\t", 2, str);
    else {
        r |= kputs(h->target_name[c->mtid], str);
        r |= kputc_('\t', str);
    }
    r |= kputll(c->mpos + 1, str); r |= kputc_('\t', str); // mate pos
    r |= kputll(c->isize, str); r |= kputc_('\t', str); // template len
    if (c->l_qseq) { // seq and qual
        uint8_t *s = bam_get_seq(b);
        if (ks_resize(str, str->l+2+2*c->l_qseq) < 0) goto mem_err;
        char *cp = str->s + str->l;

        // Sequence, 2 bases at a time
        nibble2base(s, cp, c->l_qseq);
        cp[c->l_qseq] = '\t';
        cp += c->l_qseq+1;

        // Quality
        s = bam_get_qual(b);
        i = 0;
        if (s[0] == 0xff) {
            cp[i++] = '*';
        } else {
            add33((uint8_t *)cp, s, c->l_qseq); // cp[i] = s[i]+33;
            i = c->l_qseq;
        }
        cp[i] = 0;
        cp += i;
        str->l = cp - str->s;
    } else r |= kputsn_("*\t*", 3, str);

    s = bam_get_aux(b); // aux
    end = b->data + b->l_data;

    while (end - s >= 4) {
        r |= kputc_('\t', str);
        if ((s = (uint8_t *)sam_format_aux1(s, s[2], s+3, end, str)) == NULL)
            goto bad_aux;
    }
    r |= kputsn("", 0, str); // nul terminate
    if (r < 0) goto mem_err;

    return str->l;

 bad_aux:
    errmsg("Corrupted aux data for read %.*s flag %d",
                  b->core.l_qname, bam_get_qname(b), b->core.flag);
    errno = EINVAL;
    return -1;

 mem_err:
    errmsg("Out of memory");
    errno = ENOMEM;
    return -1;
}

/*********************** END SAM IO ***********************/
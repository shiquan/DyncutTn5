#include "utils.h"
#include "number.h"
#include "kseq.h"
#include "kstring.h"
#include "thread_pool.h"
#include <zlib.h>

KSEQ_INIT(gzFile, gzread)

static char * program_name =  "TN5dyncut";

#define BASEA 1
#define BASEC 2
#define BASEG 4
#define BASET 8

unsigned char *base_tab_init()
{
    unsigned char *t = malloc(sizeof(char)*256);
    memset(t, 0, sizeof(char)*256);
    t['A'] = BASEA;
    t['C'] = BASEC;
    t['G'] = BASEG;
    t['T'] = BASET; 
    t['a'] = BASEA;
    t['c'] = BASEC;
    t['g'] = BASEG;
    t['t'] = BASET;
    return t;
}

// Transposase recognition sequences
// 19bp Mosaic Ends: CTGTCTCTTATACACATCT
const char *me = "CTGTCTCTTATACACATCT";

#define MINI_SKIP 1
#define TRIMMED   2
#define POLLUTION 3

static int check_name(kstring_t name1, kstring_t name2)
{
    char *s1 = name1.s;
    char *s2 = name2.s;
    size_t n;
    for(n = 0; n < name1.l - 1; ++n, ++s1, ++s2) {
        if (*s1 != *s2) {
            return 0;
        }
    }
    return 1;
}
int countbits(uint64_t x)
{
    int i;
    int l = 0;
    for ( i = 0; i < 16; ) {
        if ( ((x>>(i*4))&0xf) > 0 ) l++;
    }
    return l;
}
int usage()
{
    fprintf(stderr,
            "A program to dynamic cut TN5 mosaic ends on the reads.\n"
            "Usage: TN5dyncut [options] read1.fq.gz [read2.fq.gz]\n"
            "\n"
            "Options:\n"
            "  -1 [output1.fq.gz]   Output fastq 1 file or smart pairing file (ignore fastq 2).\n"
            "  -2 [output2.fq.gz]   Output fastq 2 file.\n"
            "  -S                   Trim PE reads in SE mode.\n"
            "  -m [1]               Allowed mismatches on the adaptor.\n"            
            "  -l [20]              Minimal fragment to keep.\n"
            "  -tail [0]            Trimmed both ends if no adaptor detected.\n"
            "  -adaptor [seq]       Adaptor sequences. Default is 19bp mosaic ends.\n"
            "  -report [report.txt] Export report summary.\n"
            "  -d                   Drop reads if reversed adaptor detected.\n"
            "  -t [1]               Threads.\n"
            "\n"
            "Notes:\n"
            "The program is designed to trim TN5 transposase introduced mosaic ends and\n"
            "adaptor populations. But it should also be work for other kinds of transposases\n"
            "by define a specific adaptor sequence.\n"
            "Demo 1:\n"
            " TN5dyncut -tail 5 reads1.fq.gz reads2.fq.gz | bwa mem -Cpt8 ref.fa - \\\n"
            " | samtools view -Sb -|sambamba sort /dev/stdin -|samtools rmdup /dev/stdin aln.bam\n"
            "Demo 2:\n"
            " TN5dyncut -tail 5 reads1.fq.gz reads2.fq.gz -1 out1.fq.gz -2 out2.fq.gz &&\\\n"
            " bowtie2 -x ref.fa -1 out1.fq.gz -2 out2.fq.gz | samtools view -Sb - | ...\n"
        );
    return 1;
}

struct encode {
    int l;
    uint64_t x;
};
struct encode *str_encode(const char *s, unsigned char *tab)
{
    int i;
    int l = strlen(s);
    struct encode *x = malloc(sizeof(*x));
    x->x = 0;
    x->l = 0;
    // only encode first 16 bases or shorter
    x->l = l/4 > 16 ? 16 : l/4;
    int mask = (1ULL<<l*4)-1;
    for (i = 0; i < x->l; ++i) {
        x->x &= ((uint64_t)tab[s[i]] << i*4);
    }
    return x;
}
struct bseq {
    int flag; // flag for skip reasons
    // read 1
    char *n0, *s0, *q0;
    int l0;
    // read 2
    char *s1, *q1;
    int l1;
};
struct bseq_pool {
    struct bseq *s;
    int n, m;
    struct args *opts;
};
struct bseq_pool *bseq_pool_init()
{
    struct bseq_pool *p = malloc(sizeof(*p));
    memset(p, 0, sizeof(*p));
    return p;
}
void bseq_pool_destroy(struct bseq_pool *p)
{
    int i;
    for ( i = 0; i < p->n; ++i ) {
        struct bseq *b = &p->s[i];
        if (b->l0) {
            free(b->n0);
            free(b->s0);
            if (b->q0) free(b->q0);
        }
        if (b->l1) {
            //free(b->n1);
            free(b->s1);
            if (b->q1) free(b->q1);
        }
    }
    if (p->m > 0) free(p->s);
    free(p);
}
struct bseq_pool *bseq_read(kseq_t *k1, kseq_t *k2, int chunk_size, int pe)
{
    struct bseq_pool *p = bseq_pool_init();
    int size = 0;
    
    if ( pe == 0 ) {
        while (kseq_read(k1) >= 0 ) {
            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct bseq));
            }
            struct bseq *s = &p->s[p->n];
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            size += s->l0;
            p->n++;
            if ( size >= chunk_size ) break;
        }
    }
    else {
        while (kseq_read(k1) >= 0 && kseq_read(k2) >= 0 ) {
            if ( check_name(k1->name, k2->name) ) error("Inconsistance paired read names. %s vs %s.", k1->name.s, k2->name.s);
            //if ( k1->seq.l != k2->seq.l ) error("Inconsistant PE read length, %s.", k1->name.s);
            
            struct bseq *s;
            if (p->n >= p->m) {
                p->m = p->m ? p->m<<1 : 256;
                p->s = realloc(p->s, p->m*sizeof(struct bseq));
            }
            s = &p->s[p->n];
            s->n0 = strdup(k1->name.s);
            s->s0 = strdup(k1->seq.s);
            s->q0 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l0 = k1->seq.l;
            size += s->l0;

            //s->n1 = strdup(k1->name.s);
            s->s1 = strdup(k1->seq.s);
            s->q1 = k1->qual.l? strdup(k1->qual.s) : 0;
            s->l1 = k1->seq.l;
            size += s->l1;
            
            p->n++;
            if ( size >= chunk_size ) break;            
        }
    }
    if ( p->n == 0 ) {
        bseq_pool_destroy(p);
        return NULL;
    }
    return p;
}

struct args {
    const char *input_fname1;
    const char *input_fname2;
    const char *output_fname1;
    const char *output_fname2;
    const char *report_fname;
    
    int n_thread;
    int se_mode;
    int mismatch;
    int tail_3;
    int mini_frag;
    int drop_poll;    
    const char *adaptor;
    unsigned char *base_tab;

    gzFile r1_fp;
    gzFile r2_fp;
    gzFile r1_out;
    gzFile r2_out;
    FILE *report_fp;

    kseq_t *k1;
    kseq_t *k2;

    int is_pe;
    int chunk_size;
    struct encode *me_or_ada;
} args = {
    .input_fname1 = NULL,
    .input_fname2 = NULL,
    .output_fname1 = NULL,
    .output_fname2 = NULL,
    .report_fname = NULL,
    .n_thread = 1,
    .se_mode = 0,
    .mismatch = 1,
    .mini_frag = 18,
    .tail_3 = 0,
    .drop_poll = 0,
    .adaptor = NULL,
    .base_tab = NULL,
    .k1 = NULL,
    .k2 = NULL,
    .report_fp = NULL,
    .r1_fp = NULL,
    .r2_fp = NULL,
    .r1_out = NULL,
    .r2_out = NULL,
    .me_or_ada = NULL,
    .chunk_size = 10000000,
    .is_pe = 0,
};
void memory_release()
{
    free(args.base_tab);
    free(args.me_or_ada);
    if ( args.k1 ) {
        kseq_destroy(args.k1);
        gzclose(args.r1_fp);
    }
    if ( args.k2 ) {
        kseq_destroy(args.k2);
        gzclose(args.r2_fp);
    }
    if ( args.r1_out ) gzclose(args.r1_out);
    if ( args.r2_out ) gzclose(args.r2_out);
    if ( args.report_fp ) fclose(args.report_fp);    
}
int parse_args(int argc, char **argv)
{
    if ( argc == 1 ) {
        fprintf(stderr, "* TN5dyncut reads1.fq.gz reads2.fq.gz\n"
                "* Use -h for more information.\n"
            );
        return 1;
    }
    
    int i;
    const char *mis = 0;
    const char *t3 = 0;
    const char *mini = 0;
    const char *thread = 0;
    
    for ( i = 1; i < argc; ) {

        const char *a = argv[i++];
        const char **var = 0;
        
        if ( strcmp(a, "-h") == 0 || strcmp(a, "--help") == 0 ) return usage();

        if ( strcmp(a, "-1") == 0 ) var = &args.output_fname1;
        else if ( strcmp(a, "-2") == 0 ) var = &args.output_fname2;
        else if ( strcmp(a, "-adaptor") == 0 ) var = &args.adaptor;
        else if ( strcmp(a, "-m") == 0 ) var = &mis;
        else if ( strcmp(a, "-l") == 0 ) var = &mini;
        else if ( strcmp(a, "-t") == 0 ) var = &thread;
        else if ( strcmp(a, "-tail") == 0 ) var = &t3;
        else if ( strcmp(a, "-S") == 0 ) {
            args.se_mode = 1;
            continue;
        }
        else if ( strcmp(a, "-d") == 0 ) {
            args.drop_poll = 1;
            continue;
        }

        if ( var != 0 ) {
            if ( i == argc ) error("Miss an argument after %s.", a);
            *var = argv[i++];
            continue;
        }

        if ( a[0] == '-' && a[1] ) error("Unknown parameter. %s", a);
        
        if ( args.input_fname1 == 0 ) {
            args.input_fname1 = a;
            continue;
        }

        if ( args.input_fname2 == 0 ) {
            args.input_fname2 = a;
            continue;
        }

        error("Unknown argument: %s, use -h see help information.", a);
    }

    if (mis) args.mismatch = str2int((char*)mis);
    if (mini) args.mini_frag = str2int((char*)mini);
    if (thread) args.n_thread = str2int((char*)thread);
    if (t3) args.tail_3 = str2int((char*)t3);

    args.base_tab = base_tab_init();
    
    args.me_or_ada = args.adaptor == NULL ? str_encode(me, args.base_tab) : str_encode(args.adaptor, args.base_tab);

    if ( args.input_fname1 == NULL && (!isatty(fileno(stdin))) )
        args.input_fname1 = "-";
    if ( args.input_fname1 == NULL ) error("Fastq file(s) must be set!");
    
    args.r1_fp = gzopen(args.input_fname1, "r");
    if ( args.r1_fp == NULL ) error("Failed to open %s.", args.input_fname1);

    if ( args.input_fname2 == NULL ) {
        args.is_pe = 0;        
    }
    else {
        args.r2_fp = gzopen(args.input_fname2, "r");
        if ( args.r2_fp == NULL ) error("Failed to open %s.", args.input_fname2);
        args.is_pe = 1;
    }
    
    return 0;
}
int find_sequence_adaptor(char *s, struct encode *a, int m, unsigned char *tab)
{
    int l, n = 0;
    int i;
    l = strlen(s);    
    uint64_t x = 0;
    int mask = (1ULL << a->l*4) -1;
    for ( i = 0; i < l; ++i) {
        int c = tab[s[i]];
        x = (x<<4|c)&mask;
        if ( ++n >= a->l) {
            if ( x == a->x ||countbits(x&a->x) >= a->l - m) return i;
        }
    }
    // tails
    int mis = m;
    int x1;
    for ( i = 1; i < a->l; ++i ) {
        mis = mis <= 0 ? 0 : --mis;
        x >>=4;
        int l1 = a->l - i;
        mask = (1ULL<<(a->l-i)*4)-1;
        x1 = a->x & mask;
        if ( x==x1 || countbits(x&x1)>= l1 -m ) return l - a->l + i;
    }
    return -1;
}   
void *trim_core(void *_p, int idx)
{
    int i;
    struct bseq_pool *p = (struct bseq_pool*)_p;
    struct args *opts = p->opts;
    for ( i = 0; i < p->n; ++p ) {
        struct bseq *b = &p->s[i];
        int l0, l1 = -1;
        l0 = find_sequence_adaptor(b->s0, opts->me_or_ada, opts->mismatch, opts->base_tab);
        
        if ( opts->is_pe ) {
            l1 = find_sequence_adaptor(b->s1, opts->me_or_ada, opts->mismatch, opts->base_tab);
            
            // consider PE reads are not same length, treat as SEs
            if ( opts->se_mode ) {
                if ( l1 == -1 ) {
                    if ( opts->tail_3 > 0 && b->l1 > opts->tail_3 ) b->l1 -= opts->tail_3;
                }
                else {
                    b->l1 = l1;
                }
            }
            else {
                if ( l0 == -1 ) {                
                    if (l1 > 0 ) {
                        l0 = l0 > l1 ? l1 : l0;
                        b->l1 = l1;
                    }
                    else {
                        if ( opts->tail_3 > 0 && b->l1 > opts->tail_3 ) b->l1 -= opts->tail_3;
                    }
                }
                else {
                    if (l0 < l1 ) l1 = l0;
                    else l0 = l1;
                }
            }
        }
        
        if ( l0 == -1 ) {            
            if ( opts->tail_3 > 0 && b->l0 > opts->tail_3 ) b->l0 -= opts->tail_3;
        }
        else {
            b->l0 = l0;
            if ( b->l0 < opts->mini_frag ) b->flag = MINI_SKIP;
            else b->flag = TRIMMED;
        }
    }
    return p;
}
int write_out(struct bseq_pool *p)
{
    int i;
    int mode;
    struct args *opts = p->opts;
    if ( opts->is_pe == 1 ) {
        if (opts->r1_out && opts->r2_out ) mode = 1;
        else if (opts->r1_out ) mode = 2;
        else mode = 3;
    }
    else {
        if (opts->r1_out == NULL ) mode = 4;
        else mode = 5;
    }
         
    for ( i = 0; i < p->n; ++i ) {
        struct bseq *b = &p->s[i];
        if ( b->flag != MINI_SKIP ) {
            kstring_t str1 = {0,0,0};
            kstring_t str2 = {0,0,0};
            kputs(b->n0, &str1);            
            kputc('\n', &str1);
            kputsn(b->s0, b->l0, &str1);
            kputc('\n', &str1);
            if (b->q0) {
                kputc('+', &str1);
                kputc('\n', &str1);
                kputsn(b->q0, b->l0, &str1);
                kputc('\n', &str2);
            }
            if ( opts->is_pe ) {
                kputs(b->n0, &str2);
                kputc('\n', &str2);
                kputsn(b->s1, b->l1, &str2);
                kputc('\n', &str2);
                if (b->q1) {
                    kputc('+', &str2);
                    kputc('\n', &str2);
                    kputsn(b->q1, b->l1, &str2);
                    kputc('\n', &str2);
                }                
            }

            if ( mode == 1 ) {
                gzwrite(opts->r1_out, str1.s, str1.l);
                gzwrite(opts->r2_out, str2.s, str2.l);                
            }
            else if ( mode == 2 ) {
                gzwrite(opts->r1_out, str1.s, str1.l);
                gzwrite(opts->r1_out, str2.s, str2.l);
            }
            else if ( mode == 3 ) {
                puts(str1.s);
                puts(str2.s);
            }
            else if ( mode == 4 ) {
                puts(str1.s);
            }
            else if ( mode == 5 ) {
                gzwrite(opts->r1_out, str1.s, str1.l);
            }
            free(str1.s);
            if (str2.m) free(str2.s);
        }
    }
    return 0;
}
int trim_adap_light()
{
    do {
        struct bseq_pool *p = bseq_read(args.k1, args.k2, args.chunk_size, args.k2 != NULL);
        if (p == NULL) break;
        int idx;
        p->opts = &args;
        trim_core(p, idx);
        write_out(p);
        bseq_pool_destroy(p);
    }
    while(1);
    return 0;
}

int main(int argc, char **argv)
{
    if ( parse_args(argc, argv) == 1 )
        return 1;
    if ( args.n_thread == 1 ) trim_adap_light();
    else {
        args.n_thread -= 1;
        struct thread_pool *p = thread_pool_init(args.n_thread);
        struct thread_pool_process *q = thread_pool_process_init(p, args.n_thread*2, 0);
        struct thread_pool_result *r;

        for (;;) {
            struct bseq_pool *b = bseq_read(args.k1, args.k2, args.chunk_size, args.k2 != NULL);
            if (b == NULL) break;
            b->opts = &args;
            int block;
            do {
                block = thread_pool_dispatch2(p, q, trim_core, b, 1);
                if ((r = thread_pool_next_result(q))) {
                    struct bseq_pool *d = (struct bseq_pool*)r->data;
                    write_out(d);
                }
                thread_pool_delete_result(r,1);
            }
            while(block == -1);
        }
        thread_pool_process_flush(q);
        while (( r = thread_pool_next_result(q))) {
            struct bseq_pool *d = (struct bseq_pool*)r->data;
            write_out(d);
            thread_pool_delete_result(r,1);
        }
        thread_pool_process_destroy(q);
        thread_pool_destroy(p);
    }
    memory_release();
    return 0;
}

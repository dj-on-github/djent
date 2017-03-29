
/*
    djent - A reimplementation of Fourmilab's ent with several improvements. 
    
    Copyright (C) 2017  David Johnston

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    -----
    
    Contact. David Johnston dj@deadhat.com
*/

/* 0 for no messages. */
#define DEBUG 10
 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#ifdef _WIN32
#include "vsdjent/stdafx.h"
#include "ya_getopt/ya_getopt.h"
#else
#include <unistd.h> 
#include <getopt.h>
#define  errno_t int
#endif


#define QUEUESIZE 4096
#define BUFFSIZE  2048

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

unsigned char buffer[BUFFSIZE];
unsigned char queue[QUEUESIZE];
unsigned int queue_start;     /* FIFO pointers */
unsigned int queue_end;
unsigned int queue_size;

unsigned int current_byte;
unsigned int bits_used_from_byte;
unsigned int got_byte;
int64_t      current_symbol;
unsigned int bits_in_current_symbol;
int outcount;

uint64_t symbol_count;
uint64_t mean_total;

int terse;
uint64_t filebytes;

int opt;
int i;
unsigned int symbol_length;
int textmode;
int print_occurence;
int fold;

int use_stdin;
char *filename;
FILE *fp;
int terse_index;
int not_eof;
int escape;
int64_t symbol;

double ent;

uint64_t occurance_size;
uint64_t *occurance_count;
uint64_t occurance_total;

double chisq;
double chisq_sum;
double *chisq_prob;

uint64_t mp;
uint64_t monty_total_count;
uint64_t monty_inside_count;

double radiussquared;
double position_x;
double position_y;
double montepi;
uint64_t monte[6];


uint64_t count1;
uint64_t count0;

uint64_t count00;
uint64_t count01;
uint64_t count10;
uint64_t count11;

uint64_t symbol_mask;

uint64_t t1;
uint64_t t2;
uint64_t t3;
    
int scc_first;
uint64_t first_symbol;
uint64_t scc_previous;
uint64_t scc_count;
int scc_wrap;

double    result_mean;
uint64_t result_chisq_count;
double   result_chisq_distribution;
double   result_chisq_percent;
double  result_entropy;
double  result_pi;
double  result_pierr;
double  result_compression;
double  result_scc;


void update_monte_carlo(unsigned char symbol);

void display_usage() {
	fprintf(stderr, "Usage: djent [-b] [-l <n>] [-c] [-u] [-h] [-f] [-t] [filename] [filename2] ...\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Compute statistics of random data.\n");
	fprintf(stderr, "  Author: David Johnston, dj@deadhat.com\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "  -l <n>    --symbol_length=<n> Treat incoming data symbols as bitlength n. Default is 8.\n");
	fprintf(stderr, "  -b        --binary            Treat incoming data as binary. Equivalent to -l 1\n");
	fprintf(stderr, "  -c        --occurence         Print symbol occurence counts\n");
	fprintf(stderr, "  -w        --scc_wrap          Treat data as cyclical in SCC\n");
	fprintf(stderr, "  -f        --fold              Fold uppercase letters to lower case\n");
	fprintf(stderr, "  -t        --terse             Terse output\n");
	fprintf(stderr, "  -h or -u  --help              Print this text\n");
}

uint64_t ipow(uint64_t base, uint64_t exp)
{
	uint64_t result = 1;
	while (exp)
	{
		if (exp & 1)
			result *= base;
		exp >>= 1;
		base *= base;
	}

	return result;
}

/* Chi Square P value computation */

double igf(double s, double z)
{
    double sc;
    double sum = 1.0;
    double nom = 1.0;
    double denom = 1.0;
    int k;
    
    if(z < 0.0) return 0.0;
    
    sc = (1.0 / s);
    
    sc = sc * pow(z, s);
    sc = sc * exp(-z);
 
    for(k = 0; k < 200; k++) {
	    nom *= z;
	    s++;
	    denom *= s;
	    sum += (nom / denom);
    }
 
    return sum * sc;
}

double chisqr(double crit, int df)
{
    double k;
    double x;
    double p;
    
    if((crit < 0) || (df < 1)) return (0.0);
    
    k = df * 0.5;
    x = crit * 0.5;
    if(df == 2) return exp(-1.0 * x);
    
    p = igf(k, x);
    
    if(isnan(p) || isinf(p) || p <= 1e-8) return 1e-14;

	p = p / tgamma(k);
    return (1.0 - p);
}

/* The queue
 *
 * This implements a FIFO into which bytes are pushed from a file and 
 * from which symbols (of the chosen size) are pulled from the other end.
 * Data from the file is read into buffer and that data is used to fill the
 * input side of the queue. The queue is twice as big as the buffer so the 
 * buffer read is done when the queue is less than half full.
 * It treats bits within bytes as big endian (I.E. MSB arrived first from ES).
 * There will be an option to switch to little endian at some point.
 */
 
 
void init_byte_queue() {
    int i;
    
    /* printf("Init Byte Queue\n"); */
    queue_start = 0;
    queue_end = 0;
    queue_size = 0;
    
    got_byte = 0;
    current_byte = 0;
    bits_used_from_byte = 0;
    current_symbol = 0;
    bits_in_current_symbol = 0;
    
    for (i=0;i<QUEUESIZE;i++) queue[i]=0;
    
    symbol_mask = ipow((uint64_t)2,(uint64_t)symbol_length)-1;
}

int fill_byte_queue(FILE *fp) {
    int len;
    int space;
    int i;
    
    space = QUEUESIZE-queue_size; /* Dont pull more data than needed */
    if (space > BUFFSIZE) space = BUFFSIZE;
    
    /* ("  queue: space=%d\n",space); */
    len = fread(buffer, (size_t)1,(size_t)space, fp);
    if (len==0) {
        /* printf("  queue: len = %d\n",len); */
        return 0;
    }
    
    for (i=0;i<len;i++) {
        queue[(queue_end+i) % QUEUESIZE] = buffer[i];
		update_monte_carlo(buffer[i]);
    }
    
    filebytes += len;
    
    queue_size += len;
    queue_end = ((queue_end + len) % QUEUESIZE);
    
    /* printf("  Queue : size %d,  start %d,  end %d\n",queue_size, queue_start, queue_end); */
    return len;
}

/* pull symbol length bits off the start of the queue */
int64_t get_symbol(uint64_t symbol_length) {

    unsigned int temp;
    
    current_symbol = 0;
    
    /* Get a byte if we don't have one */
    if (got_byte == 0) {
        if (((queue_size*8) < symbol_length)) return -1; /* Uh oh. Empty */
        
        current_byte = queue[queue_start];
        queue_start = (queue_start+1) % QUEUESIZE;
        queue_size -= 1;
        bits_used_from_byte = 0;
        got_byte=1;
    }
    
    /* Move bits from current byte pulled from queue to current symbol */
    if (symbol_length == 1) { /* Optimize for the single bit size case */
        current_symbol = (current_byte & 0x80) >> 7;
        current_byte <<= 1;
        bits_used_from_byte++;
        if (bits_used_from_byte == 8) {
            got_byte = 0;
            bits_used_from_byte = 0;
        }
        return current_symbol;
    } else if (symbol_length == 8) { /* optimize for the byte size case */
        current_symbol = current_byte;
        got_byte = 0;
        return current_symbol;
    } else {  /* Symbol Length != 8 or 1, do it bit by bit */
              /* Later maybe optimize when > 7 bits needed */
        /* Take upper symbol_length bits */
        bits_in_current_symbol = 0;
        do {
            temp = (current_byte & 0x80) >> 7;
            current_byte = (current_byte << 1) & 0xff;
            bits_used_from_byte++;
            if (bits_used_from_byte == 8) {
                got_byte = 0;
                bits_used_from_byte = 0;
            }
            current_symbol = ((current_symbol << 1) | temp) & symbol_mask;
            bits_in_current_symbol++;
            
            /* fetch a new byte from queue if we aren't done yet */
            if (got_byte == 0) {
                if (((queue_size*8) < symbol_length)) return -1; /* Uh oh. Empty */
        
                current_byte = queue[queue_start];
                queue_start = (queue_start+1) % QUEUESIZE;
                queue_size -= 1;
                bits_used_from_byte = 0;
                got_byte=1;
            }
        } while (bits_in_current_symbol < symbol_length);
        return current_symbol;
    }
}


/* The initialize routines for the various metrics */

void init_mean() {
    outcount = 0;
    mean_total = 0;
};

void init_entropy() {
    ent = 0.0;
};

void init_occurances() {
    uint64_t i;
    
    occurance_total = 0;
    if (symbol_length > 32) {
        fprintf(stderr,"Error, symbol length cannot be longer than 32 bits for occurence count table\n");
        exit(1);
    }
    occurance_size = ipow(2,symbol_length);
    occurance_count = (uint64_t *) malloc (sizeof(uint64_t)*occurance_size);
    /* printf("mallocating %lld bytes\n", (sizeof(uint64_t)*occurance_size));
     */
    if (occurance_count == NULL) {
        #ifdef _WIN32
        fprintf(stderr,"Error, unable to allocate %lld bytes of memory for the occurence count\n",(sizeof(uint64_t)*occurance_size));
        #elif __llvm__
        fprintf(stderr,"Error, unable to allocate %lld bytes of memory for the occurence count\n",(sizeof(uint64_t)*occurance_size));
        #elif __linux__
        fprintf(stderr,"Error, unable to allocate %ld bytes of memory for the occurence count\n",(sizeof(uint64_t)*occurance_size));
        #endif
        exit(1);
    }

    for (i=0;i<occurance_size;i++) occurance_count[i] = 0;
};

void init_chisq() {
    chisq = 0.0;
    chisq_prob = (double *) malloc (sizeof(double)*occurance_size);
    /* printf("mallocating %lld bytes for chisq probability table\n", (sizeof(double)*occurance_size));
    */
    if (chisq_prob == NULL) {
        exit(1);
    }
    for (i=0;i<occurance_size;i++) chisq_prob[i] = 0.0;
};

void init_filesize() {
};

void init_monte_carlo() {
    mp = 0;
    monty_total_count = 0;
    monty_inside_count = 0;
	radiussquared = (256.0 * 256.0 * 256.0) - 1;
	radiussquared = radiussquared*radiussquared;

};

void init_compression() {
    /* nothing to do here */
};

void init_scc() {
    t1 = 0;
    t2 = 0;
    t3 = 0;
    
    scc_first = 1;
    scc_previous = 0;
    scc_count = 0;
    first_symbol = 0;
    
    count00=0;
    count01=0;
    count10=0;
    count11=0;
};           
                
/* The update routines for the various metrics */
        
void update_mean(uint64_t symbol) {
    mean_total += symbol;
    /* printf(" mean_total = %lld,  count=%lld\n",mean_total,symbol_count); */
};

void update_entropy(uint64_t symbol) {
	/* nothin to do here */
};

void update_occurances(uint64_t symbol) {
    occurance_count[symbol]++;
    occurance_total++;
};

void update_chisq(uint64_t symbol) {
	/* Nothing to do here */
};

void update_filesize(uint64_t symbol) {
	/* Nothing to do here */
};

void update_monte_carlo(unsigned char symbol) {
	int mj;

	monte[mp++] = symbol;

	if (mp > 5) {
		mp = 0;
		monty_total_count++;
		position_x = 0;
		position_y = 0;
		for (mj = 0; mj < 3; mj++) {
			position_x = (position_x * 256.0) + monte[mj];
			position_y = (position_y * 256.0) + monte[3 + mj];
		}
		if (((position_x * position_x) + (position_y *  position_y)) <= radiussquared) {
			monty_inside_count++;
		}
	}
};

void update_compression(uint64_t symbol) {
 /* nothing to do here */
};

void update_scc(uint64_t symbol) {
    /* We need two symbols to start, so skip the first symbol */
    scc_count++;
    
    if (scc_first==1) {
        scc_first = 0;
        first_symbol = symbol;
    } else {
        t1 += (scc_previous * symbol);
    }
    t2 += symbol*symbol;
    t3 += symbol;
    
    /* printf("symbol %02X, count=%llu t1= %llX,  t2= %llx, t3= %llx\n",symbol,scc_count,t1,t2,t3); */
    scc_previous = symbol;
};

/* The finalization routines for the various metrics */
        
void finalize_mean() {
    double mean;
    mean = (double)mean_total/(double)symbol_count; 

	result_mean = mean;
	return;

    if (terse==1) printf("%f,",mean);
    else printf("   Mean = %f\n",mean);
};

void finalize_entropy() {
	unsigned int eloop;
	ent = 0.0;
	for (eloop = 0; eloop < occurance_size; eloop++) {
		if (chisq_prob[eloop] > 0.0) {
			ent += (chisq_prob[eloop] * log10(1.0 / chisq_prob[eloop]) *  3.32192809488736234787);
		}
	}

	result_entropy = ent;
	return;

	if (terse == 1) printf("%f,", ent);
	else printf("   Shannon Entropy = %f\n", ent);
};

void finalize_occurances() {
};

void finalize_chisq() {
    uint64_t i;
    double diff;
    double chisq_final_prob;
    
    double expected;
    expected = (double)occurance_total / (double)occurance_size;
    for (i=0; i < occurance_size; i++) {
        diff = (double)(occurance_count[i]) - expected;
        chisq_prob[i] = ((double)occurance_count[i])/occurance_total;
        chisq      += (diff*diff)/expected;
        chisq_sum  += (double)(i * occurance_count[i]);
    }
    
    /*chisq_final_prob = pochisq(chisq, (occurance_size-1)); */
    chisq_final_prob = chisqr(chisq, (occurance_size-1));
    
	result_chisq_count = occurance_total;
	result_chisq_distribution = chisq;
	result_chisq_percent = chisq_final_prob * 100;

	return;

    if (terse==1) {
        printf(" ,%1.2f ",chisq_final_prob * 100); 
    } else {
        #ifdef _WIN32
        printf("   Chi square: symbol count=%llu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", occurance_total, chisq, chisq_final_prob * 100);
        #elif __llvm__
        printf("   Chi square: symbol count=%llu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", occurance_total, chisq, chisq_final_prob * 100);
        #elif __linux__
        printf("   Chi square: symbol count=%lu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", occurance_total, chisq, chisq_final_prob * 100);
        #endif
    }

};

void finalize_filesize() {
};

void finalize_monte_carlo() {
	double pierr;
	double montepi;

	montepi = 4.0 * (((double)monty_inside_count) / monty_total_count);

	pierr = (fabs(M_PI - montepi) / M_PI)*100.0;

	result_pi = montepi;
	result_pierr = pierr;

	return;

	if (terse == 1) {
		printf(" ,%f ", montepi);
	}
	else {
		printf("   Monte Carlo value for Pi is %f (error %1.2f percent).\n", montepi, pierr);
	}

};

void finalize_compression() {
	double compression;

	compression = (100.0 * (symbol_length - ent)) / symbol_length;

	result_compression = compression;

	return;

	if (terse == 1) {
		printf(", %f ", compression);
	}
	else {
		printf("   Optimal compression would compress by %f percent\n", compression);
	}
};

void finalize_scc() {
    double scc;
    int64_t top;
    int64_t bottom;
    
    if (scc_wrap==1) {
            t1 += (scc_previous * first_symbol);
            t2 += first_symbol*first_symbol;
            t3 += first_symbol;
    } else {
        scc_count -= 1;
    }
    
    /* need signed arithmetic because we are subtracting */
    top = (int64_t)(scc_count * t1) - (int64_t)(t3*t3);
    bottom = (int64_t)(scc_count * t2) - (int64_t)(t3*t3);
    scc = (double)top/(double)bottom;


	result_scc = scc;
	return;

    if (terse==1) {
        printf("%08f",scc);
    } else {
        printf("   Serial Correlation = %f\n",scc);
    }
};

/********
* main() is mostly about parsing and qualifying the command line options.
*/

int main(int argc, char** argv)
{
    /* Defaults */
    symbol_length = 8;
    textmode = 1;
    print_occurence = 0;
    fold = 0;
    terse = 0;
    use_stdin = 1;
    fp = NULL;
    terse_index = 0;
    escape = 0;
    scc_wrap = 0;

	#define ERRSTRINGSIZE 256
    #ifdef _WIN32
	errno_t err;
	char errstring[ERRSTRINGSIZE];
    #endif
    int filenumber = 0;
    
    char optString[] = "bcwfthul:";
    int longIndex;
    static const struct option longOpts[] = {
    { "symbol_length", required_argument, NULL, 'l' },
    { "binary", no_argument, NULL, 'b' },
    { "occurence", no_argument, NULL, 'c' },
    { "fold", no_argument, NULL, 'f' },
    { "scc_wrap", no_argument, NULL, 'w' },
    { "terse", no_argument, NULL, 't' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
    };

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'b':
                symbol_length = 1; /* binary mode treats newlines as data */
                textmode = 0;
                break;
                
            case 'l':
                symbol_length = atoi(optarg); /* -b -l <n> is valid. symbol length
                                               * will be <n> and newlines will be treated
                                               * as data.
                                               */
                break;
                
            case 'c':
                print_occurence = 1;
                break;
            
            case 'w':
                scc_wrap = 1;
                break;
                
            case 'f':
                fold = 1;
                break;
            
            case 't':
                terse = 1;
                break;
            
            case 'u':    
            case 'h':   /* fall-through is intentional */
            case '?':
                display_usage();
                exit(0);
                 
            default:
                /* You won't actually get here. */
                break;
        }
         
        opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    } // end while
    

	/* Range check the var args */

	escape = 0;
	if (symbol_length < 1) {
		printf("Error: Symbol length %d must not be 0 or negative. \n",symbol_length);
		escape = 1;
	}

	if (escape==1) {
		exit(1);
	}

    init_byte_queue();
    
    /* Loop through the filenames */
	if (optind==argc) {
		use_stdin = 1;
	}
	else {
		use_stdin = 0;
	}

	terse_index = 0;
	/* if (use_stdin == 0) { */
	filenumber = optind;

	do {
		terse_index++;
		filebytes = 0;
        /* printf("OPTIND %d, filenumber %d, ARGC %d\n",optind,filenumber,argc); */
		if (use_stdin==1) {
			use_stdin = 1;
			if (textmode != 1) freopen(NULL, "rb", stdin);
			fp = stdin;
		}
		else {
			filename = argv[filenumber];

			if (textmode == 1) {
				if (terse == 0) printf(" opening %s as hex text\n", filename);
                #ifdef _WIN32
				if ((err = fopen_s(&fp, filename, "r")) != 0) {
					strerror_s(errstring, ERRSTRINGSIZE, err);
					fprintf(stderr, "Error : Unable to open file %s, %s\n", filename, errstring);
					exit(1);
				}
                #else
				fp = fopen(filename,"r");
                if (fp == NULL) {
					fprintf(stderr, "Error : Unable to open file %s\n", filename);
					exit(1);
                }
                #endif
			}
			else {
				if (terse == 0) printf(" opening %s as binary\n", filename);
                #ifdef _WIN32
				if ((err = fopen_s(&fp, filename, "rb")) != 0) {
					strerror_s(errstring, ERRSTRINGSIZE, err);
					fprintf(stderr, "Error : Unable to open file %s, %s\n", filename, errstring);
					exit(1);
				}
                #else
				fp = fopen(filename,"rb");
                if (fp == NULL) {
					fprintf(stderr, "Error : Unable to open file %s\n", filename);
					exit(1);
                }
                #endif
				/*  fp = fopen(filename,"rb"); */
				/* printf("           %x\n",(unsigned int)fp);*/
			}

			if (terse == 0) printf(" Symbol Size(bits) = %d\n", symbol_length);

		}


		/* Print terse header if necessary */
		if ((terse == 1) && (terse_index == 1)) {
			printf("0,  File-bytes,    Entropy,     Chi-square,  Mean,        Monte-Carlo-Pi, Serial-Correlation, Filename\n");
		}

		/* Initialize the metrics */
		symbol_count = 0;

		init_mean();
		init_entropy();
		init_occurances();
		init_chisq();
		init_filesize();
		init_monte_carlo();
		init_compression();
		init_scc();

		/* Now process the file fp */
		/* Since we have multiple possible symbol sizes, first pull the bytes
			* then put them in a queue which behaves like a bitwise queue, then
			* pull the symbols from the queue.
			*/
		do {
			if (queue_size == 0) {
				not_eof = fill_byte_queue(fp); /* get bytes from file into queue */
				if (not_eof == 0) break;
			}
			symbol = get_symbol(symbol_length);      /* Pull a symbol from the queue */
			symbol_count++;

			/* Finish up if no symbols left in queue */
			if (symbol == -1) break;

			/* Then update the algorithms using the symbol */
			update_mean(symbol);
			update_entropy(symbol);
			update_occurances(symbol);
			update_chisq(symbol);
			update_filesize(symbol);
			/* Monte Carlo is different, it works on bytes, not symbols
				* So we call the update from within the fill_byte_queue routine
				*/
				/* update_monte_carlo(symbol); */

			update_compression(symbol);
			update_scc(symbol);
		} while (1 == 1);

		finalize_mean();
		finalize_occurances();
		finalize_chisq();
		finalize_entropy();
		finalize_filesize();
		finalize_monte_carlo();
		finalize_compression();
		finalize_scc();

		if (terse == 1) {
            #ifdef _WIN32
			printf("%d,%12ld,%12f,%12f,%12f,%12f,   %12f,           %s\n", terse_index, filebytes, result_entropy, result_chisq_percent, result_mean, result_pi, result_scc, filename);
            #elif __llvm__
			printf("%d,%12lld,%12f,%12f,%12f,%12f,   %12f,           %s\n", terse_index, filebytes, result_entropy, result_chisq_percent, result_mean, result_pi, result_scc, filename);
            #elif __linux__
			printf("%d,%12ld,%12f,%12f,%12f,%12f,   %12f,           %s\n", terse_index, filebytes, result_entropy, result_chisq_percent, result_mean, result_pi, result_scc, filename);
            #endif
		}
		else {
            printf("   Shannon IID Entropy = %f bits per symbol\n",result_entropy);
		    printf("   Optimal compression would compress by %f percent\n", result_compression);
            #ifdef _WIN32
            printf("   Chi square: symbol count=%llu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", result_chisq_count, result_chisq_distribution, result_chisq_percent);
            #elif __llvm__
            printf("   Chi square: symbol count=%llu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", result_chisq_count, result_chisq_distribution, result_chisq_percent);
            #elif __linux__
            printf("   Chi square: symbol count=%lu, distribution=%1.2f, randomly exceeds %1.2f percent of the time\n", result_chisq_count, result_chisq_distribution, result_chisq_percent);
            #endif
            printf("   Mean = %f\n",result_mean);
		    printf("   Monte Carlo value for Pi is %f (error %1.2f percent).\n", result_pi, result_pierr);
            printf("   Serial Correlation = %f\n",result_scc);
		}

		filenumber++;
	} while ((filenumber < argc) && (use_stdin != 1));

	/* Find out what the various compilers give us
    #ifdef __llvm__
        printf("llvm\n");
    #endif

    #ifdef __clang__
        printf("clang\n");
    #endif
    
    #ifdef __gcc__
        printf("gcc\n");
    #endif
    #ifdef __linux__
        printf("linux\n");
    #endif
	#ifdef _WIN32
		printf("win32\n");
	#endif
    */
	return 0;

}



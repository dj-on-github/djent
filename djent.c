
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

#define __STDC_FORMAT_MACROS
#include <inttypes.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
/* #include <regex.h>*/

#ifndef NO_GMP
#include <mpfr.h>
#endif

#ifdef _WIN32
/* #include "vsdjent/stdafx.h" */
#include "ya_getopt/ya_getopt.h" /* NOTE: VS2015 goes not have getopt. Put ya_getopt in the directory. From here https://github.com/kubo/ya_getopt */
#else
#include <unistd.h> 
#include <getopt.h>
#define  errno_t int
#endif

#define MAX_ERROR_MSG 0x1000
#define QUEUESIZE 4096
#define BUFFSIZE  2048

#ifndef M_PI
#define M_PI 3.1415926535897932384626
#endif

double pochisq(double ax,int df);

unsigned char buffer[BUFFSIZE];
unsigned char buffer2[BUFFSIZE+4];
unsigned char queue[QUEUESIZE];
unsigned int queue_start;     /* FIFO pointers */
unsigned int queue_end;
size_t queue_size;

unsigned int buffer2_size;

unsigned int current_byte;
unsigned int bits_used_from_byte;
unsigned int got_byte;
int64_t      current_symbol;
unsigned int bits_in_current_symbol;
int outcount;
uint64_t scc_fifo[256];
uint64_t scc_first_lagn[256];

uint64_t symbol_count;
uint64_t mean_total;

int terse;
int ent_exact;
int suppress_header;
uint64_t filebytes;

double voltage;
double temperature;
unsigned char deviceid[256];
unsigned char process[256];
unsigned char processing[256];

int opt;
unsigned int symbol_length;
int hexmode;
int print_occurrence;
int print_longest;
int fold;
int lagn;
int byte_reverse;
int word_reverse;
int parse_filename;

int use_stdin;
char *filename;
FILE *fp;
int terse_index;
int not_eof;
int64_t symbol;

char inputlistfilename[256];
int using_inputlistfile;

double ent;

uint64_t occurrence_size;
uint64_t *occurrence_count;
uint64_t occurrence_total;
int no_occurrence_space;

uint64_t longest_size;
uint64_t *longest_count;
uint64_t longest_total;
int no_longest_space;

uint64_t longest_last_symbol;
uint64_t longest_run;
uint64_t longest_longest;
uint64_t longest_longest_symbol;

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

uint64_t aeqb_count;
uint64_t mean_count;
double   other_scc;

double    result_mean;
uint64_t result_chisq_count;
double   result_chisq_distribution;
double   result_chisq_percent;
double  result_entropy;
double  result_min_entropy;
uint32_t result_min_entropy_symbol;
double  result_pi;
double  result_pierr;
double  result_compression;
double  result_scc;
double  result_longest_pvalue;

const unsigned char byte_reverse_table[] = {
  0x00,0x80,0x40,0xC0,0x20,0xA0,0x60,0xE0,0x10,0x90,0x50,0xD0,0x30,0xB0,0x70,0xF0, 
  0x08,0x88,0x48,0xC8,0x28,0xA8,0x68,0xE8,0x18,0x98,0x58,0xD8,0x38,0xB8,0x78,0xF8, 
  0x04,0x84,0x44,0xC4,0x24,0xA4,0x64,0xE4,0x14,0x94,0x54,0xD4,0x34,0xB4,0x74,0xF4, 
  0x0C,0x8C,0x4C,0xCC,0x2C,0xAC,0x6C,0xEC,0x1C,0x9C,0x5C,0xDC,0x3C,0xBC,0x7C,0xFC, 
  0x02,0x82,0x42,0xC2,0x22,0xA2,0x62,0xE2,0x12,0x92,0x52,0xD2,0x32,0xB2,0x72,0xF2, 
  0x0A,0x8A,0x4A,0xCA,0x2A,0xAA,0x6A,0xEA,0x1A,0x9A,0x5A,0xDA,0x3A,0xBA,0x7A,0xFA,
  0x06,0x86,0x46,0xC6,0x26,0xA6,0x66,0xE6,0x16,0x96,0x56,0xD6,0x36,0xB6,0x76,0xF6, 
  0x0E,0x8E,0x4E,0xCE,0x2E,0xAE,0x6E,0xEE,0x1E,0x9E,0x5E,0xDE,0x3E,0xBE,0x7E,0xFE,
  0x01,0x81,0x41,0xC1,0x21,0xA1,0x61,0xE1,0x11,0x91,0x51,0xD1,0x31,0xB1,0x71,0xF1,
  0x09,0x89,0x49,0xC9,0x29,0xA9,0x69,0xE9,0x19,0x99,0x59,0xD9,0x39,0xB9,0x79,0xF9, 
  0x05,0x85,0x45,0xC5,0x25,0xA5,0x65,0xE5,0x15,0x95,0x55,0xD5,0x35,0xB5,0x75,0xF5,
  0x0D,0x8D,0x4D,0xCD,0x2D,0xAD,0x6D,0xED,0x1D,0x9D,0x5D,0xDD,0x3D,0xBD,0x7D,0xFD,
  0x03,0x83,0x43,0xC3,0x23,0xA3,0x63,0xE3,0x13,0x93,0x53,0xD3,0x33,0xB3,0x73,0xF3, 
  0x0B,0x8B,0x4B,0xCB,0x2B,0xAB,0x6B,0xEB,0x1B,0x9B,0x5B,0xDB,0x3B,0xBB,0x7B,0xFB,
  0x07,0x87,0x47,0xC7,0x27,0xA7,0x67,0xE7,0x17,0x97,0x57,0xD7,0x37,0xB7,0x77,0xF7, 
  0x0F,0x8F,0x4F,0xCF,0x2F,0xAF,0x6F,0xEF,0x1F,0x9F,0x5F,0xDF,0x3F,0xBF,0x7F,0xFF
};


void update_monte_carlo(unsigned char symbol);

void display_usage() {
	fprintf(stderr, "Usage: djent [-brRpcuhds] [-l <n>] [-i <input file list filename>] [filename] [filename2] ...\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Compute statistics of random data.\n");
	fprintf(stderr, "  Author: David Johnston, dj@deadhat.com\n");
	fprintf(stderr, "\n");

	fprintf(stderr, "  -i <filename>  --inputfilelist=<filename> Read list of filenames from <filename>\n");
	fprintf(stderr, "  -p             --parse_filename           Extract CID, Process, Voltage and Temperature from filename.\n");
	fprintf(stderr, "                                            The values will be included in the output.\n");
	fprintf(stderr, "  -l <n>         --symbol_length=<n>        Treat incoming data symbols as bitlength n. Default is 8.\n");
	fprintf(stderr, "  -b             --binary                   Treat incoming data as binary. Default bit length will be -l 1\n");
	fprintf(stderr, "  -r             --byte_reverse             Reverse the bit order in incoming bytes\n");
	fprintf(stderr, "  -R             --word_reverse             Reverse the byte order in incoming 4 byte words\n");
	fprintf(stderr, "  -c             --occurrence               Print symbol occurrence counts\n");
	fprintf(stderr, "  -C             --longest                  Print symbol longest run counts\n");
	fprintf(stderr, "  -w             --scc_wrap                 Treat data as cyclical in SCC\n");
	fprintf(stderr, "  -n <n>         --lagn=<n>                 Lag gap in SCC. Default=1\n");
	fprintf(stderr, "  -f             --fold                     Fold uppercase letters to lower case\n");
	fprintf(stderr, "  -t             --terse                    Terse output\n");
	fprintf(stderr, "  -e             --ent_exact                Exactly match output format of ent\n");
	fprintf(stderr, "  -s             --suppress_header          Suppress the header in terse output\n");
	fprintf(stderr, "  -h or -u       --help                     Print this text\n");

    fprintf(stderr, "\n Notes\n");
    fprintf(stderr,   "   * By default djent is in hex mode where it reads ascii hex data and converts it to binary to analyze.\n");
    fprintf(stderr,   "     In hex mode, the symbol length defaults to 8, so normal hex files can be treated as a representation\n");
    fprintf(stderr,   "     of bytes. The symbol length can be changed to any value between 1 and 32 bits using the -l <n> option.\n");
    fprintf(stderr,   "   * With the -b option djent switches to binary reads in each byte as binary with a symbol length of 1.\n");
    fprintf(stderr,   "   * To analyze ascii text instead of hex ascii, you need djent to treat each byte as a separate symbol, so\n");
    fprintf(stderr,   "     use binary mode with a symbol length of 8. I.E. djent -b -l 8 <filename>\n");
    fprintf(stderr,   "   * By default djent treats the MSB of each byte as the first. This can be switched so that djent treats\n");
    fprintf(stderr,   "     the LSB as the first bit in each byte using the -r option.\n");
    fprintf(stderr,   "   * Terse output is requested using -t. This outputs in CSV format. The first line is the header. If\n");
    fprintf(stderr,   "     multiple files are provided, there will be one line of CSV output per file in addition to the header.\n");
    fprintf(stderr,   "     The CSV header can be suppressed with -s.\n");
    fprintf(stderr,   "   * To analyze multiple files, just give multiple file names on the command line. To read data in from\n");
    fprintf(stderr,   "     the command line, don't provide a filename and pipe the data in. <datasource> | djent\n");
    fprintf(stderr,   "   * The parse filename option =p picks takes four patterns from the filename to include in the output,\n");
    fprintf(stderr,   "     This is so that it is easy to plot test conditions that are commonly encoded in a filename.\n");
    fprintf(stderr,   "     Fields are delimited by uderscores. The four patters for CID, process, Voltage and Temperature are:\n");
    fprintf(stderr,   "     _CID-<componentID>_ , _PROC-<process info>_, _<x>p<y>V_ and _<x>p<y>C_ . 'p' is the decimal point.\n");
    fprintf(stderr,   "   * To compute the statistics, djent builds a frequency table of the symbols. This can be displayed\n");
    fprintf(stderr,   "     using the -c option. The size of this table is what limits the the maximum symbol size. For each\n");
    fprintf(stderr,   "     of the 2^n symbols, a 64 bit entry in a table is created. So for n=32, that's 32GBytes so the ability\n");
    fprintf(stderr,   "     to handle large symbol sizes is limited by the available memory and the per process allocation limit.\n");
    fprintf(stderr,   "   * The serial correlation coefficient is not wrap around by default, meaning that it does not compare\n");
    fprintf(stderr,   "     the last value in the data with the first. To get wrap around behaviour, use the -w option.\n");
    fprintf(stderr,   "   * The Lag-N correlation coefficient can be computed by using the -n <n> option. This causes the SCC\n");
    fprintf(stderr,   "     computation to compare each Xth symbol with the (X+n)th symbol instead of the (X+1)th symbol.\n");
    fprintf(stderr,   "     If you use wrap around with Lag-N, then the wrap around will reach n bits further into the start\n");
    fprintf(stderr,   "     of the sequence.\n");
    fprintf(stderr,   "   * The byte reverse option -r reverses the order of bits within each byte. The word reverse option -R\n");    
    fprintf(stderr,   "     reverses the order of bytes within each 32 bit word, from 3,2,1,0 to 0,1,2,3. Both -R and -r can\n");    
    fprintf(stderr,   "     be used together. Using -R with a data that isn't a multiple of 32 bits long will get padded with\n");    
    fprintf(stderr,   "     zeros, which may not be what you want. A padding warning will be sent to STDERR.\n");    
    fprintf(stderr,   "   * Instead of providing data file names on the command line, djent can be told to read a list of files\n");
    fprintf(stderr,   "     from a text file. The file must have one filename per line. Lines beginning with # will be ignored.\n");
    fprintf(stderr,   "     Use the -i <filename> option to request that djent reads the file list from <filename>.\n");

    fprintf(stderr, "\n Examples\n");
    fprintf(stderr,   "   Print this help\n");
    fprintf(stderr,   "     djent -h\n\n");
    fprintf(stderr,   "   Analyze hex file from stdin\n");
    fprintf(stderr,   "     cat datafile.hex | djent\n\n");
    fprintf(stderr,   "   Analyze binary file\n");
    fprintf(stderr,   "     djent -b datafile.bin\n\n");
    fprintf(stderr,   "   Analyze several files with CSV output\n");
    fprintf(stderr,   "     djent -t data1.hex data2.hex data3.hex\n\n");
    fprintf(stderr,   "   Analyze ascii symbols - Read in binary and set symbol size to 8.\n");
    fprintf(stderr,   "     djent -b -l 8  textfile.txt\n\n");
    fprintf(stderr,   "   Analyze binary file with parsable filename.\n");
    fprintf(stderr,   "     djent -b -t -p  rawdata_CID-X23_PROC-TTFT_1p2V_25p0C_.bin\n");

}


// Return the probability of the longest run of heads being less than or equal to n
// in a sequence of r uniform coin tosses. Use MPFR to avoid overflows.
double longest_run_cdf(unsigned int ui_n,unsigned int ui_r) { // n=longest run. r = length of data sequence
    double answer;
    mpfr_set_default_prec(1024);
	mpfr_t n;
	mpfr_t r;
    mpfr_t topa;
    mpfr_t bottoma;
    mpfr_t first;
    mpfr_t topb;
    mpfr_t nplusone;
    mpfr_t bottomb;
    mpfr_t nplus2over2;
    mpfr_t second;
    mpfr_t mpfans;
    mpfr_set_default_prec(1024);

    mpfr_init(topa);
    mpfr_init(nplusone);
    mpfr_init(bottoma);
    mpfr_init(first);
    mpfr_init(topb);
    mpfr_init(bottomb);
    mpfr_init(nplus2over2);
    mpfr_init(second);
    mpfr_init(mpfans);
    mpfr_init_set_ui(n,ui_n,MPFR_RNDN);
    mpfr_init_set_ui(r,ui_r,MPFR_RNDN);

    mpfr_add_ui(topa,r,1,MPFR_RNDN);

    mpfr_add_ui(nplusone,n,1,MPFR_RNDN);


    mpfr_exp2(bottoma,nplusone,MPFR_RNDN);
    mpfr_sub(bottoma,bottoma,n,MPFR_RNDN);
    mpfr_sub_ui(bottoma,bottoma,2,MPFR_RNDN);

    mpfr_div(first,topa,bottoma,MPFR_RNDN);
    mpfr_neg(first,first,MPFR_RNDN);

    mpfr_exp(first,first,MPFR_RNDN);

    // Second
    mpfr_exp2(topb,nplusone,MPFR_RNDN);
    mpfr_sub_ui(topb,topb,1,MPFR_RNDN);

    mpfr_exp2(bottomb,nplusone,MPFR_RNDN);

    mpfr_add_ui(nplus2over2,n,2,MPFR_RNDN);
    mpfr_div_ui(nplus2over2,nplus2over2,2,MPFR_RNDN);

    mpfr_sub(bottomb,bottomb,nplus2over2,MPFR_RNDN);

    mpfr_div(second,topb,bottomb,MPFR_RNDN);

    //Final
    mpfr_mul(mpfans,first,second,MPFR_RNDN);
    answer = mpfr_get_d(mpfans,MPFR_RNDN);

    mpfr_clear(topa);
    mpfr_clear(nplusone);
    mpfr_clear(bottoma);
    mpfr_clear(first);
    mpfr_clear(topb);
    mpfr_clear(bottomb);
    mpfr_clear(nplus2over2);
    mpfr_clear(second);
    mpfr_clear(mpfans);
    mpfr_clear(n);
    mpfr_clear(r);

    return answer;


}

int count_lines_in_file(char *filename) {
     FILE *fp = fopen(filename,"r");
     int ch=0;
     int lines=0;

     if (fp == NULL) return 0;
     lines++;
     while ((ch = fgetc(fp)) != EOF) {
         if ((char)ch == '\n') lines++;
     }
     fclose(fp);
     return lines;
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

double zcdf(double z) {
    double w;
    double x;
    double y;
    double result;

    if (z == 0.0) return 0.5;

    y = fabs(z)/2.0;

    if (y >= 3.0) return 0.0;
    
    if (y < 1.0) {
        w = y * y;
        x =         0.000124818987;
        x = x * w - 0.001075204047;
        x = x * w + 0.005198775019;
        x = x * w - 0.019198292004;
        x = x * w + 0.059054035642;
        x = x * w - 0.151968751364;
        x = x * w + 0.319152932694;
        x = x * w - 0.531923007300;
        x = x * w + 0.797884560593;
        x = x * 2.0 * y;
    } else {
        y -= 2.0;
        x =        -0.000045255659;
        x = x * y + 0.000152529290;
        x = x * y - 0.000019538132;
        x = x * y - 0.000676904986;
        x = x * y + 0.001390604284;
        x = x * y - 0.000794620820;
        x = x * y - 0.002034254874;
        x = x * y + 0.006549791214;
        x = x * y - 0.010557625006;
        x = x * y + 0.011630447319;
        x = x * y - 0.009279453341;
        x = x * y + 0.005353579108;
        x = x * y - 0.002141268741;
        x = x * y + 0.000535310849;
        x = x * y + 0.999936657524;
    }


    if (z > 0.0) {
        result = (x/2.0)+0.5;
    } else {
        result = (0.5 - (x/2.0));
    }

    return result;
}

#define LOG_SQRT_PI 0.5723649429247000870717135 /* log (sqrt (pi)) */
#define I_SQRT_PI   0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define BIGX        20.0         /* max value to represent exp (x) */
#define ex(x)       (((x) < -BIGX) ? 0.0 : exp(x))

double chisqp(double ax, size_t df) {
    double x;
    double a;
    double y;
    double s;
    double e;
    double c;
    double z;
    int dfeven;
    
    dfeven=0;
    if ((df % 2)==0) dfeven = 1;

    x = ax;

    if (x <= 0.0 || df < 1) return 1.0;

    a = x/2.0;

    if (df > 1)  y = ex(-a);

    if (dfeven == 1) s = y;
    else s = 2.0 * zcdf(-sqrt(x));

    if (df > 2) {
        x = (df - 1.0)/2.0;
        if (dfeven==1) z = 1.0;
        else z = 0.5;

        if (a > BIGX) {
            if (dfeven==1) e = 0.0;
            else e = LOG_SQRT_PI;
            
            c = log(a);
            
            while (z <= x) {
                e = log(z) + e;
                s += ex(c * z - a - e);
                z += 1.0;
            }
            return (s);
        } else {
        if (dfeven==1) e = 1.0;
        else e = (I_SQRT_PI / sqrt(a));
        c = 0.0;
        while (z <= x) {
            e = e * (a / z);
            c = c + e;
            z += 1.0;
            }
        return (c * y + s);
        }
    } else {
        return s;
    }
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

int ishex(unsigned char c) {
    int result;
    result = 0;
    if      (c == '0') result = 1;
    else if (c == '1') result = 1;
    else if (c == '2') result = 1;
    else if (c == '3') result = 1;
    else if (c == '4') result = 1;
    else if (c == '5') result = 1;
    else if (c == '6') result = 1;
    else if (c == '7') result = 1;
    else if (c == '8') result = 1;
    else if (c == '9') result = 1;
    else if (c == 'a') result = 1;
    else if (c == 'b') result = 1;
    else if (c == 'c') result = 1;
    else if (c == 'd') result = 1;
    else if (c == 'e') result = 1;
    else if (c == 'f') result = 1;
    else if (c == 'A') result = 1;
    else if (c == 'B') result = 1;
    else if (c == 'C') result = 1;
    else if (c == 'D') result = 1;
    else if (c == 'E') result = 1;
    else if (c == 'F') result = 1;

    return result;
}

int ishexorx(unsigned char c) {
    int result;
    result = 0;
    if (ishex(c) == 1) result = 1;
    else if (c == 'x') result = 1;

    return result;
}

unsigned char hexpair[2];

int hexstate;

void init_hex2bin() {
    hexpair[0]=0x00;
    hexpair[1]=0x00;
    
    hexstate = 0;
}

/* 
 * This converts input hex text to binary. It uses a little
 * state machine to pull in 2 characters then convert them
 * to a byte. The state machine state is maintained across
 * calls so we do not lose values at read buffer boundaries.
 *
 * In the processing of the second character, an x will be
 * accepted if the first character is 0, so we get '0x'.
 * This allows the 0x prefixes to be eliminated without
 * accidentally treating the 0 as part of the data.
 */
 
size_t hex2bin(unsigned char *buffer, size_t len) {
    int outpos = 0;
    int scanpos = 0;
    unsigned char c;
    int byte;
    int nybble;
    /* Fetch characters until we get a hex 1.
     * Shift it into hexpair
     * If we have a valid hex pair put it in the buffer as binary
     * If we have 0x, drop it, including the 0.
     * If we have non hex, drop it.
     */

    do {
        if (hexstate == 0) {
            c = buffer[scanpos];
            if (ishex(c) == 1) {
                hexpair[0] = c;
                hexstate = 1;
            }
            scanpos++;
        }
        else if (hexstate == 1) {
            c = buffer[scanpos];
            if (((ishexorx(c) == 1) && (hexpair[0]=='0')) || (ishex(c)==1)){
                hexpair[1] = c;
                hexstate = 2;
            }
            scanpos++;
        }
        else if (hexstate == 2) {
            if ((hexpair[0]=='0') && (hexpair[1]=='x')) {
                hexstate = 0;
            }
            else { /* we have a valid hex pair */
                nybble = 0;                
                if ((((int)hexpair[0])>47) && (((int)hexpair[0])<58)){ /* 0-9 */
                    nybble = (int)hexpair[0] - 48;
                }
                else if ((((int)hexpair[0])>64) && (((int)hexpair[0])<71)){ /* A-F */
                    nybble = (int)hexpair[0] - 55;
                }
                else if ((((int)hexpair[0])>96) && (((int)hexpair[0])<103)){ /* a-f */
                    nybble = (int)hexpair[0] - 87;
                }

                nybble = nybble << 4;

                if ((((int)hexpair[1])>47) && (((int)hexpair[1])<58)){ /* 0-9 */
                    byte = nybble + (int)hexpair[1] - 48;
                }
                else if ((((int)hexpair[1])>64) && (((int)hexpair[1])<71)){ /* A-F */
                    byte = nybble + (int)hexpair[1] - 55;
                }
                else if ((((int)hexpair[1])>96) && (((int)hexpair[1])<103)){ /* a-f */
                    byte = nybble + (int)hexpair[1] - 87;
                }

                buffer[outpos++] = (unsigned char)byte;
                hexstate = 0;
            }
        }              
    } while (scanpos <= len);

    return outpos; /* return the number of bytes converted */
}

size_t fill_byte_queue(FILE *fp) {
    size_t len;
    size_t space;
    unsigned int i;
    unsigned int j;
    size_t total_len;
    unsigned int buff2_remaining;

    total_len = 0;
    /* Pull in a loop until there is less than BUFFSIZE space in thequeue */
    do {
        space = QUEUESIZE-queue_size; /* Dont pull more data than needed */
        if (space > BUFFSIZE) space = BUFFSIZE;
        
        /* ("  queue: space=%d\n",space); */
        len = fread(buffer, (size_t)1,(size_t)space, fp);
        if (len==0) {
            /* printf("  queue: len = %d\n",len); */
            return total_len;
        }

        /* Convert hex buffer to binary if we are in hex mode */
        if (hexmode == 1) len = hex2bin(buffer,len); 
        
        /* Fold upper case to lower */
        if (fold==1) {
            for (i=0;i<len;i++) {
                buffer[i]=tolower(buffer[i]);
            }
        }

        /* If we aren't doing word reverse, just move the buffer to the queue */ 
        if (word_reverse == 0) {
            /*  Transfer buffer to queue */
            for (i=0;i<len;i++) {
                if (byte_reverse == 1) {
                    queue[(queue_end+i) % QUEUESIZE] = byte_reverse_table[buffer[i]];
                } else {
                    queue[(queue_end+i) % QUEUESIZE] = buffer[i];
                }
                
                /* Call the monte carlo update that operated over bytes, not symbols */
                update_monte_carlo(buffer[i]);
            }
            
            filebytes += len;
            
            queue_size += len;
            queue_end = ((queue_end + len) % QUEUESIZE);

            total_len += len;
        } else {  /* word_reverse == 1  so use the word reverse buffer */
            
            /* Transfer buffer to word reverse queue */
            /*printf(" Transfer buffer to word reverse queue \n");*/
            for (i=0;i<len;i++) {
                buffer2[buffer2_size++]=buffer[i];
            }

            /* Pad out to a 4 byte boundary if we are at the end */
            /* printf(" Pad out to a 4 byte boundary if we are at the end \n");*/
            /*if ((buffer2_size % 4 != 0)) {
                printf(" Padding! \n");
                for (i=0;i<(4-buffer2_size);i++) {
                    buffer2[buffer2_size++] = 0x00;
                }
                fprintf(stderr,"Warning: Padded %d extra zeroes to make 4 byte boundary for word reverse\n",(4-buffer2_size));
            }*/

            /* Transfer it back out to the queue in blocks of 4 bytes in reverse */
            /* printf(" Transfer it back out to the queue in blocks of 4 bytes in reverse\n");*/
            /*if ((buffer2_size % 4) != 0) {
                printf("ERROR: buffer2_size = %d\n",buffer2_size);
                printf("ERROR: buffer size not on 4 byte boundary"); 
                exit(1);
            }*/
            buff2_remaining = buffer2_size;
            i = 0;
            do {
                /* printf( "    i = %d,  buffer2_size = %d \n",i,buffer2_size);*/
                for (j=0;j<4;j++) {
                    /* printf( "    j = %d\n",j); */
                    if (byte_reverse == 1) {
                        queue[(queue_end+j) % QUEUESIZE] = byte_reverse_table[buffer2[(i*4)+(3-j)]];
                        update_monte_carlo(byte_reverse_table[buffer2[(i*4)+(3-j)]]);
                    } else {
                        queue[(queue_end+j) % QUEUESIZE] = buffer2[(i*4)+(3-j)];
                        update_monte_carlo(buffer2[(i*4)+j]);
                    }
                }
                buff2_remaining -= 4;
                buffer2_size -= 4;
                queue_size += 4;
                queue_end = (queue_end+4) % QUEUESIZE;
                total_len += 4;
                i++;
            } while (buffer2_size>3);

            /* Pad any leftover */
            if (buffer2_size != 0) {
                for (j=0;j<4;j++) {
                    if (j < buffer2_size) {
                        queue[(queue_end+j) % QUEUESIZE] = 0x00;
                        update_monte_carlo(0x00);
                    } else {
                        if (byte_reverse == 1) {
                            queue[(queue_end+j) % QUEUESIZE] = byte_reverse_table[buffer2[(i*4)+(3-j)]];
                            update_monte_carlo(byte_reverse_table[buffer2[(i*4)+(3-j)]]);
                        } else {
                            queue[(queue_end+j) % QUEUESIZE] = buffer2[(i*4)+(3-j)];
                            update_monte_carlo(buffer2[(i*4)+(3-j)]);
                        }
                    }
                }
                fprintf(stderr,"Warning: Padded %d extra zeroes to make 4 byte boundary for word reverse\n",buffer2_size);
            }
            buffer2_size = 0;
        }
    } while ((QUEUESIZE-queue_size) > BUFFSIZE);
    return total_len;
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

void init_occurrences() {
    uint64_t i;
    
    no_occurrence_space = 0;
    
    occurrence_total = 0;
    if (symbol_length > 32) {
        fprintf(stderr,"Error, symbol length cannot be longer than 32 bits for occurrence count table\n");
        exit(1);
    }
    occurrence_size = ipow(2,symbol_length);
    fflush(stdout);
    occurrence_count = (uint64_t *) malloc (sizeof(uint64_t)*occurrence_size);
    /* printf("mallocating %lld bytes\n", (sizeof(uint64_t)*occurrence_size));
     */
    if (occurrence_count == NULL) {
        #ifdef _WIN32
        fprintf(stderr,"Warning, unable to allocate %lld bytes of memory for the occurrence count\n",(sizeof(uint64_t)*occurrence_size));
        #elif __llvm__
        fprintf(stderr,"Warning, unable to allocate %lld bytes of memory for the occurrence count\n",(sizeof(uint64_t)*occurrence_size));
        #elif __linux__
        fprintf(stderr,"Warning, unable to allocate %ld bytes of memory for the occurrence count\n",(sizeof(uint64_t)*occurrence_size));
        #endif
        no_occurrence_space = 1;
    }

    for (i=0;i<occurrence_size;i++) occurrence_count[i] = 0;
};

void init_longest() {
    uint64_t i;
    
    no_longest_space = 0;
    
    longest_total = 0;
    if (symbol_length > 32) {
        fprintf(stderr,"Error, symbol length cannot be longer than 32 bits for longest count table\n");
        exit(1);
    }
    longest_size = ipow(2,symbol_length);
    fflush(stdout);
    longest_count = (uint64_t *) malloc (sizeof(uint64_t)*occurrence_size);
    /* printf("mallocating %lld bytes\n", (sizeof(uint64_t)*occurrence_size));
     */
    if (longest_count == NULL) {
        #ifdef _WIN32
        fprintf(stderr,"Warning, unable to allocate %lld bytes of memory for the longest run table\n",(sizeof(uint64_t)*longest_size));
        #elif __llvm__
        fprintf(stderr,"Warning, unable to allocate %lld bytes of memory for the longest run table\n",(sizeof(uint64_t)*longest_size));
        #elif __linux__
        fprintf(stderr,"Warning, unable to allocate %ld bytes of memory for the longest run table\n",(sizeof(uint64_t)*longest_size));
        #endif
        no_longest_space = 1;
    }

    for (i=0;i<longest_size;i++) longest_count[i] = 0;

    longest_last_symbol=0;
    longest_run=0;
    longest_longest=0;
    longest_longest_symbol=0;
};

void init_chisq() {
    int i;
    chisq = 0.0;
    chisq_prob = (double *) malloc (sizeof(double)*occurrence_size);
    /* printf("mallocating %lld bytes for chisq probability table\n", (sizeof(double)*occurrence_size));
    */
    if (chisq_prob == NULL) {
        exit(1);
    }
    for (i=0;i<occurrence_size;i++) chisq_prob[i] = 0.0;
};

void init_filesize() {
    /* Nothing to do here */
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

void init_otherscc() {
    aeqb_count = 0;
    mean_count = 0;
}

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

void update_occurrences(uint64_t symbol) {
    occurrence_count[symbol]++;
    occurrence_total++;
};

void update_longest(uint64_t symbol) {
    if (symbol == longest_last_symbol) {
        longest_run++;
        if (longest_run > longest_count[symbol]) {
            longest_count[symbol] = longest_run;
        }
        if (longest_run > longest_longest) {
            longest_longest = longest_run;
            longest_longest_symbol = symbol;
        }
    } else {
        longest_run=1;
        longest_last_symbol=symbol;
    }
}

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
    int i;
    if (lagn==1) {
        /* We need lagn+1 symbols to start, so skip the first symbol(s) */
        scc_count++;
    
        if (scc_first==1) {
            scc_first = 0;
            first_symbol = symbol;
        } else {
            t1 += (scc_previous * symbol);
            if (scc_previous == symbol) aeqb_count += 1; // Other SCC
        }
        mean_count += symbol; //Other SCC
        t2 += symbol*symbol;
        t3 += symbol;
        
    
        /* printf("symbol %02X, count=%llu t1= %llX,  t2= %llx, t3= %llx\n",symbol,scc_count,t1,t2,t3); */
        scc_previous = symbol;
    } else { /* lagn > 1 */
        scc_count++;
    
        if (scc_count <= lagn) {
            scc_fifo[scc_count-1]=symbol;
            scc_first_lagn[scc_count-1]=symbol;
        } else {
            t1 += (scc_fifo[0] * symbol);
            if (scc_fifo[0] == symbol) aeqb_count += 1;
            for(i=0;i<lagn;i++) {
                scc_fifo[i]=scc_fifo[i+1];
            }
            mean_count += symbol;
            scc_fifo[lagn]=symbol;
            t2 += symbol*symbol;
            t3 += symbol;           
        }
        /* printf("symbol %02X, count=%llu t1= %llX,  t2= %llx, t3= %llx\n",symbol,scc_count,t1,t2,t3); */
         
    }
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
	for (eloop = 0; eloop < occurrence_size; eloop++) {
		if (chisq_prob[eloop] > 0.0) {
			ent += (chisq_prob[eloop] * log10(1.0 / chisq_prob[eloop]) *  3.32192809488736234787);
		}
	}

	result_entropy = ent;
	return;

	if (terse == 1) printf("%f,", ent);
	else printf("   Shannon Entropy = %f\n", ent);
};

void finalize_occurrences() {
    unsigned int i;
    unsigned int maxc;
    unsigned int maxsymbol;
    double maxp;
    double maxp_ent;

    /* Find the most frequent symbol */
    maxc=0;
    maxsymbol=0;
    for (i=0;i<occurrence_size;i++) {
        if (occurrence_count[i] > maxc) {
            maxc = occurrence_count[i];
            maxsymbol = i;
        }
    }

    //printf("maxc: %f\n",(double)maxc);
    //printf("occurance_size: %f\n",(double)occurrence_size);
    //printf("occurance_total: %f\n",(double)occurrence_total);
    maxp = ((double)maxc)/((double)occurrence_total);
    //printf("maxp: %f\n",maxp);
    maxp_ent = (-log10(maxp)/log10(2))/symbol_length;
    //printf("maxp_ent: %f\n",maxp_ent);
    result_min_entropy = maxp_ent;
    result_min_entropy_symbol = maxsymbol;

	if (terse != 1) {
        printf("   Min Entropy (by max occurrence of symbol %x) = %f\n", maxsymbol, maxp_ent);
    }
};

void finalize_longest() {
    result_longest_pvalue = longest_run_cdf((unsigned int)longest_longest, (unsigned int)symbol_count); 
}

void finalize_chisq() {
    uint64_t i;
    double diff;
    double chisq_final_prob;
    
    double expected;
    expected = (double)occurrence_total / (double)occurrence_size;
    for (i=0; i < occurrence_size; i++) {
        diff = (double)(occurrence_count[i]) - expected;
        chisq_prob[i] = ((double)occurrence_count[i])/occurrence_total;
        chisq      += (diff*diff)/expected;
        chisq_sum  += (double)(i * occurrence_count[i]);
    }
   
    chisq_final_prob = chisqp(chisq, (occurrence_size-1)); 
	result_chisq_count = occurrence_total;
	result_chisq_distribution = chisq;
	result_chisq_percent = chisq_final_prob * 100;

	return;
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

};

void finalize_compression() {
	double compression;

	compression = (100.0 * (symbol_length - ent)) / symbol_length;

	result_compression = compression;

	return;

};

void finalize_scc() {
    double scc;
    int64_t top;
    int64_t bottom;
    int i;

    double paeqb;
    double bias;

    if (scc_wrap==1) {
        if (lagn==1) {
            t1 += (scc_previous * first_symbol);
            t2 += first_symbol*first_symbol;
            t3 += first_symbol;
        } else {
            for (i=0;i<lagn;i++) {
                t1 += (scc_fifo[i] * scc_first_lagn[i]);
                t2 += scc_first_lagn[i]*scc_first_lagn[i];
                t3 += scc_first_lagn[i];
            }
        }
    } else {
        scc_count -= lagn;
    }

    /* need signed arithmetic because we are subtracting */
    top = (int64_t)(scc_count * t1) - (int64_t)(t3*t3);
    bottom = (int64_t)(scc_count * t2) - (int64_t)(t3*t3);
    scc = (double)top/(double)bottom;

	result_scc = scc;

    // This computation is to try to use the A=B count
    // The bias masks the serial correlation
    // Hence the * 1 - (2*abs(bias-0.5)) part.

    bias = (double)t3/(double)scc_count;
    paeqb = (double)aeqb_count/(double)scc_count;

    // Conversion assuming no bias.
    other_scc = ((2*paeqb)-1); 
    // Adjust using the bias - more bias pulls SCC towards 0.
    other_scc = other_scc * pow((1.0-(2.0*fabs(bias-0.5))),2);

    // Debugging output for A==B computation.
    //printf("OtherSCC\n     BIAS = %f\n     paeqb = %f\n     other_scc = %f\n",bias,paeqb,other_scc); 
    //printf("     (1.0-(2.0*fabs(bias-0.5))) = %f\n",(1.0-(2.0*fabs(bias-0.5))));
    //printf("     full other scc = %f\n",(((2*paeqb)-1)*  (1.0-(2.0*fabs(bias-0.5)))));
	return;
};

/* Visual Studio C doesnt have a regex library. So this does the
 * pattern search instead so I can compile on windows, linux and macos.
 */


/* look for vpattern in str. Return the match to found. Return True if found */
int find_vpattern(char *str,char *found) {
    size_t len;
    int i;
    int start;
    int end;
    int pos;
    int state;
    int done;
    char c;
    len = strlen(str);
    start = 0;
    end = 0;
    done = 0;

    /* A little state machine to match the _<int>p<int>V_ pattern */
    pos = 0;
    state = 1;
    done = 0;
    do {
        c = str[pos];
        if (state == 1) { /* _ */
            if (c=='_') {
                state++;
                start=pos;
            }
            pos++;
        } else if (state == 2) { /* first int */
            if (isdigit((char)c)) {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 3) { /* rest of int */
            if (isdigit((char)c)) {
                ; /* stay here */
            } else if ((c=='p') || (c=='.')) { /* decimal point */
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 4) { /* first int */
            if (isdigit((char)c)) {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 5) { /* rest of int */
            if (isdigit((char)c)) {
                ; /* stay here */
            } else if (c=='V') { /* V */
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 6) { /* _ */
            if ((c=='_') || (c=='.')) { // Allow 1p0V.bin  instead of 1p0V_.bin.
                done = 1;
                end = pos;
            } else {
                state = 1;
            }
            pos++;
        }

    } while ((pos < len) && (done == 0));
    
    if (done == 0) return 0;

    for(i=start;i<=end;i++) {
        found[i-start]=str[i];
    }
    found[i-start] = 0x00;
    return 1;
   
}

/* look for tpattern in str. Return the match to found. Return True if found */
int find_tpattern(char *str,char *found) {
    size_t len;
    int i;
    int start;
    int end;
    int pos;
    int state;
    int done;
    char c;
    len = strlen(str);
    start = 0;
    end = 0;
    done = 0;

    /* A little state machine to match the _<int>p<int>C_ pattern */
    pos = 0;
    state = 1;
    done = 0;
    do {
        c = str[pos];
        if (state == 1) { /* _ */
            if (c=='_') {
                state++;
                start=pos;
            }
            pos++;
        } else if (state == 2) { /* first int */
            if (isdigit((char)c) || ((char)c == '-')) {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 3) { /* rest of int */
            if (isdigit((char)c)) {
                ; /* Stay here */
            } else if ((c=='p') || (c=='.')) { /* decimal point */
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 4) { /* first int */
            if (isdigit((char)c)) {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 5) { /* rest of int */
            if (isdigit((char)c)) {
                ; /* Stay here */
            } else if (c=='C') { /* C */
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 6) { /* _ */
            if ((c=='_') || (c=='.')) { // Allow 10p0C.bin  instead of 10p0C_.bin.
                done = 1;
                end = pos;
            } else {
                state = 1;
            }
            pos++;
        }

    } while ((pos < len) && (done == 0));
    
    if (done == 0) return 0;

    for(i=start;i<=end;i++) {
        found[i-start]=str[i];
    }
    found[i-start] = 0x00;
    return 1;
   
}


/* look for cidpattern in str. Return the match to found. Return True if found */
int find_cidpattern(char *str,char *found) {
    size_t len;
    int i;
    int start;
    int end;
    int pos;
    int state;
    int done;
    char c;
    len = strlen(str);
    start = 0;
    end = 0;
    done = 0;

    /* A little state machine to match the _<int>p<int>C_ pattern */
    pos = 0;
    state = 1;
    done = 0;
    do {
        c = str[pos];
        if (state == 1) { /* _ */
            if (c=='_') {
                state++;
                start=pos;
            }
            pos++;
        } else if (state == 2) {
            if (c=='C') state++;
            else state = 1;
            pos++;         
        } else if (state == 3) {
            if (c=='I') state++;
            else state = 1;
            pos++;         
        } else if (state == 4) {
            if (c=='D') state++;
            else state = 1;
            pos++;         
        } else if (state == 5) {
            if (c=='-') state++;
            else state = 1;
            pos++;         
        } else if (state == 6) { /* first char of ID */
            if (c != '_') {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 7) { /* rest of ID */
            if (c != '_') {
                ; /* Stay here */
            } else { /* _ */
                done = 1;
                end = pos;
            } 
            pos++;
        }

    } while ((pos < len) && (done == 0));
    
    if (done == 0) return 0;

    for(i=start;i<=end;i++) {
        found[i-start]=str[i];
    }
    found[i-start] = 0x00;
    return 1;
   
}

/* look for procpattern in str. Return the match to found. Return True if found */
int find_procpattern(char *str,char *found) {
    size_t len;
    int i;
    int start;
    int end;
    int pos;
    int state;
    int done;
    char c;
    len = strlen(str);
    start = 0;
    end = 0;
    done = 0;

    /* A little state machine to match the _PROC-<name>_ pattern */
    pos = 0;
    state = 1;
    done = 0;
    do {
        c = str[pos];
        if (state == 1) { /* _ */
            if (c=='_') {
                state++;
                start=pos;
            }
            pos++;
        } else if (state == 2) {
            if (c=='P') state++;
            else state = 1;
            pos++;         
        } else if (state == 3) {
            if (c=='R') state++;
            else state = 1;
            pos++;         
        } else if (state == 4) {
            if (c=='O') state++;
            else state = 1;
            pos++;         
        } else if (state == 5) {
            if (c=='C') state++;
            else state = 1;
            pos++;         
        } else if (state == 6) {
            if (c=='-') state++;
            else state = 1;
            pos++;         
        } else if (state == 7) { /* first char of ID */
            if (c != '_') {
                state++;
            } else {
                state = 1;
            }
            pos++;
        } else if (state == 8) { /* rest of ID */
            if (c != '_') {
                ; /* Stay here */
            } else { /* _ */
                done = 1;
                end = pos;
            } 
            pos++;
        }

    } while ((pos < len) && (done == 0));
    
    if (done == 0) return 0;

    for(i=start;i<=end;i++) {
        found[i-start]=str[i];
    }
    found[i-start] = 0x00;
    return 1;
   
}

void parse_the_filename(char *filename) {

    char match[256];
    int i;

    if (find_vpattern(filename,match)) {
        for (i=0;i<strlen(match);i++) {
            if (match[i]=='p') match[i] = '.';
        }
        sscanf(match,"_%lfV_",&voltage);
    } else {
        fprintf(stderr,"Regex error scanning for _<num>p<num>V_:\n");
        voltage = 0.0;
    }


    if (find_tpattern(filename,match)) {
        for (i=0;i<strlen(match);i++) {
            if (match[i]=='p') match[i] = '.';
        }
        sscanf(match,"_%lfC_",&temperature);
    } else {
        fprintf(stderr,"Regex error scanning for _<num>p<num>C_:\n");
        temperature = 0.0;
    }

    if (find_cidpattern(filename,match)) {
        match[strlen(match)-1]=0x00;
        sscanf(match,"_CID-%s",(char *)&deviceid);
    } else {
        fprintf(stderr,"Regex error scanning for _CID-<ID>__:\n");
    }

    if (find_procpattern(filename,match)) {
        match[strlen(match)-1]=0x00;
        sscanf(match,"_PROC-%s",(char *)&process);
    } else {
        fprintf(stderr,"Regex error scanning for _PROC-<ID>__:\n");
    }
}


/********
* main() is mostly about parsing and qualifying the command line options.
*/

int main(int argc, char** argv)
{
    int i;

    /* Defaults */
    symbol_length = 8;
    hexmode = 1;
    print_occurrence = 0;
    print_longest = 0;
    fold = 0;
    terse = 0;
    use_stdin = 1;
    fp = NULL;
    terse_index = 0;
    scc_wrap = 0;
    lagn = 1;
    using_inputlistfile = 0;
    suppress_header = 0;
    byte_reverse = 0;
    parse_filename = 0;
    word_reverse = 0; 
    buffer2_size = 0;
    ent_exact = 0;
	#define ERRSTRINGSIZE 256
    #ifdef _WIN32
	errno_t err;
	char errstring[ERRSTRINGSIZE];
    #endif
    int filenumber = 0;
    
    int got_symbol_length=0;
    
    char optString[] = "bprRcCwftehusi:n:l:";
    int longIndex;
    static const struct option longOpts[] = {
    { "symbol_length", required_argument, NULL, 'l' },
    { "binary", no_argument, NULL, 'b' },
    { "byte_reverse", no_argument, NULL, 'r' },
    { "word_reverse", no_argument, NULL, 'R' },
    { "occurrence", no_argument, NULL, 'c' },
    { "fold", no_argument, NULL, 'f' },
    { "parse_filename", no_argument, NULL, 'p' },
    { "inputlistfile", required_argument, NULL, 'i' },
    { "scc_wrap", no_argument, NULL, 'w' },
    { "lagn", required_argument, NULL, 'n' },
    { "terse", no_argument, NULL, 't' },
    { "ent_exact", no_argument, NULL, 'e' },
    { "suppress_header", no_argument, NULL, 's' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
    };

    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    while( opt != -1 ) {
        switch( opt ) {
            case 'b':
                if (got_symbol_length == 0) {
                    symbol_length = 1; /* binary mode treats newlines as data */
                }
                hexmode = 0;
                break;
                
            case 'l':
                symbol_length = atoi(optarg); /* -b -l <n> is valid. symbol length
                                               * will be <n> and newlines will be treated
                                               * as data.
                                               */
                got_symbol_length = 1;
                
                break;
            case 'i':
                strncpy(inputlistfilename,optarg,255);
                using_inputlistfile = 1;
                break;
 
            case 'c':
                print_occurrence = 1;
                break;
            
            case 'C':
                print_longest = 1;
                break;
 
            case 'p':
                parse_filename = 1;
                break;

            case 'r':
                byte_reverse = 1;
                break;

            case 'R':
                word_reverse = 1;
                break;
                            
            case 'w':
                scc_wrap = 1;
                break;
                
            case 'n':
                lagn = atoi(optarg);
                break;
                                
            case 'f':
                fold = 1;
                break;
            
            case 't':
                terse = 1;
                break;
           
            case 'e':
                ent_exact = 1;
                break;
 
            case 's':
                suppress_header = 1;
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

    if ((fold==1) && (symbol_length != 8)) {
            fprintf(stderr,"Error: Fold must be used with 8 bit word size\n");
            exit(1);
    }
        
	if (symbol_length < 1) {
		fprintf(stderr,"Error: Symbol length %d must not be 0 or negative. \n",symbol_length);
		exit(1);
	}

    init_byte_queue();
    
    /* Loop through the filenames */
	if ((optind==argc) && (using_inputlistfile==0)) {
		use_stdin = 1;
	}
	else {
		use_stdin = 0;
	}

    if ((parse_filename==1) && (use_stdin==1)) {
        fprintf(stderr,"Error: Can't parse filename when using stdin for input\n");
        exit(1);
    }

	terse_index = 0;
	/* if (use_stdin == 0) { */
	filenumber = optind;

    char *filelist;
    int lines;
    int lineno;
    FILE *ifp;
    int filenamecount = 0;
    char * res;
    char line[256];

    filelist = (char*)0;
    /* build the list of filenames from the input list file */
    if (using_inputlistfile==1) {
        lines = count_lines_in_file(inputlistfilename);
        ifp = fopen(inputlistfilename,"r");
        if (ifp==NULL) {
            fprintf(stderr,"Error: Cannot open %s for reading\n",inputlistfilename);
            exit(1);
        }
        
        filelist = (char *)malloc(sizeof(char *)*256*lines);
        
        if (filelist==NULL) {
            fprintf(stderr,"Error: Cannot allocate memory for filename list from input file list file %s\n",inputlistfilename);
            exit(1);
        }
         
        for (lineno=0;lineno<lines;lineno++) {
            res = fgets(line, 256, ifp);
            if (res != NULL) {
                /* mute the newlines from the file*/
                for (i = 0;i<256;i++) {
                    if (line[i]=='\n') line[i]=0;
                }

                /* Grab the file names, skipping the ones beginning with # */
                strncpy(&(filelist[256*filenamecount]),line,256);
                /*printf(" Scanning file list, got :%s:\n",&(filelist[256*filenamecount]));*/
                if (line[0]!='#') filenamecount++; /* ignore lines beginning with # */ 
            }
        }
        fclose(ifp);
        if (filenamecount==0) {
            fprintf(stderr,"Error: Did not file any filenames in input file list file %s\n",inputlistfilename);
            exit(1);
        }
        filenumber = 0; 
    }

	do {
		terse_index++;
		filebytes = 0;
        /* printf("OPTIND %d, filenumber %d, ARGC %d\n",optind,filenumber,argc); */
		if (use_stdin==1) {
			use_stdin = 1;
			if (hexmode != 1) freopen(NULL, "rb", stdin);
			fp = stdin;
		}
		else {
            if (using_inputlistfile==0) {
			    filename = argv[filenumber];
            } else {
                /* Get file from input file list file */ 
                filename = &(filelist[256*filenumber]);
                /*printf("FILENUMBER = %d , Filename = %s\n",filenumber,filename);*/
            }
            if (parse_filename==1) parse_the_filename(filename);

			if (hexmode == 1) {
				if ((terse == 0) && (ent_exact == 0))printf(" opening %s as hex text\n", filename);
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
				if ((terse == 0) && (ent_exact==0)) printf(" opening %s as binary\n", filename);
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

			if ((terse == 0) && (ent_exact == 0))printf(" Symbol Size(bits) = %d\n", symbol_length);

		}
/*
   0,  File-bytes,    Entropy, Min_entropy, MinEntropy-Symbol,     Chi-square,  Mean,        Monte-Carlo-Pi, Serial-Correlation, Filename
   1,      210109,    7.942742,    2.349834,55,    0.000000,  127.497723,    3.240505,       0.002354,           pt1a.bin
*/
		/* Print terse header if necessary */
		if ((terse == 1) && (terse_index == 1) && (suppress_header==0)) {
            if (ent_exact == 1) {
                if (symbol_length==1) {
                    printf("0,File-bits,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n");
                } else {
                    printf("0,File-bytes,Entropy,Chi-square,Mean,Monte-Carlo-Pi,Serial-Correlation\n");
                }
            }
            else if (parse_filename==1) {
			    printf("   0,  File-bytes,     CID, Process, Voltage,    Temp,     Entropy,  MinEntropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation, Filename, Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue\n");
            } else {
			    printf("   0,  File-bytes,    Entropy,  Min_entropy, MinEntropy-Symbol,  Chi-square,        Mean, Monte-Carlo-Pi, Serial-Correlation, Filename, Longest-Run-Symbol, Longest-Run-Length, Longest-Run-PValue\n");
		    }
        }

		/* Initialize the metrics */
		symbol_count = 0;

        fflush(stdout);
		init_mean();
		init_entropy();
		init_occurrences();
		init_longest();
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
				not_eof = (int)fill_byte_queue(fp); /* get bytes from file into queue */
				if (not_eof == 0) break;
			}
			symbol = get_symbol(symbol_length);      /* Pull a symbol from the queue */
			symbol_count++;

			/* Finish up if no symbols left in queue */
			if (symbol == -1) break;

			/* Then update the algorithms using the symbol */
			update_mean(symbol);
			update_entropy(symbol);
			if (no_occurrence_space == 0) update_occurrences(symbol);
			if (no_longest_space == 0) update_longest(symbol);
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
		if (no_occurrence_space == 0) finalize_occurrences();
        if (no_longest_space == 0) finalize_longest();
		finalize_chisq();
		finalize_entropy();
		finalize_filesize();
		finalize_monte_carlo();
		finalize_compression();
		finalize_scc();

		if (terse == 1) {
            if (ent_exact==1) {
                if (symbol_length == 1) {
                printf("%d,%"PRIu64",%f,%f,%f,%f,%f\n",terse_index,filebytes*8,result_entropy,result_chisq_distribution,result_mean,result_pi,result_scc);
                } else {
                printf("%d,%"PRIu64",%f,%f,%f,%f,%f\n",terse_index,filebytes,result_entropy,result_chisq_distribution,result_mean,result_pi,result_scc);
                }
            }
            else if ((parse_filename==1) && (symbol_length==1)) {
                printf("%4d,%12"PRIu64",%8s,%8s,%8.2f,%8.2f,%12f,%12f,%"PRIu32",%12f,%12f,%15f,   %16f, %s, %"PRIx64", %"PRIu64", %f\n", terse_index, filebytes, deviceid,process,voltage,temperature,result_entropy, result_min_entropy,result_min_entropy_symbol, result_chisq_percent, result_mean, result_pi, result_scc, filename, longest_longest_symbol,longest_longest,result_longest_pvalue);

		    } else if ((parse_filename==0) && (symbol_length==1)) {
                printf("%4d,%12"PRIu64",%11f, %12f,%18"PRIx32",%12f,%12f,%15f,       %12f, %s, %"PRIx64", %"PRIu64", %f\n", terse_index, filebytes, result_entropy, result_min_entropy,result_min_entropy_symbol, result_chisq_percent, result_mean, result_pi, result_scc, filename, longest_longest_symbol,longest_longest,result_longest_pvalue);
            }

            else if ((parse_filename==1) && (symbol_length!=1)) {
                printf("%4d,%12"PRIu64",%8s,%8s,%8.2f,%8.2f,%12f,%12f,%18"PRIx32",%12f,%12f,%15f,   %16f, %s, %"PRIx64", %"PRIu64", (null)\n", terse_index, filebytes, deviceid,process,voltage,temperature,result_entropy, result_min_entropy,result_min_entropy_symbol, result_chisq_percent, result_mean, result_pi, result_scc, filename, longest_longest_symbol,longest_longest);
		    
            } else if ((parse_filename==0) && (symbol_length!=1)) {
                printf("%4d,%12"PRIu64",%11f, %12f,%18"PRIx32",%12f,%12f,%15f,       %12f, %s, %"PRIx64", %"PRIu64", (null)\n", terse_index, filebytes, result_entropy, result_min_entropy,result_min_entropy_symbol, result_chisq_percent, result_mean, result_pi, result_scc, filename, longest_longest_symbol,longest_longest);
            }
        }
		else {

            /* Output the occurrence count if requested */
            if ((print_occurrence==1) && (no_occurrence_space == 0) ) {
                double fraction;
                for (i=0; i<occurrence_size;i++) {
                    fraction = (double)occurrence_count[i]/(double)occurrence_total;
                    #ifdef _WIN32
                    printf("   Value %4d , frequency=%llu , fraction=%f\n", i, occurrence_count[i], fraction);
                    #elif __llvm__
                    printf("   Value %4d , frequency=%llu , fraction=%f\n", i, occurrence_count[i], fraction);
                    #elif __linux__
                    printf("   Value %4d , frequency=%lu , fraction=%f\n", i, occurrence_count[i], fraction);
                    #endif
                }
            }
            
            /* Output the occurrence count if requested */
            if ((print_longest==1) && (no_longest_space == 0) ) {
                for (i=0; i<occurrence_size;i++) {
                    #ifdef _WIN32
                    printf("   Symbol %x , Longest Run=%"PRIu64"\n", i, longest_count[i]);
                    #elif __llvm__
                    printf("   Symbol %x , Longest Run=%"PRIu64"\n", i, longest_count[i]);
                    #elif __linux__
                    printf("   Symbol %x , Longest Run=%"PRIu64"\n", i, longest_count[i]);
                    #endif
                }
            }

            /* Output the non terse results */
            if (ent_exact == 1) {
                /* Make it look like this:
                 *
                 * Entropy = 4.676598 bits per byte.
                 *
                 * Optimum compression would reduce the size
                 * of this 57737 byte file by 41 percent.
                 * 
                 * Chi square distribution for 57737 samples is 1499055.27, and randomly
                 * would exceed this value less than 0.01 percent of the times.
                 *
                 * Arithmetic mean value of data bytes is 71.6317 (127.5 = random).
                 * Monte Carlo value for Pi is 4.000000000 (error 27.32 percent).
                 * Serial correlation coefficient is 0.515629 (totally uncorrelated = 0.0).
                 */
                
                printf("Entropy = %f bits per byte.\n\n",(result_entropy *(8.0/symbol_length)));
                printf("Optimum compression would reduce the size\n");
                printf("of this %d byte file by %d percent\n\n",(int)filebytes,(int)result_compression);
                printf("Chi square distribution for %d samples is %f, and randomly\n",(int)result_chisq_count,result_chisq_distribution);
                printf("would exceed this value less than %f percent of the times.\n\n",result_chisq_percent);
                printf("Arithmetic mean value of data bytes is %f (127.5 = random).\n",result_mean);
                printf("Monte Carlo value for Pi is %f (error %f percent).\n",result_pi,result_pierr);
                printf("Serial correlation coefficient is %f (totally uncorrelated = 0.0).\n",result_scc);
            }
            else {
                if (parse_filename==1) {
                printf("   Device ID   : %s\n",deviceid);
                printf("   Process     : %s\n",process);
                printf("   Voltage     : %0.2lfV\n",voltage);
                printf("   Temperature : %0.2lfC\n",temperature);
                }
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
                printf("   Longest Run Symbol = %"PRIx64". Run Length = %"PRIu64"\n",longest_longest_symbol,longest_longest);
                if (symbol_length == 1) printf("   Probabilty of longest run being <= %"PRIu64" = %f\n",longest_longest,result_longest_pvalue);
                //printf("SCC by A=B Count is %f (totally uncorrelated = 0.0).\n",other_scc);
            }
		}

        /* Free the per-loop mallocs */		
        free(occurrence_count);
        free(chisq_prob);
    
        if (fp != stdin) fclose(fp); 
		filenumber++;
	} while (
                    ((filenumber < argc) && (use_stdin != 1) && (using_inputlistfile!=1)) /* still going through argv */
                ||
                    ((filenumber < filenamecount) && (using_inputlistfile==1)) /* still going through file list file */
            );
    /* Free the per-run malloc */
    free(filelist);

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



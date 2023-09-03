

/* Visual Studio C doesnt have a regex library. So this does the
 * pattern search instead so I can compile on windows, linux and macos.
 */


#include <inttypes.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#include "filename_parse.h"

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



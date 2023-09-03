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

#include <inttypes.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>


#ifndef NO_GMP
#include <mpfr.h>
#endif

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

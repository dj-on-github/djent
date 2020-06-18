#include <math.h>
#include <stdio.h>
#include <stdarg.h>
#include <mpfr.h>

// Return the probability of the longest run of heads being less than or equal to n
// in a sequence of r uniform coin tosses.
double P(mpfr_t n,mpfr_t r) { // n=longest run. R = length of data sequence
    double answer;
    mpfr_t topa;
    mpfr_t bottoma;
    mpfr_t first;
    mpfr_t topb;
    mpfr_t nplusone;
    mpfr_t bottomb;
    mpfr_t nplus2over2;
    mpfr_t second;
    mpfr_t mpfrans;
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
    return answer;
}

int main() {
    unsigned int i;
    mpfr_t n;
    mpfr_t r;
    double answer;

    mpfr_set_default_prec(1024);
    mpfr_init(n);
    mpfr_init(r);

    for(i=2;i<32;i++) {
        mpfr_set_ui(n,i,MPFR_RNDN);
        mpfr_set_ui(r,1024*1024*8,MPFR_RNDN);
        
        answer=P(n,r);
        printf("n=%d, r=%d, P=%f\n",i,1024*1024*8,answer);
    }

    mpfr_clear(n);
    mpfr_clear(r);
    return 1;
}
     

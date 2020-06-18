#!/usr/bin/env python

from __future__ import print_function
from __future__ import division

import math
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision=1024

# Return the probability of the longest run of heads being n
# in a sequence of r coin tosses.
def P(n,r): # n=longest run. R = length of data sequence
    topa = -mpfr(r+1)
    bottoma = (mpfr(2)**mpfr(n+1))-n-2
    first = gmpy2.exp(topa/bottoma)

    topb = (mpfr(2)**mpfr(n+1))-1
    bottomb = (mpfr(2)**mpfr(n+1))-(mpfr(n+2)/mpfr(2))
    second = topb/bottomb

    answer = first*second
    return answer

for n in range(2,30):
    mn = mpfr(n)
    r = mpfr(256)
    answer = P(mn,r)

    print("n={n}, r={r}, P={answer}".format(**locals()))

     

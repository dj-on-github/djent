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

#ifndef NO_GMP
#include <mpfr.h>
#endif

// Return the probability of the longest run of heads being less than or equal to n
// in a sequence of r uniform coin tosses. Use MPFR to avoid overflows.
double longest_run_cdf(unsigned int ui_n,unsigned int ui_r);

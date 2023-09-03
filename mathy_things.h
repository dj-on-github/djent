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

#define LOG_SQRT_PI 0.5723649429247000870717135 /* log (sqrt (pi)) */
#define I_SQRT_PI   0.5641895835477562869480795 /* 1 / sqrt (pi) */
#define BIGX        20.0         /* max value to represent exp (x) */
#define ex(x)       (((x) < -BIGX) ? 0.0 : exp(x))

uint64_t ipow(uint64_t base, uint64_t exp);
double zcdf(double z);
double chisqp(double ax, size_t df);


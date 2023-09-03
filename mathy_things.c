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


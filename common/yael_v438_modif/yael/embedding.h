/*
Copyright © INRIA 2009-2014.
Authors: Matthijs Douze & Herve Jegou
Contact: matthijs.douze@inria.fr  herve.jegou@inria.fr

This file is part of Yael.

    Yael is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Yael is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Yael.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __embedding_h
#define __embedding_h


/* Dimensionality after polynomial embedding of a d-dimensional vector */
#define EMB_POLY2_D(d)  ((d)*((d)+1)/2)
#define EMB_POLY3_D(d)  (((d)*(d)*(d)+3*(d)*(d)+2*(d))/6)


/* polynomial embedding of degree 2: cos(a,b) -> cos(a,b)^2 
   d is dimensionality of x. 
   y should be pre-allocated with dimensionality d*(d+1)/2     */
void emb_poly2 (const float * x, float * y, int d);

/* Same as emb_poly2a, but add to existing vector y */
void emb_poly2a (const float * x, float * y, int d);

/* polynomial embedding of degree 3: cos(a,b) -> cos(a,b)^3 */
void emb_poly3 (const float * x, float * y, int d);

/* Same as emb_poly3a,, but add instead of mapping */
void emb_poly3a (const float * x, float * y, int d);

/* Fourier modulation with F frequencies */
void ang_modulate (const float * an, const float * x, const float theta, float * y, int d, int F);


#endif
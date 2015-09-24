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

#ifndef SPECTRAL_CLUSTERING_H_INCLUDED
#define SPECTRAL_CLUSTERING_H_INCLUDED

/*! @addtogroup clustering
 *  @{  */


/* perform a spectral clustering of the dataset v, 
   as proposed in [Ng Jordan Weiss 01]               */
double spectral_clustering (int d, int n, int k, double sigma, int niter,
			    const float * v, int nt, int seed, int nredo,
			    int * assign, int * nassign);


/*! @} */
#endif

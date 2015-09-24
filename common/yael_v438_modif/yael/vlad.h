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

#ifndef VLAD_H_INCLUDED
#define VLAD_H_INCLUDED

#include <stdlib.h>



/*! VLAD descriptor implementation
 *
 * @param k               nb of centroids
 * @param d               local descriptor dimension
 * @param centroids(d,k)  local descriptor centroids 
 * @param n               nb of local image descriptors
 * @param v(d,n)          local image descriptors
 * @param desc(d,k)       global VLAD descriptor for the image (output) 
 */
void vlad_compute (int k, int d, const float *centroids, 
                   int n, const float *v,
                   float *desc);

/*! weighted VLAD descriptor 
 * @param weights(n)      weight applied to each local descriptor before accumulation
 */
void vlad_compute_weighted (int k, int d, const float *centroids, 
                            int n, const float *v, const float *weights, 
                            float *desc);


/*! like vlad_compute, but compute vlads on subsets of the n descriptors. 
 *
 * To compute colorvlads, just set d=3, n=nb of pixels in image, v=all
 * pixels in image, and have the subsets correspond to the pixel
 * indices of local descriptor ellipses.
 *
 * @param n_subset                                nb of subsets
 * @param subset_indexes(subset_ends[n_subset-1]) subset i is defined as subset_indexes[subset_ends[i-1]...subset_ends[i]-1]
 * @param subset_ends                             ends (and beginnings) of the subset lists
 * @param desc(d,k,n_subset)                      VLAD descriptor for all subsets
 */
void vlad_compute_subsets(int k, int d, const float *centroids, 
                          int n, const float *v,
                          int n_subset,
                          const int *subset_indexes, 
                          const int *subset_ends,
                          float *desc); 


/*! Compute bag-of-features (not normalized, just count)
 *
 * @param k               nb of centroids
 * @param d               local descriptor dimension
 * @param centroids(d,k)  local descriptor centroids 
 * @param n               nb of local image descriptors
 * @param v(d,n)          local image descriptors
 * @param desc(d,k)       global VLAD descriptor for the image (output) 
 */
void bof_compute (int k, int d, const float *centroids, 
		  int n, const float *v, int *desc);


/*! Same as bof_compute, but with multiple assignment, and might be multi-threaded */
void bof_compute_ma (int k, int d, const float *centroids, 
		     int n, const float *v, int *desc, 
		     int ma, float alpha, int nt);


/*! like vlad_compute_subsets, but compute BOFs on subsets of the n descriptors 
 * instead of VLADs (intended for bag-of-colors computation)
 *
 * @param desc(k,n_subset)                      BOF descriptor for all subsets
 */
void bof_compute_subsets(int k, int d, const float *centroids, 
                          int n, const float *v,
                          int n_subset,
                          const int *subset_indexes, 
                          const int *subset_ends,
                          float *desc); 



#endif

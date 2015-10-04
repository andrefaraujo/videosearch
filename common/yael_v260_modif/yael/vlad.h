/*
Copyright © INRIA 2010-2011. 
Authors: Matthijs Douze & Herve Jegou 
Contact: matthijs.douze@inria.fr  herve.jegou@inria.fr

This software is a computer program whose purpose is to provide 
efficient tools for basic yet computationally demanding tasks, 
such as find k-nearest neighbors using exhaustive search 
and kmeans clustering. 

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
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

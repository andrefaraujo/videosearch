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

#ifndef KMEANS_H_INCLUDED
#define KMEANS_H_INCLUDED


/*---------------------------------------------------------------------------*/
/*! @addtogroup clustering
 *  @{  */


/* flags for kmeans */
#define KMEANS_QUIET                    0x10000
#define KMEANS_INIT_BERKELEY            0x20000
#define KMEANS_NORMALIZE_CENTS          0x40000
#define KMEANS_INIT_RANDOM              0x80000
#define KMEANS_INIT_USER               0x100000

#define KMEANS_L1                       0x200000
#define KMEANS_CHI2                     0x400000

/*! Compute the k-means centroids. 
 *
 * @param v(d, n)           vectors to cluster
 * @param centroids(d, k)   output centroids (input centroids here if KMEANS_INIT_USER)
 * @param flags             a set of computation parameters: 
 *                     - flags & 0xffff : use this many threads to compute 
 *                     - flags & KMEANS_QUIET: suppress kmeans output
 *                     - flags & KMEANS_INIT_RANDOM: random initialization 
 *                     - flags & KMEANS_NORMALIZE_CENTS: normalize centroids to L2=1 after they are computed
 *                     - flags & KMEANS_INIT_USER: the user gives the initialization
 *                     - flags & KMEANS_L1: L1 distance kmeans
 *                     - flags & KMEANS_CHI2: chi-squared distance kmeans
 *                          -> provided by user with parameter centroids_out)
 * @param seed              random seed for intialization (used only if !=0)
 * @param redo              perform clustering this many times and keep clusters with smallest quantization error
 * @param dis(n)            squared distance to assigned centroid of each input vector (may be NULL)
 * @param assign(n)         index of assigned centroid in 0..k-1 (may be NULL) 
 * @param nassign(k)        nb of vectors assigned to each centroid (may be NULL)
 *
 * @return final quantization error 
 */
float kmeans (int d, int n, int k, int niter, 
	      const float * v, int flags, long seed, int redo, 
	      float * centroids, float * dis, 
	      int * assign, int * nassign);

/*--- Following functions are for forward compatibility (and may be removed in the future) ---*/

/*! simplified call */
float* clustering_kmeans (int n, int d,
                          const float *points, int k, int nb_iter_max, 
                          double normalize);

/*! Same as kmeans, but generate in addition the assignment
 *  performed on the input set
 */
float* clustering_kmeans_assign (int n, int d,
				 const float *points, int k, int nb_iter_max, 
				 double normalize, 
				 int ** clust_assign_out);



float* clustering_kmeans_assign_with_score (int n, int d,
                                            const float *points, int k, int nb_iter_max,
                                            double normalize, int n_thread, double *score_out,
                                            int ** clust_assign_out);


/*! @} */
#endif

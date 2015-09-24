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

#include <assert.h> 
#include <stdlib.h>
#include <math.h>

#include "nn.h"
#include "eigs.h"
#include "vector.h"
#include "matrix.h"
#include "kmeans.h"
#include "spectral_clustering.h"

/* compute the Gaussian kernel */
void gaussian_kernel (int d, int n1, int n2, 
		      const float * v1,  const float * v2, 
		      double sigma, float * ker)
{
  int i;
  /*  double nr = 1 / (sigma * sqrt (3.141592653589793)); not used */
  double nr2 = 1 / (2 * sigma * sigma);

  /* first compute the L2 distances */
  compute_cross_distances (d, n1, n2, v1, v2, ker);

  for (i = 0 ; i < n1 * n2 ; i++)
    ker[i] = exp (-ker[i] * nr2); 
}


/* this spectral clustering implementation would advantageously 
   use the sparsity of the kernel matrix and arpack. */
double spectral_clustering (int d, int n, int k, double sigma, int niter,
			    const float * v, int nt, int seed, int nredo,
			    int * assign, int * nassign)
{
  int i, j;

  /* first compute the kernel */
  float * ker = fvec_new (n * n);
  gaussian_kernel (d, n, n, v, v, sigma, ker);

  fmat_print (ker, n, n);


  /* set diagonal to zero and perform the normalization */
  float * deg = fvec_new_0 (n);
  for (i = 0 ; i < n ; i++) {
    ker[i * n + i] = 0;
    for (j = 0 ; j < n ; j++)
      deg[i] += ker[i * n + j];
  }

  for (i = 0 ; i < n ; i++)
    deg[i] = 1 / sqrt (deg[i]);

  for (i = 0 ; i < n ; i++)
    for (j = 0 ; j < n ; j++)
      ker[i * n + j] *= deg[i] * deg[j];

  /* eigenvalue decomposition */
  float * eigval = fvec_new (n);
  float * eigvec = fvec_new (n * n);
  int ret = eigs_sym (n, ker, eigval, eigvec); 
  free (deg);
  free (ker);
  assert (ret==0);
  

  /* keep only the first eigenvectors */
  eigs_reorder (n, eigval, eigvec, 1);

  /* create the "embedded" vectors to normalize and cluster */
  float * e = fvec_new (n * k);

  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < k ; j++)
      e[i * k + j] = eigvec[j * n + i];
    fvec_normalize (e, k, 2);
  }

  free (eigval);
  free (eigvec);

  double err = kmeans (k, n, k, niter, e, nt, seed, nredo, 
		       NULL, NULL, assign, nassign);

  free (e);
  return err;
}

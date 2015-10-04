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

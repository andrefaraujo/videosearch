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

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "vector.h"
#include "sorting.h"
#include "matrix.h"
#include "kmeans.h"
#include "nn.h"
#include "machinedeps.h"



static void random_init(long d, int n, int k, const float * v, int * sel, 
                        unsigned int seed) {
  int *perm=ivec_new_random_perm_r(n,seed);
  ivec_cpy(sel,perm,k);
  free(perm);
}

static double drand_r(unsigned int *seed) {
  return rand_r(seed)/((double)RAND_MAX + 1.0);
}

/* the kmeans++ initialization (see wikipedia) */
static void kmeanspp_init (long d, int n, int k, const float * v, 
			   int * sel, int verbose, unsigned int seed, int n_thread)
{
  /* select the first centroid and set the others unitialized*/

  long i, j;
  /*
  for(i=0;i<k;i++) sel[i]=i;
  return;
  */
  float * disbest = fvec_new_set (n, HUGE_VAL);
  float * distmp = fvec_new (n);

  sel[0] = rand_r(&seed) % k;

  for (i = 1 ; i < k ; i++) {
    int newsel = sel[i - 1];
    
    if (verbose && i % 10 == 0) {
      printf ("%d/%d\r", (int)i, k); 
      fflush (stdout);
    }

    if(0) { /* simple and slow */
      for (j = 0 ; j < n ; j++) {
        distmp[j] = fvec_distance_L2sqr (v + d * j, v + d * newsel, d);
        if (distmp[j] < disbest[j])
          disbest[j] = distmp[j];
      }
    } 
    else { /* complicated and fast */
      compute_distances_1_thread(d, n, v + d * newsel, v, distmp, n_thread); 

      for (j = 0 ; j < n ; j++) 
        if (distmp[j] < disbest[j]) 
	  disbest[j] = distmp[j];
    }
    
    /* convert the best distances to probabilities */
    memcpy (distmp, disbest, n * sizeof (*distmp));
    fvec_normalize (distmp, n, 1);
    double rd = drand_r(&seed);
    
    for (j = 0 ; j < n - 1 ; j++) {
      rd -= distmp[j];
	if (rd < 0)
	  break;
    }
    
    sel[i] = j;
  }
  if (verbose)
    printf("\n");
  free (disbest);
  free (distmp);
}

/* Minimize 
     sum( (x - a[i])^2 / (x + a[i]), i = 0..n-1) 

in x. Return minimal x st. x > 0. 

Assumes a[i] > 0 for all i.

Uses Newton's method on the derivative of the expression, since the
function is convex.
*/

float minimize_sum_chi2(const float *a, int n) {

  if(n == 0) return 0.0/0.0;
  if(n == 1) return a[0]; 

  float x, prev_x, mag; 
  int i; 
  
  /* typical magnitude  */
  mag = 0; 
  for(i = 0; i < n; i++) mag += a[i];
  mag /= n;

  /* start at 0 */
  x = 0; 

  /* start loops */

  int niter = 0; 

  do {

    float d, dd; /* derivative and second derivative */

    d = dd = 0;

    for(i = 0; i < n ; i++) {
      
      float sum2 = (x + a[i]) * (x + a[i]);
      float sum3 = (x + a[i]) * sum2;
      
      if(sum2 == 0) continue;

      d += (x - a[i]) * (x + 3 * a[i]) / sum2;
      
      dd += 8 * a[i] * a[i] / sum3;
    }
    
    if( d == 0 ) /* ok, found minimum */
      break;

    if( dd == 0 ) {
      fprintf(stderr, "warn: minimize_sum_chi2 with zero second derivative, a = [");       
      for(i = 0; i < n; i++) fprintf(stderr, "%g ", a[i]); 
      fprintf(stderr, "]\n");
    }
    
    prev_x = x; 
    
    x -= d / dd; 
    
    if(x < 0) x = prev_x / 2.0; 
    
    if(++niter==1000) {
      fprintf(stderr, "warn: minimize_sum_chi2 reached 1000 iterations, a=[" ); 
      for(i = 0; i < n; i++) fprintf(stderr, "%g ", a[i]); 
      fprintf(stderr, "]\n");      
      break;
    }
    
    assert(isfinite(x));  

  } while(fabs(x - prev_x) > 1e-4 * mag); 
 

  return x;
}



/* Manage the empty clusters. Return the number of re-assigned clusters */
static int kmeans_reassign_empty (int d, int n, int k, float * centroids,
				  int * assign, int * nassign, unsigned int seed)
{
  int c, j, nreassign = 0;
  float * proba_split = fvec_new (k);
  float * vepsilon = fvec_new (d);


  for (c = 0 ; c < k ; c++)
    proba_split[c] = (nassign[c] < 2 ? 0 : nassign[c]*nassign[c] - 1);

  fvec_normalize (proba_split, k, 1);

  for (c = 0 ; c < k ; c++) {
    if(nassign[c]==0) {
      nreassign++;

      double rd = drand_r(&seed);
    
      /* j is the cluster that is selected for splitting */
      for (j = 0 ; j < k - 1 ; j++) {
	rd -= proba_split[j];
	if (rd < 0)
	  break;
      }
      fvec_cpy (centroids + c * d, centroids + j * d, d);
    
      /* generate a random perturbation vector, taking into account vector norm */
      double s = fvec_norm (centroids + j * d, d, 2) * 0.0000001;
      fvec_randn_r (vepsilon, d, rand_r(&seed));
      fvec_mul_by (vepsilon, d, s);
      fvec_add (centroids + j * d, vepsilon, d);
      fvec_sub (centroids + c * d, vepsilon, d);

      /* set the probability to choose cluster j to 0 for next re-assignment */
      proba_split[j] = 0;
      fvec_normalize (proba_split, k, 1);
    }
  }
  free (proba_split);
  free (vepsilon);

  return nreassign;
}


/* the core of the kmeans function (no initialization) */
static int kmeans_core (int d, int n, int k, int niter, int nt, int flags, int verbose, 
			 float * centroids, const float * v, 
                         unsigned int seed,
			 int * assign, int * nassign,
			 float * dis, 
			 double * qerr_out, long * iter_tot)
{
  int i, j, iter;
  int nreassign;
  int ret = 0; 

  /* the quantization error */
  double qerr = HUGE_VAL, qerr_old;

  float *tmp_v = NULL; 
  int *tmp_cumsum = NULL;
  
  if(flags & (KMEANS_L1 | KMEANS_CHI2)) {
    tmp_v = fvec_new((long)n * d);
    tmp_cumsum = ivec_new(k); 
  }                          

  int tot_nreassign=0;

  for (iter = 1 ; iter <= niter ; iter++) {
    (*iter_tot)++;

    /* Assign point to cluster and count the cluster size */

    knn_full_thread (flags & KMEANS_L1 ? 1 : 
                     flags & KMEANS_CHI2 ? 3 : 2, 
                     n, k, d, 1, centroids, v, NULL, assign, dis, nt);

    
    /* compute the number of points assigned to each cluster and a 
       probability to select a given cluster for splitting */
    ivec_0 (nassign, k);
    for (i = 0 ; i < n ; i++)
      nassign[assign[i]]++;

    if(flags & (KMEANS_L1 | KMEANS_CHI2)) {
      ivec_cpy(tmp_cumsum, nassign, k);
      ivec_cumsum(tmp_cumsum, k); 

      /* vector i will be written at column tmp_cumsum[assign[i]] - 1
         in tmp_v, a row major matrix  */
      
      for(i = 0; i < n; i++) {
        int write_to = --tmp_cumsum[assign[i]]; 
        for(j = 0; j < d; j++) 
          tmp_v[write_to + j * n] = v[i * d + j];        
      }
      
      /* centroids are medians */

      for(i = 0; i < k; i++) 
        for(j = 0; j < d; j++) 
          if(flags & KMEANS_L1) 
            centroids[i * d + j] = fvec_median(tmp_v + tmp_cumsum[i] + j * n, nassign[i]); 
          else /* if(flags & KMEANS_CHI2) */
            centroids[i * d + j] = minimize_sum_chi2(tmp_v + tmp_cumsum[i] + j * n, nassign[i]); 

    } else {

      /* update the centroids */
      fvec_0 (centroids, d * k);
      
      for (i = 0 ; i < n ; i++) {
        assert(assign[i] >= 0 || !"Something wrong in input. Maybe there are NaNs?");
        fvec_add (centroids + assign[i] * d, v + i * d, d);
      }
      
      /* normalize */
      for (i = 0 ; i < k ; i++) {          
        fvec_mul_by (centroids + i * d, d, 1.0 / nassign[i]);
      }
    } 

    if(flags & KMEANS_NORMALIZE_CENTS) 
      for (i = 0 ; i < k ; i++) 
        fvec_normalize(centroids + i * d, d, flags & KMEANS_L1 ? 1.0 : 2.0);

    /* manage empty clusters and update nassign */
    nreassign = kmeans_reassign_empty (d, n, k, centroids, assign, nassign, rand_r(&seed));
    if (nreassign > 0 && verbose)
      fprintf (stderr, "# kmeans warning: %d empty clusters -> split\n", nreassign);

    tot_nreassign+=nreassign;

    if(tot_nreassign>n/100 && tot_nreassign>1000) {
      fprintf (stderr,"# kmeans: reassigned %d times, abandoning\n", tot_nreassign);
      ret = -1;
      goto out; 
    }      

    /* compute the quantization error */
    qerr_old = qerr;
    qerr = fvec_sum (dis, n);

    if (qerr_old == qerr && nreassign == 0)
      break;
    if (verbose) {
      printf (" -> %.3f", qerr / n);
      fflush(stdout);
    }
  }
  if (verbose)      
    printf ("\n");

  *qerr_out = qerr;

 out:

  free(tmp_cumsum);
  free(tmp_v); 
  return ret;
}


float kmeans (int di, int n, int k, int niter, 
	      const float * v, int flags, long seed_in, int redo, 
	      float * centroids_out, float * dis_out, 
	      int * assign_out, int * nassign_out)
{
  long i, run, iter_tot = 0, d=di;

  int nt = flags & 0xffff;
  if (nt == 0) nt = 1;

  int verbose = !(flags & KMEANS_QUIET);

  niter = (niter == 0 ? 1000000 : niter);

  /* look at which variables have to be returned */
  int isout_centroids = (centroids_out == NULL ? 0 : 1);
  int isout_dis = (dis_out == NULL ? 0 : 1);
  int isout_assign = (assign_out == NULL ? 0 : 1);
  int isout_nassign = (nassign_out == NULL ? 0 : 1);

  /* if flags KMEANS_INIT_USER is activated, centroids_out contains initial centroids */
  int is_user_init = ((flags & KMEANS_INIT_USER) > 0 ? 1 : 0);
  if (is_user_init) {
    assert (centroids_out != NULL);
    redo = 1;   /* no randomness if initialization is provided by user */
  }

  /* the centroids */
  float * centroids = fvec_new (k * (size_t) d);
  
  /* store the distances from points to their nearest centroids */
  float * dis = fvec_new (n);
  
  /* the centroid indexes to which each vector is assigned */
  int * assign = ivec_new (n);
  
  /* the number of points assigned to a cluster */
  int * nassign = ivec_new (k);
  
  /* the total quantization error */
  double qerr, qerr_best = HUGE_VAL;
  
  /* for the initial configuration */
  int * selected = ivec_new (k);
  
  assert(k <= n || !"better to have fewer clusters than points");

  if (seed_in == 0) 
    seed_in = lrand48();

  unsigned int seed=seed_in;
  int core_ret=0;

  for (run = 0 ; run < redo ; run++) {
    if(verbose)
      printf ("<><><><> kmeans / run %d <><><><><>\n", (int)run);

    if (is_user_init) {
      fvec_cpy (centroids, centroids_out, d * k);
    } else {
      if (flags & KMEANS_INIT_BERKELEY) {
	int nsubset = n;
      
	if (n > k * 8 && n > 8192) { 
	  nsubset = k * 8; 
	  if(verbose) 
	    printf ("Restricting k-means++ initialization to %d points\n", nsubset);
	}
	kmeanspp_init (d, nsubset, k, v, selected, verbose, rand_r(&seed), nt);
      } else {
	random_init(d,n,k,v,selected,rand_r(&seed));
      }
      for (i = 0 ; i < k ; i++) 
	fvec_cpy (centroids + i * d, v + selected[i] * d, d);
    }

    core_ret=kmeans_core (d, n, k, niter, nt, flags, verbose, 
                          centroids, v, rand_r(&seed), assign, nassign, dis, 
                          &qerr, &iter_tot);
    
    if(core_ret<0) 
      break;


    /* If this run is the best encountered, save the results */
    if (qerr < qerr_best) {
      qerr_best = qerr;

      if (isout_centroids) 
	memcpy (centroids_out, centroids, k * d * sizeof (*centroids));
      if (isout_dis) 
	memcpy (dis_out, dis, n * sizeof (*dis));
      if (isout_assign) 
	memcpy (assign_out, assign, n * sizeof (*assign));
      if (isout_nassign)
	memcpy (nassign_out, nassign, k * sizeof (*nassign));
    }
  }

  if(verbose && core_ret>=0) {
    printf ("Total number of iterations: %d\n", (int)iter_tot);
    printf ("Unbalanced factor of last iteration: %g\n",ivec_unbalanced_factor(nassign,k));
  }
  
  /* free the variables that are not returned */
  free (selected);
  free (centroids);
  free (dis);
  free (assign);
  free (nassign);

  if(core_ret<0) 
    return -1;
  else 
    return qerr_best / n; 
}


/*---------- Functions for forward compatibility ----------*/

float *clustering_kmeans_assign_with_score (int n, int di,
                                            const float *points, int k, int nb_iter_max, 
                                            double normalize, 
                                            int n_thread,
                                            double *score, int **clust_assign_out)
{

  long d=di; /* to force 64-bit address computations */ 
  float *centroids = fvec_new (k * d);
  int *ca=clust_assign_out ? ivec_new(n) : NULL;
  int nredo=1;
  
  if(nb_iter_max/100000!=0) {
    nredo=nb_iter_max/100000;
    nb_iter_max=nb_iter_max % 100000;    
/*    printf("redo: %d iter: %d\n",nredo,nb_iter_max); */
  }   

  float ret=kmeans(di,n,k,nb_iter_max,points,n_thread | KMEANS_INIT_RANDOM,0,nredo,centroids,NULL,ca,NULL);

  if(ret>=0) {
    if(clust_assign_out) *clust_assign_out=ca;
    return centroids;
  } else {
    free(centroids);
    free(ca);
    *clust_assign_out=NULL;
    return NULL;
  }

}

float *clustering_kmeans_assign (int n, int d,
                                 const float *points, int k, int nb_iter_max,
                                 double normalize, int **clust_assign_out)
{
   return clustering_kmeans_assign_with_score (n, d, points, k, nb_iter_max,
                                               normalize, count_cpu(), NULL, clust_assign_out);
}

float *clustering_kmeans (int n, int d,
                          const float *points, int k, int nb_iter_max,
                          double normalize)
{

  int *clust_assign;

  float *centroids = clustering_kmeans_assign_with_score (n, d, points, k, 
		     nb_iter_max, normalize, count_cpu(), NULL, &clust_assign);
  free (clust_assign);

  return centroids;
}



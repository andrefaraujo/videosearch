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


#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>


#include "vector.h"
#include "matrix.h"
#include "kmeans.h"
#include "nn.h"
#include "gmm.h"
#include "sorting.h"
#include "machinedeps.h"

#include <sys/time.h>



/* Estimation of a Gaussian mixture (diagonal covariance matrix)
     k              number of mixture
     d              vector dimension
     g              gmm structure, namely:
     g->w   (k)     weights of the mixture
     g->mu  (k*d)   the centroids (mean of the mixture)
     g->sigma (k*d) the  diagonal of the covariance matrix
*/   


/* Initialize a new GMM structure */
static gmm_t * gmm_new (int d, int k)
{
  gmm_t * g = (gmm_t *) malloc (sizeof(*g));
  g->d=d;
  g->k=k;
  g->w = fvec_new (k);
  g->mu = fvec_new (k * d);
  g->sigma = fvec_new (k * d);

  return g;
}

/* Free an existing GMM structure */
void gmm_delete (gmm_t * g)
{
  free(g->w);
  free(g->mu);
  free(g->sigma);
  free(g);
}





#define real float
#define integer FINTEGER

int sgemm_ (char *transa, char *transb, integer * m, integer *
            n, integer * k, real * alpha, const real * a, integer * lda,
            const real * b, integer * ldb, real * beta, real * c__,
            integer * ldc);


int sgemv_(char *trans, integer *m, integer *n, real *alpha, 
           const real *a, integer *lda, const real *x, integer *incx, real *beta, real *y, 
           integer *incy);

#undef real
#undef integer

/* compute sum and diagonal of covariance matrix of a set of points (v) weighted by probabilities (p) */
static void compute_sum_dcov(int ni,int ki,int di,
                             const float *v,const float *mu_old,const float *p,
                             float *mu,float *sigma,float *w) {
  long i,j,l;
  FINTEGER n=ni,k=ki,d=di;

  for (j = 0 ; j < k ; j++) {
    double dtmp = 0;
    for (i = 0 ; i < n ; i++) dtmp += p[i * k + j];
    w[j] = dtmp;
  }

  float zero=0,one=1;
  sgemm_("Not transposed","Transposed",&d,&k,&n,&one,v,&d,p,&k,&zero,mu,&d);
  
  float *v2=fvec_new_cpy(v,n*(long)d);
  fvec_sqr(v2,n*(long)d);
  
  sgemm_("Not transposed","Transposed",&d,&k,&n,&one,v2,&d,p,&k,&zero,sigma,&d);
  free(v2);
  
  for (j = 0 ; j < k ; j++) {
    float *sigma_j=sigma+j*d;
    const float *mu_old_j=mu_old+j*d;
    const float *mu_j=mu+j*d;
    for(l=0;l<d;l++) {
      sigma_j[l]+=mu_old_j[l]*(mu_old_j[l]*w[j]-2*mu_j[l]);
    }
  }    

}

/* see threaded version below */
static void compute_sum_dcov_thread(int n,int k,int d,
                                    const float *v,const float *mu_old,const float *p,
                                    float *mu,float *sigma,float *w,
                                    int n_thread);


static float min_sigma=1e-10;

/* estimate the GMM parameters */
static void gmm_compute_params (int n, const float * v, const float * p, 
                                gmm_t * g,
                                int flags,                         
                                int n_thread)
{
  long i, j;

  long d=g->d, k=g->k;
  float * vtmp = fvec_new (d);
  float * mu_old = fvec_new_cpy (g->mu, k * d);
  float * w_old = fvec_new_cpy (g->w, k);

  fvec_0 (g->w, k);
  fvec_0 (g->mu, k * d);
  fvec_0 (g->sigma, k * d);

  if(0) {
    /* slow and simple */
    for (j = 0 ; j < k ; j++) {
      double dtmp = 0;
      for (i = 0 ; i < n ; i++) {
        /* contribution to the gaussian weight */
        dtmp += p[i * k + j];
        /* contribution to mu */
        
        fvec_cpy (vtmp, v + i * d, d);
        fvec_mul_by (vtmp, d, p[i * k + j]);
        fvec_add (g->mu + j * d, vtmp, d);
        
        /* contribution to the variance */
        fvec_cpy (vtmp, v + i * d, d);
        fvec_sub (vtmp, mu_old + j * d, d);
        fvec_sqr (vtmp, d);
        fvec_mul_by (vtmp, d, p[i * k + j]);
        fvec_add (g->sigma + j * d, vtmp, d);
        
      }
      g->w[j] = dtmp;
    }

  } else {
    /* fast and complicated */

    if(n_thread<=1) 
      compute_sum_dcov(n,k,d,v,mu_old,p,g->mu,g->sigma,g->w);
    else
      compute_sum_dcov_thread(n,k,d,v,mu_old,p,g->mu,g->sigma,g->w,n_thread);
  }

  if(flags & GMM_FLAGS_1SIGMA) {
    for (j = 0 ; j < k ; j++) {
      float *sigma_j=g->sigma+j*d;
      double var=fvec_sum(sigma_j,d)/d;
      fvec_set(sigma_j,d,var);
    }
  }

  long nz=0;
  for(i=0;i<k*d;i++) 
    if(g->sigma[i]<min_sigma) {
      g->sigma[i]=min_sigma;
      nz++;
    }

  if(nz) printf("WARN %ld sigma diagonals are too small (set to %g)\n",nz,min_sigma);

  for (j = 0 ; j < k ; j++) {
    fvec_div_by (g->mu + j * d, d, g->w[j]);
    fvec_div_by (g->sigma + j * d, d, g->w[j]);
  }

  assert(finite(fvec_sum(g->mu, k*d)));

  fvec_normalize (g->w, k, 1);

  printf ("w = ");
  fvec_print (g->w, k);
  double imfac = k * fvec_sum_sqr (g->w, k);
  printf (" imfac = %.3f\n", imfac);

  free (vtmp);
  free (w_old);
  free (mu_old);
}



double static inline sqr (double x)
{
  return x * x;
}


/* support function to compute log(a+b) from log(a) and log(b) 
   without loss of precision */
static double log_sum (double log_a, double log_b)
{
  double log_s;
  if (log_a < log_b)
    log_s = log_b + log (1 + exp (log_a - log_b));
  else
    log_s = log_a + log (1 + exp (log_b - log_a));
  return log_s;
}

#define CHECKFINITE(a) if(!finite(a)) {fprintf(stderr,"!!!! gmm_compute_p: not finite " #a "=%g at line %d\n",a,__LINE__); abort(); }; 


static void compute_mahalanobis_sqr(int n,long k,long d,
                                    const float *mu,const float *sigma,
                                    const float *v,
                                    float *p) {
  FINTEGER di=d,ki=k,ni=n; /* for blas functions */
  long i, j, l;
    
  float *mu2_sums=fvec_new(k);
  
  for (j = 0 ; j < k ; j++) {
    double dtmp = 0;
    for (l = 0 ; l < d ; l++) 
      dtmp += sqr(mu[j * d + l]) / sigma[j * d + l];      
    mu2_sums[j]=dtmp;
  }
  
  for (i = 0 ; i < n ; i++) 
    for (j = 0 ; j < k ; j++) 
      p[i * k + j]=mu2_sums[j];
  
  free(mu2_sums);
  
  float *v2=fvec_new(d*n);
  for (i = 0 ; i < n*d ; i++) 
    v2[i]=v[i]*v[i];
  
  float *inv_sigma=fvec_new(k*d);
  for (i = 0 ; i < k*d ; i++) 
    inv_sigma[i]=1.0/sigma[i];
  
  float one=1;
  
  sgemm_("Transposed","Not transposed",&ki,&ni,&di,&one,inv_sigma,&di,v2,&di,&one,p,&ki);
  
  free(v2);
  
  float *mu_sigma=inv_sigma;
  for (i = 0 ; i < k*d ; i++) 
    mu_sigma[i]=mu[i]/sigma[i];
  
  float minus_two=-2;
  
  sgemm_("Transposed","Not transposed",&ki,&ni,&di,&minus_two,mu_sigma,&di,v,&di,&one,p,&ki);  
  
  free(mu_sigma);      

}

/* compute p(ci|x). Warning: also update det */

void gmm_compute_p (int n, const float * v, 
                    const gmm_t * g, 
                    float * p,
                    int flags)
{
  if(n==0) return; /* sgemm doesn't like empty matrices */

  long i, j, l;
  double dtmp;
  long d=g->d, k=g->k;


  float * logdetnr = fvec_new(k);

  for (j = 0 ; j < k ; j++) {
    logdetnr[j] = -d / 2.0 * log (2 * M_PI);
    for (i = 0 ; i < d ; i++)
      logdetnr[j] -= 0.5 * log (g->sigma[j * d + i]);
  }

  /* compute all probabilities in log domain */

  /* compute squared Mahalanobis distances (result in p) */

  if(0) { /* simple & slow */
    for (i = 0 ; i < n ; i++) {
      for (j = 0 ; j < k ; j++) {
        dtmp = 0;
        for (l = 0 ; l < d ; l++) {
          dtmp += sqr (v[i * d + l] - g->mu[j * d + l]) / g->sigma[j * d + l];
        }
        p[i * k + j] = dtmp;
      }
    }
  } else { /* complicated & fast */
    compute_mahalanobis_sqr(n,k,d,g->mu,g->sigma,v,p); 
  }

  /* convert distances to probabilities, staying in the log domain 
     until the very end */
  for (i = 0 ; i < n ; i++) {

    for (j = 0 ; j < k ; j++) {
      p[i * k + j] = logdetnr[j] - 0.5 * p[i * k + j];
      CHECKFINITE(p[i * k + j]);
    }

    /* at this point, we have p(x|ci) -> we want p(ci|x) */
    

    if(flags & GMM_FLAGS_NO_NORM) {     /* compute the normalization factor */

      dtmp=0;

    } else {

      dtmp = p[i * k + 0];
      
      if(flags & GMM_FLAGS_W) 
        dtmp+=log(g->w[0]);

      for (j = 1 ; j < k ; j++) {
        double log_p=p[i * k + j];

        if(flags & GMM_FLAGS_W) 
          log_p+=log(g->w[j]);

        dtmp = log_sum (dtmp, log_p);
      }

      /* now dtmp contains the log of sums */
    } 

    for (j = 0 ; j < k ; j++) {
      double log_norm=0;

      if(flags & GMM_FLAGS_W)
        log_norm=log(g->w[j])-dtmp;
      else
        log_norm=-dtmp;

      p[i * k + j] = exp (p[i * k + j] + log_norm);
      CHECKFINITE(p[i * k + j]);
    }

    //    printf ("p[%d] = ", i);
    //    fvec_print (p + i * k, k);
  }

  free(logdetnr);

}



void gmm_handle_empty(int n, const float *v, gmm_t *g, float *p) {
  long d=g->d, k=g->k;
  
  long nz=fvec_count_occurrences(p,k*n,0);
  printf("nb of 0 probabilities: %ld / (%ld*%d) = %.1f %%\n",
         nz,k,n,nz*100.0/(k*n));         

  int i,j;
  float *w=fvec_new_0(k);
  for (i = 0 ; i < n ; i++) 
    for (j = 0 ; j < k ; j++) 
      w[j]+=p[j+i*k];
      
  int bigprime=1000003;

  for (j = 0 ; j < k ; j++) if(w[j]==0) {
    printf("center %d is empty....",j);
    fflush(stdout);
    int j2;

    j2=j;
    for(i=0;i<k;i++) {
      j2=(j2+bigprime)%k; 
      if(w[j2]>0) break;
    }
    assert(i<k || !"could not find centroid to split, veeeery bad input data");
    
    /* dimension to split: that with highest variance */
    int split_dim = fvec_arg_max (g->sigma + d * j2, d);

    /* transfer half(?) of the points from j2 -> j */
    int nt=0,nnz=0;
    for(i=0;i<n;i++) if(p[j2+i*k]>0) { 
      nnz++;
      if(v[i*d+split_dim]<g->mu[j2*d+split_dim]) {
        p[j+i*k]=p[j2+i*k];
        p[j2+i*k]=0;
        nt++;
      }
    }

    printf("split %d at dim %d (variance %g, transferred %d/%d pts)\n",                      
           j2,split_dim,g->sigma[d*j2+split_dim],nt,nnz);    

    w[j2]=-1; /* avoid further splits */
  }
  
  free(w);

}


gmm_t * gmm_learn (int di, int ni, int ki, int niter, 
	   const float * v, int nt, int seed, int nredo,
	   int flags)
{
  long d=di,k=ki,n=ni;

  int iter, iter_tot = 0;
  double old_key, key = 666;

  niter = (niter == 0 ? 10000 : niter);

  /* the GMM parameters */
  float * p = fvec_new_0 (n * k);      /* p(ci|x) for all i */
  gmm_t * g = gmm_new (d, k);

  /* initialize the GMM: k-means + variance estimation */
  int * nassign = ivec_new (n);  /* not useful -> to be removed when debugged */
  float * dis = fvec_new (n);
  kmeans (d, n, k, niter, v, nt, seed, nredo, g->mu, dis, NULL, nassign); 
  
  fflush (stderr);
  fprintf (stderr, "assign = ");
  ivec_print (nassign, k);
  fprintf (stderr, "\n");
  free (nassign);

  /* initialization of the GMM parameters assuming a diagonal matrix */
  fvec_set (g->w, k, 1.0 / k);
  double sig = fvec_sum (dis, n) / n;
  printf ("sigma at initialization = %.3f\n", sig);
  fvec_set (g->sigma, k * d, sig);
  free (dis);


  /* start the EM algorithm */
  fprintf (stdout, "<><><><> GMM  <><><><><>\n");
      
  if(flags & GMM_FLAGS_PURE_KMEANS) niter=0;

  for (iter = 1 ; iter <= niter ; iter++) {
    
    gmm_compute_p_thread (n, v, g, p, flags, nt);
    fflush(stdout);

    gmm_handle_empty(n, v, g, p);
    
    gmm_compute_params (n, v, p, g, flags, nt);
    fflush(stdout);


    iter_tot++;

    /* convergence reached -> leave */
    old_key = key;
    key = fvec_sum (g->mu, k * d);

    printf ("keys %5d: %.6f -> %.6f\n", iter, old_key, key);
    fflush(stdout);

    if (key == old_key)
      break;
  }
  fprintf (stderr, "\n");

  free(p);
  
  return g;
}

size_t gmm_fisher_sizeof(const gmm_t * g,int flags) {
  int sz=0;
  if(flags & GMM_FLAGS_W) sz+=g->k-1;
  if(flags & GMM_FLAGS_MU) sz+=g->d*g->k;
  if(flags & GMM_FLAGS_1SIGMA) sz+=g->k;
  if(flags & GMM_FLAGS_SIGMA) sz+=g->d*g->k;
  return sz;
}

void gmm_fisher(int n, const float *v, const gmm_t * g, int flags, float *dp_dlambda) {
    float *word_total_soft_assignment = NULL;
    gmm_fisher_save_soft_assgn(n, v, g, flags, dp_dlambda, 
                               word_total_soft_assignment);
}


void gmm_fisher_save_soft_assgn(int n, const float *v, const gmm_t * g, int flags,
                                float *dp_dlambda,
                                float *word_total_soft_assignment) {
  long d=g->d, k=g->k;
  float *p = fvec_new(n * k);
  long i,j,l;
  long ii=0;

  float * vp = NULL; /* v*p */
  float * sum_pj = NULL; /* sum of p's for a given j */  

  gmm_compute_p(n,v,g,p,flags | GMM_FLAGS_W);

#define P(j,i) p[(i)*k+(j)]
#define V(l,i) v[(i)*d+(l)]
#define MU(l,j) g->mu[(j)*d+(l)]
#define SIGMA(l,j) g->sigma[(j)*d+(l)]
#define VP(l,j) vp[(j)*d+(l)]

  // Save total soft assignment per centroid
  if (word_total_soft_assignment != NULL) {
      for (j=0;j<k;j++) {
          double sum=0;
          for (i=0;i<n;i++) {
              sum += P(j,i);
          }
          if (n != 0) {
              word_total_soft_assignment[j] = (float)(sum/n);
          } else {
              word_total_soft_assignment[j] = 0.0;
          }
      }
  }

  if(flags & GMM_FLAGS_W) {

    for(j=1;j<k;j++) {
      double accu=0;
      
      for(i=0;i<n;i++) 
        accu+= P(j,i)/g->w[j] - P(0,i)/g->w[0];
      
      /* normalization */
      double f=n*(1/g->w[j]+1/g->w[0]);
      
      dp_dlambda[ii++]=accu/sqrt(f);
    }
  } 
  
  if(flags & GMM_FLAGS_MU) {
    float *dp_dmu=dp_dlambda+ii;

#define DP_DMU(l,j) dp_dmu[(j)*d+(l)]
    
    if(0) { /* simple and slow */
    
      for(j=0;j<k;j++) {
        for(l=0;l<d;l++) {
          double accu=0;
          
          for(i=0;i<n;i++) 
            accu += P(j,i) * (V(l,i)-MU(l,j)) / SIGMA(l,j);
          
          DP_DMU(l,j)=accu;
        }
      }
      
    } else { /* complicated and fast */

      /* precompute  tables that may be useful for sigma too */
      vp = fvec_new(k * d);
      fmat_mul_tr(v,p,d,k,n,vp);

      sum_pj = fvec_new(k);
      for(j=0;j<k;j++) {        
        double sum=0;        
        for(i=0;i<n;i++) sum += P(j,i);        
        sum_pj[j] = sum;
      }

      for(j=0;j<k;j++) {
        for(l=0;l<d;l++)
          DP_DMU(l,j) = (VP(l,j) - MU(l,j) * sum_pj[j]) / SIGMA(l,j);
      }

    }
    /* normalization */
    if(!(flags & GMM_FLAGS_NO_NORM)) {
      for(j=0;j<k;j++) 
        for(l=0;l<d;l++) {
          float nf = sqrt(n*g->w[j]/SIGMA(l,j));
          if(nf > 0) DP_DMU(l,j) /= nf;                
        }        
    }
#undef DP_DMU
    ii+=d*k;
  }

  if(flags & (GMM_FLAGS_SIGMA | GMM_FLAGS_1SIGMA)) {

    
    if(flags & GMM_FLAGS_1SIGMA) { /* fast not implemented for 1 sigma */

      for(j=0;j<k;j++) {
        double accu2=0;
        for(l=0;l<d;l++) {
          double accu=0;
        
          for(i=0;i<n;i++) 
            accu += P(j,i) * (sqr(V(l,i)-MU(l,j)) / SIGMA(l,j) - 1) / sqrt(SIGMA(l,j));
        
          if(flags & GMM_FLAGS_SIGMA) {

            double f=flags & GMM_FLAGS_NO_NORM ? 1.0 : 2*n*g->w[j]/SIGMA(l,j);
          
            dp_dlambda[ii++]=accu/sqrt(f);
          } 
          accu2+=accu;        
        }

        if(flags & GMM_FLAGS_1SIGMA) {
          double f=flags & GMM_FLAGS_NO_NORM ? 1.0 : 2*d*n*g->w[j]/SIGMA(0,j);
          dp_dlambda[ii++]=accu2/sqrt(f);        
        }

      }  
    
    } else { /* fast and complicated */
      assert(flags & GMM_FLAGS_SIGMA);
      float *dp_dsigma = dp_dlambda + ii;

      if(!vp) {
        vp = fvec_new(k * d);
        fmat_mul_tr(v,p,d,k,n,vp);
      }

      if(!sum_pj) {
        sum_pj = fvec_new(k);
        for(j=0;j<k;j++) {        
          double sum=0;        
          for(i=0;i<n;i++) sum += P(j,i);        
          sum_pj[j] = sum;
        }
      }
      float *v2 = fvec_new(n * d);
      for(i = n*d-1 ; i >= 0; i--) v2[i] = v[i] * v[i];
      float *v2p = fvec_new(k * d);
      fmat_mul_tr(v2,p,d,k,n,v2p);
      free(v2);

#define V2P(l,j) v2p[(j)*d+(l)]
#define DP_DSIGMA(i,j) dp_dsigma[(i)+(j)*d]
      for(j=0;j<k;j++) {

        for(l=0;l<d;l++) {
          double accu;

          accu = V2P(l, j);

          accu += VP(l, j) * (- 2 * MU(l,j));

          accu += sum_pj[j] * (sqr(MU(l,j))  - SIGMA(l,j));

          /* normalization */

          double f;

          if(flags & GMM_FLAGS_NO_NORM) {
            f = pow(SIGMA(l,j), -1.5);
          } else {
            f = 1 / (SIGMA(l,j) * sqrt(2*n*g->w[j]));
          }

          DP_DSIGMA(l,j) = accu * f;

        }

      }  
      
      free(v2p);

#undef DP_DSIGMA
#undef V2P
      ii += d * k;
    }

  }
  
  assert(ii==gmm_fisher_sizeof(g,flags));
#undef P
#undef V
#undef MU
#undef SIGMA
  free(p);
  free(sum_pj);
  free(vp);
}

void gmm_fisher_point_indexed(int n, const float *v, const gmm_t * g,
                              int flags, unsigned int* pi_assgns,
                              float* pi_assgn_weights,
                              float* pi_residuals) {
  long d=g->d, k=g->k;
  float *p = fvec_new(n * k);
  long i,j,l;

  gmm_compute_p(n,v,g,p,flags | GMM_FLAGS_W);

#define P(j,i) p[(i)*k+(j)]
#define V(l,i) v[(i)*d+(l)]
#define MU(l,j) g->mu[(j)*d+(l)]
#define SIGMA(l,j) g->sigma[(j)*d+(l)]

  /*Loop over input vectors and populate information to be returned*/
  for (i=0;i<n;i++) {
    /*Find Gaussian with largest soft assignment*/
    int largest_ind = 0;
    float largest_p = 0;
    for (j=0;j<k;j++) {
        if (P(j,i) > largest_p) {
            largest_ind = j;
            largest_p = P(j,i);
        }
    }

    /*Set pi_assgns for this item*/
    pi_assgns[i] = largest_ind;

    /*Set pi_assgn_weights for this item*/
    pi_assgn_weights[i] = largest_p/sqrt(g->w[largest_ind]);

    /*Set pi_residuals for this item*/
    for (l=0;l<d;l++) {
        pi_residuals[i*d + l] = (V(l,i) - MU(l,largest_ind))/SIGMA(l,largest_ind);
    }

  }

  free(p);
}

void gmm_print(const gmm_t *g) {
  printf("gmm (%d gaussians in %d dim)=[\n",g->k,g->d);
  int i,j;
  for(i=0;i<g->k;i++) {
    printf("   w=%g, mu=[",g->w[i]);
    for(j=0;j<g->d;j++) printf("%g ",g->mu[i*g->d+j]);
    printf("], sigma=diag([");
    for(j=0;j<g->d;j++) printf("%g ",g->sigma[i*g->d+j]);    
    printf("])\n");
  }
  printf("]\n");
}

#define WRITEANDCHECK(a,n) if(fwrite(a,sizeof(*a),n,f)!=n) {perror("gmm_write"); abort(); }


void gmm_write(const gmm_t *g, FILE *f) {
  
  WRITEANDCHECK(&g->d,1);
  WRITEANDCHECK(&g->k,1);
  WRITEANDCHECK(g->w,g->k);
  WRITEANDCHECK(g->mu,g->k*g->d);
  WRITEANDCHECK(g->sigma,g->k*g->d);
  
}

#define READANDCHECK(a,n) if(fread(a,sizeof(*a),n,f)!=n) {perror("gmm_read"); abort(); }


gmm_t * gmm_read(FILE *f) {
  int d,k;

  READANDCHECK(&d,1);
  READANDCHECK(&k,1);

  gmm_t *g=gmm_new(d,k);
  
  READANDCHECK(g->w,g->k);
  READANDCHECK(g->mu,g->k*g->d);
  READANDCHECK(g->sigma,g->k*g->d);
  
  return g;
}



/*********************** threaded versions **********************/



typedef struct {
  long n;
  const float * v;
  const gmm_t * g;
  float * p;
  int do_norm;
  int n_thread;   
} compute_p_params_t;

/* n sliced */
static void compute_p_task_fun (void *arg, int tid, int i) {
  compute_p_params_t *t=arg;
  long n0=i*t->n/t->n_thread;
  long n1=(i+1)*t->n/t->n_thread;
  
  gmm_compute_p(n1-n0, t->v+t->g->d*n0, t->g, t->p+t->g->k*n0, t->do_norm);
}

void gmm_compute_p_thread (int n, const float * v, 
                           const gmm_t * g, 
                           float * p, 
                           int do_norm,
                           int n_thread) {
  compute_p_params_t t={n,v,g,p,do_norm,n_thread};
  compute_tasks(n_thread,n_thread,&compute_p_task_fun,&t);
}


typedef struct {
  long n,k,d;
  const float *v,*mu_old,*p;
  float *mu,*sigma,*w;
  int n_thread;
} compute_sum_dcov_t;

/* n sliced */
static void compute_sum_dcov_task_fun (void *arg, int tid, int i) {
  compute_sum_dcov_t *t=arg;
  long n0=i*t->n/t->n_thread;
  long n1=(i+1)*t->n/t->n_thread;
  
  compute_sum_dcov(n1-n0,t->k,t->d,t->v+t->d*n0,
                   t->mu_old,t->p+n0*t->k,
                   t->mu+i*t->d*t->k,
                   t->sigma+i*t->d*t->k,
                   t->w+t->k*i);

}



static void compute_sum_dcov_thread(int ni,int ki,int di,
                                    const float *v,const float *mu_old,const float *p,
                                    float *mu,float *sigma,float *w,
                                    int n_thread) {
  long n=ni,d=di,k=ki;

  compute_sum_dcov_t t={
    n,k,d,
    v,mu_old,p,
    fvec_new(n_thread*d*k), /* mu */
    fvec_new(n_thread*d*k), /* sigma */
    fvec_new(n_thread*k), /* w */
    n_thread
  };

  compute_tasks(n_thread,n_thread,&compute_sum_dcov_task_fun,&t);
  
  /* accumulate over n's */

  long i;
  fvec_cpy(mu,t.mu,k*d);
  fvec_cpy(sigma,t.sigma,k*d);
  fvec_cpy(w,t.w,k);
  for(i=1;i<n_thread;i++) {
    fvec_add(mu,t.mu+i*d*k,d*k);
    fvec_add(sigma,t.sigma+i*d*k,d*k);    
    fvec_add(w,t.w+i*k,k);    
  }
  free(t.mu);
  free(t.sigma);
  free(t.w);
}





/*
Copyright Â© INRIA 2009-2014.
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
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include "machinedeps.h"
#include "vector.h"
#include "nn.h"
#include "binheap.h"
#include "sorting.h"



#define NEWA(type,n) (type*)malloc((n)*sizeof(type))
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))


/**********************************************************************************
 * Distance functions  
 */


/*------------------ Blas subroutine ------------------*/

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


/*
 * computes dist2 := dist2 - 2 * descs * clusters' 
 * where 
 *   dist2    is ndesc-by-nclust 
 *   clusters is nclust-by-d
 *   descs    is ndesc-by-d
 * (all matrices stored by lines, à la C, and packed)
 */


static void add_matmul (FINTEGER d, FINTEGER na, FINTEGER nb,
                        const float *a, FINTEGER lda, 
                        const float *b, FINTEGER ldb,
                        float *dist2, FINTEGER ldd)
{
  /* ldd >= na */

  float minus_two = -2;
  float one = 1;

  sgemm_ ("Transposed", "Not trans", &na, &nb, &d,
          &minus_two, a, &lda, b, &ldb, &one, dist2, &ldd);

}


static void add_matvecmul (FINTEGER d, FINTEGER nb,
                           const float *a, 
                           const float *b, FINTEGER ldb,
                           float *dist2)
{
  /* ldd >= na */

  float minus_two = -2;
  float one = 1;                               
  FINTEGER ione = 1;

  sgemv_ ("Transposed", &d, &nb, &minus_two, b, &ldb, a, &ione, &one, dist2, &ione);

}



/* computes all distances between a line of a and a line of b. 
 *   a(na,d) by lines
 *   b(nb,d) by lines
 *  dist2[i+na*j] = || a(i,:)-b(j,:) ||^2
 */
void compute_cross_distances (int d, int na, int nb,
                              const float *a, const float *b,  
                              float *dist2) 
{
  compute_cross_distances_nonpacked (d, na, nb, a, d, b, d, dist2, na);
}


void compute_cross_distances_nonpacked (int d, int na, int nb,
                                        const float *a, int lda,
                                        const float *b, int ldb, 
                                        float *dist2, int ldd)
{
  long i, j;
  float *sum_c2 = (float *) malloc (sizeof (float) * na);

  for (i = 0; i < na; i++) {
    float s = 0;
    const float *cl = a + lda * i;
    for (j = 0; j < d; j++)
      s += cl[j] * cl[j];
    sum_c2[i] = s;
  }

  for (i = 0; i < nb; i++) {
    double sum_d2 = 0;
    const float *dl = b + ldb * i;
    for (j = 0; j < d; j++)
      sum_d2 += dl[j] * dl[j];
    float *d2l = dist2 + i * ldd;
    for (j = 0; j < na; j++)
      d2l[j] = sum_d2 + sum_c2[j];
  }

  add_matmul (d, na, nb, a, lda, b, ldb, dist2, ldd);

  free (sum_c2);
}


void compute_distances_1_nonpacked (int d, int nb,
				    const float *a, 
				    const float *b, int ldb, 
				    float *dist2)
{
  long i, j;
  double sum_c2;

  sum_c2=0;
  for (j = 0; j < d; j++)
    sum_c2 += a[j] * a[j];

  for (i = 0; i < nb; i++) {
    double sum_d2 = 0;
    const float *dl = b + ldb * i;
    for (j = 0; j < d; j++)
      sum_d2 += dl[j] * dl[j];
    dist2[i] = sum_d2 + sum_c2;
  }
 
  add_matvecmul (d, nb, a, b, ldb, dist2);

}

void compute_distances_1 (int d, int nb,
                          const float *a, 
                          const float *b,                         
                          float *dist2) 
{
  compute_distances_1_nonpacked(d,nb,a,b,d,dist2);  
}






#ifndef __SSE2__

#warning "SSE optimized distance computations not set"

#else

#include <emmintrin.h>

/* compute chi2 distance between two vectors */
static float vec_chi2 (const float *a, const float *b, int n) 
{
  int i=0;
  __v4sf accu={0,0,0,0};
  __v4sf *av=(void*)a,*bv=(void*)b;
  n/=4;

  for (i = 0 ; i < n ; i++) {
    __v4sf ai = av[i],bi=bv[i];
    __v4sf sum = ai+bi,diff=ai-bi;
    __v4sf mask = _mm_cmpneq_ps(ai,bi);
    accu+=_mm_and_ps (_mm_div_ps(_mm_mul_ps(diff,diff),sum),mask);
  }

  float *af=(void*)&accu;
      
  return af[0] + af[1] + af[2] + af[3];
}


static float vec_L1 (const float *a, const float *b, int n) 
{
  int i=0;
  __v4sf accu = {0,0,0,0};
  __v4si signbit = {0x7fffffff, 0x7fffffff, 0x7fffffff, 0x7fffffff};  /* to clear out sign bit */
  __v4sf *av=(void*)a,*bv=(void*)b;
  n/=4;

  for(i=0;i<n;i++) {
    __v4sf ai = av[i], bi = bv[i];
    __v4sf diff = ai - bi;
    __v4sf diffabs = _mm_and_ps(diff, (__m128)signbit);
    accu += diffabs;
  }

  float *af=(void*)&accu;
      
  return af[0]+af[1]+af[2]+af[3];
}


/* TODO optimize a bit more with blocks */
static void cross_distances_chi2_vec (int d, int na, int nb,
				      const float *a, int lda,
				      const float *b, int ldb,
				      float *c, int ldc) 
{  
  int i, j;
  const float *bl = b;
  float *cl=c;    

  for(j=0;j<nb;j++) {
    const float *al=a;
    for(i=0;i<na;i++) {
      cl[i]=vec_chi2(al,bl,d);;
      al+=lda;
    }
    cl+=ldc;
    bl+=ldb;
  }
}

static void cross_distances_L1_vec(int d,int na,int nb,
                              const float *a,int lda,
                              const float *b,int ldb,
                              float *c,int ldc) {  
  int i,j;
  const float *bl=b;
  float *cl=c;    

  for(j=0;j<nb;j++) {
    const float *al=a;
    for(i=0;i<na;i++) {
      cl[i]=vec_L1(al,bl,d);;
      al+=lda;
    }
    cl+=ldc;
    bl+=ldb;
  }
}

#endif  

static double sqr (double x)
{
  return x * x;
}

static void mat_product(FINTEGER d, FINTEGER na, FINTEGER nb,
                        const float *a, FINTEGER lda,
                        const float *b, FINTEGER ldb,
                        float *dist2, FINTEGER ldd) {
  float zero = 0; 
  float one = 1;

  sgemm_("Trans", "Not", &na, &nb, &d, &one, a, &lda, b, &ldb, &zero, dist2, &ldd); 

}


/* alternative distance functions */

void compute_cross_distances_alt_nonpacked (int distance_type, int d, int na, int nb,
                                            const float *a, int lda,
                                            const float *b, int ldb,
                                            float *dist2, int ldd) {
  /* special cases for optimized versions */
#ifdef __SSE2__
  if(d%4==0) {
    if(distance_type == 3) {
      cross_distances_chi2_vec(d,na,nb,a,lda,b,ldb,dist2,ldd);
      return;
    } else if(distance_type == 1) {
      cross_distances_L1_vec(d,na,nb,a,lda,b,ldb,dist2,ldd);
      return;
    }
  }
#endif  


  if(distance_type == 12) {
    compute_cross_distances_nonpacked(d, na, nb, a, lda, b, ldb, dist2, ldd);
    return;
  } 

  if(distance_type == 16) {
    mat_product(d, na, nb, a, lda, b, ldb, dist2, ldd);
    return;
  }

  int i,j,k;

  for(j=0;j<nb;j++) {
    float * dline=dist2+j*ldd;
    const float * bline=b+ldb*j;
    const float * aline=a;
    for(i=0;i<na;i++) {

      double sum=0;
      if(distance_type==1) 
        for(k=0;k<d;k++) sum+=fabs(aline[k]-bline[k]);
      else if(distance_type==2) 
        for(k=0;k<d;k++) sum+=sqr(aline[k]-bline[k]);
      else if(distance_type==3) 
        for(k=0;k<d;k++) {
          float av=aline[k],bv=bline[k];
          sum+=av+bv==0 ? 0 : sqr(av-bv)/(av+bv);
        }      
      else if(distance_type==4)
        for(k=0;k<d;k++) {
          float av=aline[k],bv=bline[k];
          float den=fabs(av+bv);
          sum+=den==0 ? 0 : sqr(av-bv)/den;
        }
      else if(distance_type == 5) {
        for(k=0;k<d;k++) {
          float av = aline[k], bv = bline[k];
          sum += MIN(av, bv); 
        }        
      }
      else if(distance_type == 6) {
        for(k=0;k<d;k++) {
          float av = aline[k], bv = bline[k];
          sum += av * bv; 
        }        
      }
      dline[i]=sum;

      aline+=lda;
    }
  }

}

void compute_cross_distances_alt (int distance_type, int d, int na, int nb,
                                  const float *a,
                                  const float *b, float *dist2) {
  compute_cross_distances_alt_nonpacked(distance_type, d, na, nb, a, d, b, d, dist2, na);
}




/**********************************************************************************
 * Elementary cluster assignment 
 */


/*
 * Computations are done by blocks (nice for cache access, distance matrix must fit in mem)
 * blocks are BLOCK_N1 * BLOCK_N2
 */

#define BLOCK_N1 256
#define BLOCK_N2 256



#define MIN(a,b) ((a)<(b) ? (a) : (b))

/* This function quantizes the vectors coords according to the codebook clusters, 
   and sets the corresponding indexes in the index vector vw accordingly    
 */

/* n1 = pts */
static void nn_single_full (int distance_type,
			    int n1, int n2, int d,
			    const float *mat2, const float *mat1, 
			    const float *vw_weights,                             
			    int *vw, float *vwdis)
{
  int step1 = MIN (n1, BLOCK_N1), step2 = MIN (n2, BLOCK_N2);

  float *dists = fvec_new (step1 * step2);

  /* divide the dataset into sub-blocks to:
   * - not make a too big dists2 output array 
   */
  
  long i1,i2,j1,j2;
  for (i1 = 0; i1 < n1; i1 += step1) {  

    int m1 = MIN (step1, n1 - i1);

    /* clear mins */

    for (j1 = 0; j1 < m1; j1++) {
      vw[j1+i1]=-1;
      vwdis[j1+i1]=1e30;
    }

    for (i2 = 0; i2 < n2 ; i2 += step2) {     
      
      int m2 = MIN (step2, n2 - i2);
      
      if(distance_type==2)       
        compute_cross_distances (d, m2, m1, mat2+i2*d, mat1+i1*d, dists);
      else
        compute_cross_distances_alt (distance_type, d, m2, m1, mat2+i2*d, mat1+i1*d, dists);

      if(vw_weights) {
        for(j1=0;j1<m1;j1++) for (j2 = 0; j2 < m2; j2++)
          dists[j1 * m2 + j2] *= vw_weights[j2 + i2];        
      }

      /* update mins */

      for(j1=0;j1<m1;j1++) {
        float *dline=dists+j1*m2;
        
        int imin=vw[i1+j1];
        float dmin=vwdis[i1+j1];

        for(j2=0;j2<m2;j2++) 
          if(dline[j2]<dmin) {
            imin=j2+i2;
            dmin=dline[j2];
          }
          
        vw[i1+j1]=imin;
        vwdis[i1+j1]=dmin;

      }      

    }  
  }

  free (dists);
}




void knn_full (int distance_type,int n1, int n2, int d, int k,
	       const float *mat2, const float *mat1,
	       const float *vw_weights,
	       int *vw, float *vwdis)
{
  assert (k <= n2);

  if(k==1) {
    nn_single_full(distance_type, n1, n2, d, mat2, mat1, vw_weights, vw, vwdis);
    return;
  }

  
  int step1 = MIN (n1, BLOCK_N1), step2 = MIN (n2, BLOCK_N2);

  float *dists = fvec_new (step1 * step2);


  /* allocate all heaps at once */
  long oneh = fbinheap_sizeof(k);
  /* oneh=(oneh+7) & ~7; */  /* round up to 8 bytes */

  char *minbuf = malloc (oneh * step1);

#define MINS(i) ((fbinheap_t*)(minbuf + oneh * i))
  
  long i1,i2,j1,j2;
  for (i1 = 0; i1 < n1; i1 += step1) {  

    int m1 = MIN (step1, n1 - i1);

    /* clear mins */
    for (j1 = 0; j1 < m1; j1++) 
      fbinheap_init(MINS(j1),k);
        

    for (i2 = 0; i2 < n2 ; i2 += step2) {     
      
      int m2 = MIN (step2, n2 - i2);
      
      
      if(distance_type==2)       
        compute_cross_distances (d, m2, m1, mat2+i2*d, mat1+i1*d, dists);
      else 
        compute_cross_distances_alt (distance_type, d, m2, m1, mat2+i2*d, mat1+i1*d, dists);    

      if(vw_weights) {
        for(j1=0;j1<m1;j1++) for (j2 = 0; j2 < m2; j2++)
          dists[j1 * m2 + j2] *= vw_weights[j2 + i2];        
      }

      /* update mins */

      for(j1=0;j1<m1;j1++) {
        float *dline=dists+j1*m2; 
        fbinheap_addn_label_range(MINS(j1),m2,i2,dline);
      }      

    }  

    for (j1 = 0; j1 < m1; j1++) {
      fbinheap_t *mh = MINS(j1);
      /* mh->k may not be same as k if there are NaNs in the distances */
      fbinheap_sort(mh, vw + (i1+j1) * k, vwdis + (i1+j1) * k);
      if(mh->k < k) {
        memset(vw + (i1+j1) * k + mh->k, 0xff, sizeof(int) * (k - mh->k));
        memset(vwdis + (i1+j1) * k + mh->k, 0xff, sizeof(float) * (k - mh->k));
      }
    }
  }

#undef MINS
  free (minbuf);
  free(dists);
}


void knn_reorder_shortlist(int n, int nb, int d, int k,
                           const float *b, const float *v,
                           int *assign,
                           float *dists) 
{
  float *subb=fvec_new(k*d);
  float *diststmp=fvec_new(k);
  int *perm=ivec_new(k);
  int *assigntmp=ivec_new(k);  
  int i,j;

  for(i=0;i<n;i++) {
    int *assigni=assign+i*k;
    float *disti=dists+i*k;

    int ki  ;
    if(1) {

      for(j=0;j<k;j++) {
        if(assigni[j]<0) break;
        memcpy(subb+j*d,b+assigni[j]*(long) d,sizeof(*subb)*d);
      }

      ki=j;

    } else {
      for(j=0;j<k;j++) 
        if(assigni[j]<0) break;
      ki=j;
      ivec_sort(assigni,ki); /* to improve access locality */
      for(j=0;j<ki;j++) {
        memcpy(subb+j*d,b+assigni[j]*(long) d,sizeof(*subb)*d);
      } 
    }


    compute_distances_1(d,ki,v+i*d,subb,diststmp);
    
    fvec_sort_index(diststmp,ki,perm);

    memcpy(assigntmp,assigni,sizeof(*assigni)*ki);
    
    for(j=0;j<ki;j++) {
      disti[j]=diststmp[perm[j]];
      assigni[j]=assigntmp[perm[j]];
    }
    
  }
  free(assigntmp);
  free(diststmp);
  free(subb);
  free(perm);
}


void knn_recompute_exact_dists(int nq, int nb, int d, int k,
			       const float *b, const float *v,
			       int label0, int *kp,
			       const int *idx, float *dis) 
{
  long q, i;

  for(q = 0; q < nq; q++) {
    const float * vq = v + d * q;
    for(i = kp[q]; i < k; i++) {
      long j = idx[q * k + i] - label0;
      assert(j >= 0); 
      if(j >= nb) break;      
      dis[q * k + i] = fvec_distance_L2sqr(vq, b + j * d, d);
    }
    kp[q] = i;
  }
}


/**********************************************************************************
 * Simple call versions
 */


double nn (int npt, int nclust, int d,
	 const float *codebook, const float *coords, int *vw) 
{
  
  /* The distances to centroids that will be returned */
  float *vwdis = fvec_new(npt);
  
  knn_full (2, npt, nclust, d, 1, codebook, coords, NULL, vw, vwdis);
  
  double toterr = fvec_sum(vwdis, npt);
  free(vwdis);

  return toterr;
}


float *knn (int npt, int nclust, int d, int k,
	    const float *codebook, const float *coords, int *vw)
{
  /* The distances to centroids that will be returned */
  float *vwdis = fvec_new(npt * k);

  knn_full (2, npt, nclust, d, k, codebook, coords, NULL, vw, vwdis);
  return vwdis;
}



/**********************************************************************************
 * Threaded versions
 */


/* a common function dispatches the calls */
typedef struct {

  /* input */

  int distance_type;  

  int nclust, d, k;
  const float *codebook;

  int npt;
  const float *points;
  const float *vw_weights;

  /* output */
  int *vw;
  float *vwdis;

  /* bookkeeping */
  int n_thread;
} nn_input_t;



static void nn_task (void *arg, int tid, int i)
{
  nn_input_t *t = arg;

  long n0 = t->npt * (long)i / t->n_thread;
  long n1 = t->npt * (long)(i + 1) / t->n_thread;

  knn_full (t->distance_type, n1 - n0, t->nclust, t->d, t->k, t->codebook,
			  t->points + n0 * t->d, t->vw_weights,
			  t->vw + n0 * t->k, t->vwdis + n0 * t->k);
}

/********** frontends */

void knn_full_thread (int distance_type, int npt, int nclust, int d, int k,
                                    const float *codebook, const float *coords,
                                    const float *vw_weights,
                                    int *vw, float *vwdis2,
                                    int n_thread) 
{
  if (npt < n_thread || n_thread == 1) {        /* too few pts */
    return knn_full (distance_type, npt, nclust, d, k, codebook, coords, vw_weights, 
                     vw, vwdis2);
  }

  nn_input_t task = { 
    distance_type,
    nclust, d, k, codebook, 
    npt, coords, vw_weights, vw, vwdis2,
    n_thread
  };

  compute_tasks (n_thread, n_thread, &nn_task, &task);

}


/***************** simplified calls */

float *knn_thread (int npt, int nclust, int d, int k,
		   const float *codebook, const float *coords,
		   int *vw, int n_thread) 
{
  float *vwdis2=fvec_new(k*npt);
  knn_full_thread (2, npt, nclust, d, k, codebook, coords, NULL, vw, vwdis2, n_thread);
  return vwdis2;
}



double nn_thread (int npt, int nclust, int d,
		  const float *codebook, const float *coords,
		  int *vw, int n_thread)
{
  float *vwdis2=fvec_new(npt);
  knn_full_thread (2, npt, nclust, d, 1, codebook, coords, NULL, vw, vwdis2, n_thread);
   
  double toterr = fvec_sum(vwdis2, npt); 

  free(vwdis2);
  return toterr;
}


/***************** cross distances */

typedef struct {
  int distance_type; /* <0 mreans optimized L2 distance */
  int d, na, nb;
  const float *a,*b;
  float *dist2;
  int nt;
  int split_a;  
} cross_distances_params_t;


/*Core funcs*/
void compute_cross_distances_task (void *arg, int tid, int i) 
{
  cross_distances_params_t *t=arg;
  if(t->split_a) {
    long i0 = t->na*(long)i/t->nt;
    long i1 = t->na*(long)(i+1)/t->nt;
    if(t->distance_type<0) 
      compute_cross_distances_nonpacked (t->d,i1-i0,t->nb,
                                         t->a+i0*t->d,t->d,
                                         t->b,t->d,
                                         t->dist2+i0,t->na);
    else 
      compute_cross_distances_alt_nonpacked (t->distance_type,t->d,i1-i0,t->nb,
                                             t->a+i0*t->d,t->d,
                                             t->b,t->d,
                                             t->dist2+i0,t->na);
  } else {
    long i0 = t->nb*(long)i/t->nt;
    long i1 = t->nb*(long)(i+1)/t->nt;
    if(t->distance_type<0) 
      compute_cross_distances_nonpacked (t->d,t->na,i1-i0,
                                         t->a,t->d,
                                         t->b+i0*t->d,t->d,
                                         t->dist2+i0*t->na,t->na);   
    else 
      compute_cross_distances_alt_nonpacked (t->distance_type,
                                             t->d,t->na,i1-i0,
                                             t->a,t->d,
                                             t->b+i0*t->d,t->d,
                                             t->dist2+i0*t->na,t->na);         
  }  
}

/*main funcs*/

void compute_cross_distances_thread (int d, int na, int nb,
                                     const float *a,
                                     const float *b, float *dist2,
                                     int nt) 
{
  cross_distances_params_t t={-1,d,na,nb,a,b,dist2,nt};
  
  int n=MAX(na,nb);
  
  if(n<nt) /* too small, no threads */
    compute_cross_distances(d,na,nb,a,b,dist2);
  else { 
    t.split_a=na>nb;    
    compute_tasks(nt,nt,&compute_cross_distances_task,&t);
  } 
}


void compute_cross_distances_alt_thread (int distance_type,int d, int na, int nb,
                                         const float *a,
                                         const float *b, float *dist2,
                                         int nt) 
{
  cross_distances_params_t t = {distance_type,d,na,nb,a,b,dist2,nt};
  
  int n=MAX(na,nb);
  
  if(n<nt) /* too small, no threads */
    compute_cross_distances_alt (distance_type,d,na,nb,a,b,dist2);
  else { 
    t.split_a=na>nb;    
    compute_tasks(nt,nt,&compute_cross_distances_task,&t);
  } 
}



#ifdef _OPENMP

#include <omp.h>

#define SET_NT  omp_set_num_threads(nt)  

#else 

#define SET_NT

/* #pragma's will be ignored */

#endif



void compute_distances_1_thread (int d, int nb,
                                 const float *a, 
                                 const float *b,                         
                                 float *dist2,
                                 int n_thread) {

  compute_distances_1_nonpacked_thread(d,nb,a,b,d,dist2,n_thread);  
}


void compute_distances_1_nonpacked_thread (int d, int nb,
                                           const float *a, 
                                           const float *b, int ldb, 
                                           float *dist2,
                                           int nt) {
  int i;
  
  SET_NT;
#pragma omp parallel 
  {
#pragma omp for 
    for(i=0;i<nt;i++) {
      int i0=i*nb/nt;
      int i1=(i+1)*nb/nt;
      compute_distances_1_nonpacked (d,i1-i0,a,
                                     b+i0*ldb,ldb,
                                     dist2+i0);
    }  

  }
}

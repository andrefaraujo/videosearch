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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "vector.h"
#include "matrix.h"
#include "sorting.h"
#include "machinedeps.h"
#include "eigs.h"

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))
#define NEWAC(type,n) (type*)calloc(sizeof(type),(n))
#define NEW(type) NEWA(type,1)


/* blas/lapack subroutines */

#define real float
#define integer FINTEGER

int sgemm_ (char *transa, char *transb, integer * m, integer *
            n, integer * k, real * alpha, const real * a, integer * lda,
            const real * b, integer * ldb, real * beta, real * c__,
            integer * ldc);

int ssyev_ (char *jobz, char *uplo, integer * n, real * a,
            integer * lda, real * w, real * work, integer * lwork,
            integer * info);


int sgeqrf_ (integer * m, integer * n, real * a, integer * lda,
             real * tau, real * work, integer * lwork, integer * info);

int slarft_ (char *direct, char *storev, integer * n, integer *
             k, real * v, integer * ldv, real * tau, real * t, integer * ldt);

int slarfb_ (char *side, char *trans, char *direct, char *storev, integer * m,
             integer * n, integer * k, real * v, integer * ldv, real * t,
             integer * ldt, real * c__, integer * ldc, real * work,
             integer * ldwork);

int ssyrk_(char *uplo, char *trans, integer *n, integer *k, 
           real *alpha, real *a, integer *lda, real *beta, real *c__, integer *
           ldc);


void sgemv_(const char *trans, integer *m, integer *n, real *alpha, 
                   const real *a, integer *lda, const real *x, integer *incx, real *beta, real *y, 
                   integer *incy);


int sgels_(char *trans, integer *m, integer *n, integer *
           nrhs, float *a, integer *lda, float *b, integer *ldb,
           float *work, integer *lwork, integer *info);


#undef real
#undef integer



/*---------------------------------------------------------------------------*/
/* Standard operations                                                       */
/*---------------------------------------------------------------------------*/

float *fmat_new (int nrow, int ncol)
{
  float *m = fvec_new (nrow * (long)ncol);
  return m;
}

float *fmat_new_0 (int nrow, int ncol)
{
  return fvec_new_0 (nrow * (long)ncol);
}


void fmat_mul_full(const float *left, const float *right,
                   int m, int n, int k,
                   const char *transp,
                   float *result) {

  fmat_mul_full_nonpacked(left, right, m, n, k, transp, 
                          (transp[0] == 'N' ? m : k), 
                          (transp[1] == 'N' ? k : n), 
                          result, m);
                       
  
}

void fmat_mul_full_nonpacked(const float *left, const float *right,
                             int mi, int ni, int ki,
                             const char *transp,
                             int ld_left, int ld_right, 
                             float *result,
                             int ld_result) {

  float alpha = 1;
  float beta = 0;
  FINTEGER m=mi,n=ni,k=ki;
  FINTEGER lda = ld_left;
  FINTEGER ldb = ld_right;
  FINTEGER ldc = ld_result; 
  
  sgemm_ ((char*)transp, (char*)(transp+1), &m, &n, &k,
          &alpha, left, &lda, right, &ldb, &beta, result, &ldc);

}


float* fmat_new_mul_full(const float *left, const float *right,
                         int m, int n, int k,
                         const char *transp) {
  float *result=fmat_new(m,n);

  fmat_mul_full(left, right, m, n, k, transp, result);

  return result;
}



void fmat_mul (const float *left, const float *right, int m, int n, int k, float *mout) {
  fmat_mul_full(left,right,m,n,k,"NN",mout);
}

void fmat_mul_tl (const float *left, const float *right, int m, int n, int k, float *mout) {
  fmat_mul_full(left,right,m,n,k,"TN",mout);
}

void fmat_mul_tr (const float *left, const float *right, int m, int n, int k, float *mout) {
  fmat_mul_full(left,right,m,n,k,"NT",mout);
}

void fmat_mul_tlr (const float *left, const float *right, int m, int n, int k, float *mout) {
  fmat_mul_full(left,right,m,n,k,"TT",mout);
}

float* fmat_new_mul (const float *left, const float *right, int m, int n, int k) {
  return fmat_new_mul_full(left,right,m,n,k,"NN");
}

float* fmat_new_mul_tl (const float *left, const float *right, int m, int n, int k) {
  return fmat_new_mul_full(left,right,m,n,k,"TN");
}

float* fmat_new_mul_tr (const float *left, const float *right, int m, int n, int k) {
  return fmat_new_mul_full(left,right,m,n,k,"NT");
}

float* fmat_new_mul_tlr (const float *left, const float *right, int m, int n, int k) {
  return fmat_new_mul_full(left,right,m,n,k,"TT");
}



void fmat_print (const float *a, int nrow, int ncol)
{
  int i, j;

  printf ("[");
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++)
      printf ("%.5g ", a[i + nrow * j]);
    if (i == nrow - 1)
      printf ("]\n");
    else
      printf (";\n");
  }
}

void fmat_print_tranposed(const float *a, int nrow, int ncol)
{
  long i, j;

  printf ("[");
  for (i = 0; i < nrow; i++) {
    for (j = 0; j < ncol; j++)
      printf ("%.5g ", a[i * ncol + j]);
    if (i == nrow - 1)
      printf ("]\n");
    else
      printf (";\n");
  }
}


/*---------------------------------------------------------------------------*/
/* Matrix manipulation functions                                             */
/*---------------------------------------------------------------------------*/


float *fmat_get_submatrix (const float *a, int nrow, 
                           int nrow_out,
                           int ncol) {
  long i;
  float *b=fmat_new(nrow_out,ncol);
  
  for(i=0;i<ncol;i++) 
    memcpy(b+i*nrow_out,a+i*nrow,nrow_out*sizeof(*a));

  return b;
}

int *imat_get_submatrix (const int *a, int nrow, 
                         int nrow_out,
                         int ncol) {
  long i;
  int *b=ivec_new(nrow_out*(long)ncol);
  
  for(i=0;i<ncol;i++) 
    memcpy(b+i*nrow_out,a+i*nrow,nrow_out*sizeof(*a));

  return b;
  
}


float *fmat_new_get_row (const float *a, int nrow, int ncol, int row)
{
  float *v = fvec_new (ncol);
  int j;
  for (j = 0 ; j < ncol ; j++) 
      v[j] = a[row + nrow * j];
  
  return v;
}

float *fmat_new_get_rows (const float *a, int d, int n,                              
			  int nrowout, const int *rows) {
  float *b=fmat_new(nrowout,n);
  long i, j;
  long ii=0;
  for(j=0;j<n;j++) {
    const float *aj = a + d*(long)j;
    for(i=0;i<nrowout;i++) 
      b[ii++]=aj[rows[i]];
  }
  
  return b;
}

void fmat_shuffle_columns(float *a, int nrow, int ncol) {
  long k,i;
  for (i = 0; i < ncol ; i++) {
    int j = i + random() % (ncol - i);
    /* swap i and j */
    float *ci=a+i*nrow;
    float *cj=a+(long)j*nrow;
    for(k=0;k<nrow;k++) {
      float tmp=ci[k];
      ci[k]=cj[k];
      cj[k]=tmp;
    }
  }
}


float *fmat_new_get_columns (const float *a, int nrow, int ncolout, const int *cols) {
  float *b = fmat_new (nrow, ncolout);
  fmat_get_columns(a, nrow, ncolout, cols, b);
  return b;
}

void fmat_get_rows_cols(const float *a, int d, 
                        int nrow, const int *rows, 
                        int ncol, const int *cols, 
                        float *out) {
  long i, j, k = 0;

  for(j = 0; j < ncol; j++) {
    const float *a_col = a + cols[j] * (long) d;
    for(i = 0; i < nrow; i++) {
      out[k++] = a_col[rows[i]];
    }
  }
  
}


void fmat_get_columns (const float *a, int d, int ncolout, const int *cols, float *b) {
  long j;
  for(j=0;j<ncolout;j++)
    memcpy(b + j * d, 
           a + cols[j] * (long) d, 
           d * sizeof(a[0]));
}


void fmat_sum_columns (const float * a, int nrow, int ncol, float * sums)
{
  long i, j;
  fvec_0 (sums, ncol);

  for(j=0 ; j<ncol ; j++)
    for (i=0 ; i<nrow ; i++) 
      sums[j] += a[nrow*j+i];
}


float *fmat_new_sum_columns (const float *a, int nrow, int ncol) 
{
  float *sums = fvec_new(ncol);

  fmat_sum_columns (a, nrow, ncol, sums);
  return sums;
}



void fmat_sum_rows (const float * a, int nrow, int ncol, float * sums)
{
  long i, j;
  fvec_0 (sums, nrow);

  for(j=0 ; j<ncol ; j++)
    for (i=0 ; i<nrow ; i++) 
      sums[i] += a[nrow*j+i];
}


float *fmat_new_sum_rows (const float *a, int nrow, int ncol) 
{
  float *sums = fvec_new(nrow);

  fmat_sum_rows (a, nrow, ncol, sums);
  return sums;
}

int fmat_remove_0_columns(float *a, int d, int n) {
  long nnz = 0, i; 
  for(i = 0; i < n; i++) {
    if(!fvec_all_0(a + d * i, d)) {
      if(nnz != i) 
        memcpy(a + d * nnz, a + d * i, sizeof(float) * d);
      nnz ++; 
    }    
  }
  return nnz;
}

void fmat_normalize_columns_l2sqr_pow(float *a, int d, int n, float pw) {
  long i; 
  if(pw == 0.5) {
    for(i = 0; i < n; i++) {
      double l2sqr = fvec_norm2sqr(a + d * i, d); 
      double norm = sqrt(l2sqr); 
      fvec_mul_by(a + d * i, d, norm); 
    }
  } else {
    for(i = 0; i < n; i++) {
      double l2sqr = fvec_norm2sqr(a + d * i, d); 
      double norm = pow(l2sqr, pw); 
      fvec_mul_by(a + d * i, d, norm); 
    }
  }
}

void fmat_normalize_columns_l2sqr_pow_robust(float *a, int d, int n, float pw, float eps) {
  long i; 
  if(pw == 0.5) {
    for(i = 0; i < n; i++) {
      double l2sqr = fvec_norm2sqr(a + d * i, d); 
      double norm = sqrt(l2sqr + eps); 
      fvec_mul_by(a + d * i, d, norm); 
    }
  } else {
    for(i = 0; i < n; i++) {
      double l2sqr = fvec_norm2sqr(a + d * i, d); 
      double norm = pow(l2sqr + eps, pw); 
      fvec_mul_by(a + d * i, d, norm); 
   }
  }
}



/*---------------------------------------------------------------------------*/
/* Special matrices                                                          */
/*---------------------------------------------------------------------------*/
float *fmat_new_rand_gauss (int nrow, int ncol)
{
  long i;
  float *m = fmat_new (nrow, ncol);

  for (i = 0; i < nrow * ncol; i++)
    m[i] = gaussrand ();

  return m;
}


/* method: we compute the QR decomposition of a matrix with Gaussian
   values */
float *random_orthogonal_basis (int di)
{ 
  FINTEGER d=di;
  long i;


  /* generate a Gaussian matrix */
  float *x = fmat_new_rand_gauss (d, d);

  float *tau = NEWA (float, d);

  {                             /* compute QR decomposition */

    /* query work size */
    float lwork_query;
    FINTEGER lwork = -1, info;
    sgeqrf_ (&d, &d, x, &d, tau, &lwork_query, &lwork, &info);
    assert (info == 0);

    lwork = (int) lwork_query;
    float *work = NEWA (float, lwork);
    sgeqrf_ (&d, &d, x, &d, tau, work, &lwork, &info);
    assert (info == 0);

    free (work);
  }

  /* Decomposition now stored in x and tau. Apply to identity to get
     explicit matrix Q */

  float *q = NEWAC (float, d * d);
  {

    float *t = NEWA (float, d * d);

    slarft_ ("F", "C", &d, &d, x, &d, tau, t, &d);

    for (i = 0; i < d; i++)
      q[i + d * i] = 1;

    float *work = NEWA (float, d * d);

    slarfb_ ("Left", "N", "F", "C",
             &d, &d, &d, x, &d, t, &d, q, &d, work, &d);

    free (t);
    free (work);
  }

  free (tau);
  free (x);
  return q;
}


/* Construct a Hadamard matrix of dimension d using the Sylvester construction.
   d should be a power of 2 */
float *hadamard (int d)
{
  assert ((d & (d - 1)) == 0 || !"d must be power of 2");

  int i, j;
  float *had = fvec_new (d * d);

  if (d == 1) {
    had[0] = 1;
    return had;
  }

  /* Compute the Hadamard matrix of dimension d / 2 */
  int dd = d / 2;
  float *had_part = hadamard (dd);

  for (i = 0; i < dd; i++)
    for (j = 0; j < dd; j++) {
      had[i * d + j] = had_part[i * dd + j];
      had[i * d + j + dd] = had_part[i * dd + j];
      had[(i + dd) * d + j] = had_part[i * dd + j];
      had[(i + dd) * d + j + dd] = -had_part[i * dd + j];
    }

  free (had_part);
  return (had);
}




/*---------------------------------------------------------------------------*/
/* Statistical matrix operations                                             */
/*---------------------------------------------------------------------------*/

float *fmat_center_columns(int d,int n,float *v) 
{
  assert(n>0);

  float *accu=fvec_new_cpy(v,d);
  long i;

  for(i=1;i<n;i++) 
    fvec_add(accu,v+i*d,d);

  fvec_div_by(accu,d,n);
  
  for(i=0;i<n;i++) 
    fvec_sub(v+i*d,accu,d);

  return accu;  
}

void fmat_subtract_from_columns(int d,int n,float *v,const float *avg) {
  long i;
  for(i=0;i<n;i++) 
    fvec_sub(v+i*d,avg,d);
}

void fmat_add_to_columns(int d,int n,float *v,const float *avg) {
  long i;
  for(i=0;i<n;i++) 
    fvec_add(v+i*d,avg,d);
}


void fmat_rev_subtract_from_columns(int d,int n,float *v,const float *avg) {
  long i;
  for(i=0;i<n;i++) 
    fvec_rev_sub(v+i*d,avg,d);

}


float *fmat_new_vstack(const float *a,int da,
                       const float *b,int db,
                       int n) {
  int i;
  float *c=fmat_new(da+db,n),*ret=c;
  for(i=0;i<n;i++) {
    memcpy(c,a,da*sizeof(float));
    c+=da;
    a+=da;
    memcpy(c,b,db*sizeof(float));
    c+=db;
    b+=db;
  }
  return ret;
}


void fmat_splat_separable(const float *a,int nrow,int ncol,
                          const int *row_assign,const int *col_assign,
                          int k,
                          float *accu) {
  long i,j;

  for(i=0;i<nrow;i++) for(j=0;j<ncol;j++) {
    accu[row_assign[i]*k+col_assign[j]]+=a[i*ncol+j];
  }

}

void fmat_splat_separable_1D(const float *a,int nrow,int ncol,
                             const int *assign,
                             float *accu) {
  long i, j;
  const float *acol = a;
  for(i = 0; i < ncol; i++) {
    if(assign[i] < 0) continue;
    float *ocol = accu + nrow * (long)assign[i];
    for(j = 0; j < nrow; j++) 
      ocol[j] += acol[j];
    acol += nrow;
  }  
}


int *imat_joint_histogram(int n,int k,int *row_assign,int *col_assign) {
  int *hist=ivec_new_0(k*k);
  int i;

  for(i=0;i<n;i++) 
    hist[row_assign[i]*(long)k+col_assign[i]]++;

  return hist;
}



/******************************************************************
 * Covariance and PCA computation
 *****************************************************************/

/* Input matrix: v(d,n) stored by rows.

   x is v data centered for each dimension 0<=j<n
   x = v - (1/n) * u * m 

   where :
   *   u(n,1) contains only 1's
   *   m = u' * v is the sum of values for each column of v 

   cov is the covariance matrix :

   cov = (1/n) x' * x
       = (1/n) v' * v - (1/n^2) * m' * m

   => no need to allocate an auxiliary array.
*/



float *fmat_new_covariance (int d, int n, const float *v, float *avg, int assume_centered)
{
  
  long i, j;

  float *cov = fvec_new_0 (d * d);
  
  if(!assume_centered) {

    float *sums = avg ? avg : fvec_new(d);
    fvec_0(sums,d);
    
    for (i = 0; i < n; i++)
      for (j = 0; j < d; j++)
        sums[j] += v[i * d + j];
    
    
    for (i = 0; i < d; i++)
      for (j = 0; j < d; j++)
        cov[i + j * d] = sums[i] * sums[j];
    
    
    if(avg)
      for(i=0;i<d;i++) avg[i]/=n;
    else
      free (sums);

  } 

  FINTEGER di=d,ni=n;

  if(0)  {
    float alpha = 1.0 / n, beta = -1.0 / (n * n);
    sgemm_ ("N", "T", &di, &di, &ni, &alpha, v, &di, v, &di, &beta, cov, &di);
  } else if(1) {
    /* transpose input matrix */
    float *vt=fvec_new(n*d);
    for(i=0;i<d;i++) 
      for(j=0;j<n;j++) 
        vt[i*n+j]=v[j*d+i];
    float alpha = 1.0 / n, beta = -1.0 / (n * n);
    
    sgemm_ ("T", "N", &di, &di, &ni, &alpha, vt, &ni, vt, &ni, &beta, cov, &di);
    
    free(vt);
  } else {
    float alpha = 1.0 / n, beta = -1.0 / (n * n);
    ssyrk_("L","N", &di, &ni, &alpha,(float*)v,&di,&beta,cov,&di);

    /* copy lower triangle to upper */

    for(i=0;i<d;i++)
      for(j=i+1;j<d;j++) 
        cov[i+j*d]=cov[j+i*d];

  }

  return cov;
}

float* fmat_new_transp (const float *a, int ncol, int nrow)
{
  int i,j;
  float *vt=fvec_new(ncol*nrow);

  for(i=0;i<ncol;i++) 
    for(j=0;j<nrow;j++) 
      vt[i*nrow+j]=a[j*ncol+i];

  return vt;
}



void fmat_pca_from_covariance(int d,const float *cov, float *singvals, float * pcamat) 
{
  float *evals=singvals;

  if(!singvals) evals=fvec_new(d);

  if(eigs_sym(d,cov,evals,pcamat)!=0) {
    free(pcamat);
    pcamat=NULL;
    goto error;
  }
  eigs_reorder(d,evals,pcamat,1); /* 1 = descending */

 error:
  if(!singvals) free(evals);
}



float *fmat_new_pca_from_covariance(int d,const float *cov, float *singvals) 
{
  float *pcamat=fvec_new(d*d);
  fmat_pca_from_covariance (d, cov, singvals, pcamat);
  return pcamat;
}




float *fmat_new_pca(int d,int n,const float *v, float *singvals) {

  float *cov=fmat_new_covariance(d,n,v,NULL,1);
  
  assert(fvec_all_finite(cov,d*d));
  
  float *evals=singvals;

  if(!singvals) evals=fvec_new(d);
  
  float *ret=fmat_new_pca_from_covariance(d,cov,evals);

  if(!singvals) free(evals);
    
  free(cov);  
  
  return ret;
}


















#ifdef _OPENMP

#include <omp.h>

#define SET_NT  omp_set_num_threads(nt)  

#else 

#define SET_NT

/* #pragma's will be ignored */

#endif

/* multithreaded matrix-vector multiply */

/* m=nb rows, n=nb cols */
void fmat_mul_v(int mi,int ni,const float*a,int ldai,
                 const float *x,
                 float *y,int nt) {
  int i;
  FINTEGER lda=ldai,n=ni,m=mi;
  
#pragma omp parallel num_threads(nt)
  {
#pragma omp for 
    for(i=0;i<nt;i++) {
      int i0=i*(long)m/nt;
      int i1=(i+1)*(long)m/nt;
      FINTEGER m1=i1-i0;
      float one=1.0,zero=0.0;
      FINTEGER ione=1;
      /* printf("%d %d\n",i,m1); */
      sgemv_("Trans",&n,&m1,&one,
             a+lda*(long)i0,&lda,x,&ione,&zero,y+i0,&ione);
      
    }
  }   

}

void fmat_mul_tv(int mi,int ni,const float*a,int ldai,
                 const float *x,
                 float *y,int nt) {
  int i,j;
  FINTEGER lda=ldai,n=ni,m=mi;
  
  float *ybuf=malloc(sizeof(float)*nt*m);
  assert(ybuf);

  if(nt>n) nt=n;
#pragma omp parallel num_threads(nt)
  {
#pragma omp for 
    for(i=0;i<nt;i++) {
      int i0=i*(long)n/nt;
      int i1=(i+1)*(long)n/nt;
      FINTEGER n1=i1-i0;
      float one=1.0,zero=0.0;
      FINTEGER ione=1;
      sgemv_("Not transposed",&m,&n1,&one,
             a+lda*(long)i0,&lda,x+i0,&ione,&zero,ybuf+i*(long)m,&ione);
      
    }  

  }  
  /* accumulate y results */
  memcpy(y,ybuf,sizeof(float)*m);
  float *yb=ybuf;
  for(i=1;i<nt;i++) {    
    yb+=m;
    for(j=0;j<m;j++) 
      y[j]+=yb[j];
  }

  free(ybuf);
}



int fmat_svd_partial_full(int n,int m,int nev,const float *a,int a_transposed,
                          float *s,float *vout,float *uout,int nt) {
  
  arpack_eigs_t *ae=arpack_eigs_begin(n,nev);
  if(!ae) return -100;
  int ret=0;
  
  int j,i;
  float *ax=NEWA(float,m);
  
  int it;

  for(it=0;;it++) {
    float *x,*y;
    ret=arpack_eigs_step(ae,&x,&y); 

    printf("arpack iteration %d ret=%d\r",it,ret);

    if(ret<0) break; /* error */

    if(ret==0) break; /* stop iteration */

    /* ret==1 */

    if(!a_transposed) {
      fmat_mul_v(m,n,a,n,x,ax,nt);
      fmat_mul_tv(n,m,a,n,ax,y,nt);
    } else {
      fmat_mul_tv(m,n,a,m,x,ax,nt);
      fmat_mul_v(n,m,a,m,ax,y,nt);
    }

    fflush(stdout);
  } 
  printf("\n");

  free(ax);

  float *v=vout ? vout : fmat_new(nev,n);

  ret=arpack_eigs_end(ae,s,v);

  if(ret>0) {
    int nconv=ret;
        
    if(s)
      for(j=0;j<nconv;j++) 
        s[j]=sqrt(s[j]);    

    if(uout) 
      for(i=0;i<nconv;i++) {
        float *u=uout+m*(long)i;
        if(!a_transposed)
          fmat_mul_v(m,n,a,n,v+n*(long)i,u,nt);
        else
          fmat_mul_tv(m,n,a,m,v+n*(long)i,u,nt);
        fvec_normalize(u,m,2);
      }               
    
  }

  if(!vout) free(v);
  
  return ret;
}

int fmat_svd_partial(int d,int n,int ns,const float *a,
                     float *singvals,float *u,float *v) {
  return fmat_svd_partial_full(d,n,ns,a,0,singvals,u,v,count_cpu());
}



float *fmat_new_pca_part(int d,int n,int nev,
                         const float *v,float *singvals) {

  if(!(nev<=d && nev<=n)) {
    fprintf(stderr,"fmat_new_pca_part: asking for too many eigenvalues (%d) wrt %d*%d data\n",nev,n,d);
    return NULL;
  }


  float *pcamat=fmat_new(d,nev);  

  
  int ret;

  if(n>=d) {
    ret=fmat_svd_partial_full(d,n,nev,v,0,singvals,pcamat,NULL,count_cpu());
  } else {
    fprintf(stderr,"fmat_new_pca_part: warn fewer learning points (%d) than dimensions (%d): transposing\n",n,d);
    
    ret=fmat_svd_partial_full(n,d,nev,v,1,singvals,NULL,pcamat,count_cpu());
  }

  if(ret<0) {
    free(pcamat); 
    pcamat=NULL;
  }

  return pcamat;
}



static int fmat_solve_ls_t_inplace(int mi, int ni, float *a, float *bx) {
  
  /* solve system */ 
  FINTEGER info;
  FINTEGER m=mi, nrhs=1, lda=mi, lwork=-1;
  FINTEGER n=ni, ldb=n;
  float work_sz;
  
  sgels_("Transposed", &m, &n, &nrhs, a, &lda, 
         bx, &ldb, &work_sz, &lwork, &info); 
 
  lwork = (long)work_sz;

  float *work = fvec_new(lwork);
 
  sgels_("Transposed", &m, &n, &nrhs, a, &lda, 
         bx, &n, work, &lwork, &info); 
  
  free(work);
  
  return info;
}

int fmat_solve_ls_t(int m, int n, const float *a, const float *b, float *x) {

  /* m < n : not enough space in x to store copy of b */

  float *aux = fvec_new(m * n + (m < n ? n : 0)); 

  float *a_copy = aux; 
  
  memcpy(a_copy, a, sizeof(float) * n * m); 
  
  float *bx = m < n ? aux + n * m : x; 

  memcpy(bx, b, sizeof(float) * n);
  
  int info = fmat_solve_ls_t_inplace(m, n, a_copy, bx); 

  if(m < n) 
    memcpy(x, bx, sizeof(float) * m);

  free(aux);

  return info;
}



/*--- Another way to do it by accumulating covariance matrice on-the-fly, using blocks of data ---*/

pca_online_t * pca_online_new (int d)
{
  pca_online_t * pca = (pca_online_t *) malloc (sizeof (pca_online_t));
  pca->d = d;
  pca->n = 0;
  pca->mu = fvec_new_0 (d);
  pca->cov = fvec_new_0 (d*(long)d);
  pca->eigvec = fvec_new (d*(long)d);
  pca->eigval = fvec_new (d);
  return pca;
}


void pca_online_delete (struct pca_online_s * pca)
{
  free (pca->mu);
  free (pca->cov);
  free (pca->eigvec);
  free (pca->eigval);
  free (pca);
}


/* Accumulate information for PCA for n input vectors */
void pca_online_accu (struct pca_online_s * pca, const float * v, long n)
{
  int d = pca->d;
  float * cov = fvec_new (d*(long)d);
  float * mu = fvec_new (d);

  fmat_sum_rows (v, d, n, mu);
  fmat_mul_tr (v, v, d, d, n, cov);

  fvec_add (pca->mu, mu, d);
  fvec_add (pca->cov, cov, d*(long)d);

  pca->n += n;

  free (cov);
  free (mu);
}


/* compute the mean and covariance matrix */
void pca_online_cov (struct pca_online_s * pca)
{
  int d = pca->d;
  int n = pca->n;

  fvec_div_by (pca->mu, d, n);
  fvec_div_by (pca->cov, d * (long)d, n);

  float * mumut = fvec_new (d*(long)d);
  fmat_mul_tr (pca->mu, pca->mu, d, d, 1, mumut);
  fvec_sub (pca->cov, mumut, d*(long)d);
  free (mumut);
  
  fvec_mul_by (pca->cov, d * (long)d, n / (double) (n-1));
  assert(fvec_all_finite(pca->cov,d*(long)d));
  pca->n = -pca->n;
}


/* compute the mean, the covariance matrix, and the eigenvectors.
   They are stored in the structure itself  */
void pca_online_complete (struct pca_online_s * pca)
{
  if (pca->n > 0)
    pca_online_cov (pca);
  fmat_pca_from_covariance (pca->d, pca->cov, pca->eigval, pca->eigvec);
}


/* compute the mean, the covariance matrix, and the eigenvectors.
   They are stored in the structure itself  */
void pca_online_complete_part (struct pca_online_s * pca, int nev)
{
  if (pca->n > 0)
    pca_online_cov (pca);

  if (nev * 2 >= pca->d)
    pca_online_complete (pca);
  else {
    int ret = eigs_sym_part(pca->d,pca->cov,nev,pca->eigval,pca->eigvec);
    assert (ret > 0);
    if(ret<nev) 
      printf("!!! only %d / %d eigenvalues converged\n",ret,nev);
    
  }
}


void pca_online_project (const pca_online_t * pca, const float * v, float * vo, int d, long n, int dout)
{
  const char trmat[2] = {'T', 'N'};
  float * vb = fvec_new_cpy (v, n*d);
  assert (d == pca->d);

  fmat_subtract_from_columns (pca->d, n, vb, pca->mu);
  fmat_mul_full (pca->eigvec, vb, dout, n, pca->d, trmat, vo);
  free (vb);
}

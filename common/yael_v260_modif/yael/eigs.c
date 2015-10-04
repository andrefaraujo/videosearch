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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "eigs.h"
#include "vector.h"
#include "sorting.h"
#include "machinedeps.h"



extern void dsyev_( char *jobz, char *uplo, FINTEGER *n, double *a, FINTEGER *lda,
        double *w, double *work, FINTEGER *lwork, FINTEGER *info );

extern void dsygv_(FINTEGER * itype, char *jobz, char *uplo, FINTEGER *n, double *a, FINTEGER *lda,
		    double *b, FINTEGER *lbd, double *w, double *work, FINTEGER *lwork, FINTEGER *info );


typedef float real;

extern void sgemv_(const char *trans, FINTEGER *m, FINTEGER *n, real *alpha, 
                   const real *a, FINTEGER *lda, const real *x, FINTEGER *incx, real *beta, real *y, 
                   FINTEGER *incy);




int eigs_sym (int di, const float * m, float * eigval, float * eigvec)
{
  int i, j;
  FINTEGER d=di;
  double * md = (double *) memalign (16, sizeof (*md) * d * d);

  /* processing is performed in double precision */
  for (i = 0 ; i < d ; i++) {
    for (j = 0 ; j < d ; j++)
      md[i * d + j] = (float) m[i * d + j];
  }

  /* variable for lapack function */
  double workopt = 0;
  FINTEGER lwork = -1, info;

  double * lambda = (double *) memalign (16, sizeof (*lambda) * d);
  dsyev_( "V", "L", &d, md, &d, lambda, &workopt, &lwork, &info );
  lwork = (int) workopt;
  double * work = (double *) memalign (16, lwork * sizeof (*work));
  dsyev_( "V", "L", &d, md, &d, lambda, work, &lwork, &info );
  
  if (info > 0) {
    fprintf (stderr, "# eigs_sym: problem while computing eigen-vectors/values info=%d\n",info);
    goto error;
  }
  /* normalize the eigenvectors, copy and free */
  double nr = 1;
  for (i = 0 ; i < d ; i++) {
    if(eigval)
      eigval[i] = (float) lambda[i];
    
    if(eigvec) 
      for (j = 0 ; j < d ; j++) 
        eigvec[i * d + j] = (float) (md[i * d + j] / nr);
  }
 error:
  free (md);
  free (lambda);
  free (work);
  return info;
}


int geigs_sym (int di, const float * a, const float * b, float * eigval, float * eigvec)
{
  int i, j;
  FINTEGER d=di;
  double * ad = (double *) memalign (16, sizeof (*ad) * d * d);
  double * bd = (double *) memalign (16, sizeof (*bd) * d * d);

  /* processing is performed in double precision */
  for (i = 0 ; i < d ; i++) 
    for (j = 0 ; j < d ; j++) {
      ad[i * d + j] = (float) a[i * d + j];
      bd[i * d + j] = (float) b[i * d + j];
    }
  
  /* variable for lapack function */
  double workopt = 0;
  FINTEGER lwork = -1, info, itype = 1;

  double * lambda = (double *) memalign (16, sizeof (*lambda) * d);
  dsygv_ (&itype, "V", "L", &d, ad, &d, bd, &d, lambda, &workopt, &lwork, &info );
  lwork = (int) workopt;
  double * work = (double *) memalign (16, lwork * sizeof (*work));
  dsygv_ (&itype, "V", "L", &d, ad, &d, bd, &d, lambda, work, &lwork, &info );
  
  if (info != 0) {
    fprintf (stderr, "# eigs_sym: problem while computing eigen-vectors/values info=%d\n",info);
    goto error;
  }

  /* normalize the eigenvectors, copy and free */
  double nr = 1;
  for (i = 0 ; i < d ; i++) {

    if(eigval)
      eigval[i] = (float) lambda[i];
    
    if(eigvec) 
      for (j = 0 ; j < d ; j++) 
        eigvec[i * d + j] = (float) (ad[i * d + j] / nr);
  }

 error:
  free (ad);
  free (bd);
  free (lambda);
  free (work);
  return info;
}



void eigs_reorder (int d, float * eigval, float * eigvec, int criterion)
{
  int i;
  int * perm = ivec_new (d);

  float * eigvalst = fvec_new (d);
  float * eigvecst = fvec_new (d * d);

  fvec_sort_index (eigval, d, perm);

  if (criterion) 
    for (i = 0 ; i < d / 2 ; i++) {
      int tmp = perm[i];
      perm[i] = perm[d - 1 - i];
      perm[d - 1 - i] = tmp;
    }

  for (i = 0 ; i < d ; i++) {
    eigvalst[i] = eigval[perm[i]];
    memcpy (eigvecst + i * d, eigvec + perm[i] * d, sizeof (*eigvecst) * d);
  }

  memcpy (eigval, eigvalst, d * sizeof (*eigval));
  memcpy (eigvec, eigvecst, d * d * sizeof (*eigvec));

  free (eigvalst);
  free (eigvecst);
  free (perm);
}




int eigs_sym_part (int ni, const float * a, int nev, float * sout, float * vout) {
  FINTEGER n=ni;
  arpack_eigs_t *ae=arpack_eigs_begin(n,nev);
  int ret=0;
  

  for(;;) {
    float *x,*y;
    ret=arpack_eigs_step(ae,&x,&y); 

    if(ret<0) break; /* error */

    if(ret==0) break; /* stop iteration */

    /* ret==1 */

    float zero=0,one=1;
    FINTEGER ione=1;
    
    sgemv_("Trans",&n,&n,&one,a,&n,x,&ione,&zero,y,&ione);    
  } 
  ret=arpack_eigs_end(ae,sout,vout);
 
  return ret;
}




#ifdef HAVE_ARPACK



typedef FINTEGER integer;
typedef FINTEGER logical;

extern void ssaupd_ (integer *ido,const char*bmat,integer *n, const char*which,integer *nev,
                     float* tol, float*resid, integer *ncv, float *v, integer *ldv, 
                     integer *iparam, integer * ipntr, float *workd, float *workl, 
                     integer *lworkl, integer *info );


extern void sseupd_ (logical *rvec, const char *howmny, logical *select, float *d    ,
                     float *z     ,integer *ldz   , float *sigma , const char*bmat,
                     integer *n       , const char*which,integer *nev, float* tol, 
                     float*resid, integer *ncv, float *v, integer *ldv, 
                     integer *iparam, integer * ipntr, float *workd, float *workl, 
                     integer *lworkl, integer *info );

struct arpack_eigs_t {
  FINTEGER n,nev;

  FINTEGER ncv;
  FINTEGER ido,info;

  FINTEGER lworkl;
  float *resid,*workd,*workl;
  float *v;
  FINTEGER *iparam,*ipntr;
  logical *select;

};

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))

arpack_eigs_t *arpack_eigs_begin(int n,int nev) {
  arpack_eigs_t *ae=NEWA(arpack_eigs_t,1);

  ae->n=n;
  ae->nev=nev;

  int ncv=ae->ncv=2*nev;  /* should be enough (see remark 4 of ssaupd doc) */

  ae->lworkl = ncv*(long)(ncv+8);
  ae->resid=NEWA(float,n);
  ae->workd=NEWA(float,3*n);
  ae->workl=NEWA(float,ae->lworkl);
  
  ae->v=NEWA(float,n*(long)ncv);
  FINTEGER *iparam=ae->iparam=NEWA(FINTEGER,11);
  ae->ipntr=NEWA(FINTEGER,11);

  ae->info=0; /* use random initial vector */
  ae->ido=0;

  iparam[0]=1;
  iparam[2]=n;
  iparam[6]=1;

  return ae;
}


int arpack_eigs_step(arpack_eigs_t *ae,
                     float **x, float **y) {
  
  const char *bmat="I",*which="LM";
  
  float tol=0;
  
  ssaupd_(&ae->ido, bmat, &ae->n, which, &ae->nev, 
          &tol, ae->resid, &ae->ncv, ae->v, &ae->n, 
          ae->iparam, ae->ipntr, ae->workd, ae->workl, &ae->lworkl,
          &ae->info);
  

  if(ae->ido==-1 || ae->ido==1) {
    *x=ae->workd+ae->ipntr[0]-1;
    *y=ae->workd+ae->ipntr[1]-1;
    return 1;
  } 
  
  *x=*y=NULL;
  if(ae->info<0) {
    printf("arpack_eigs_step: ssaupd_ error info=%d\n",ae->info);
    
    return ae->info;
  } 

  return 0; 
}


int arpack_eigs_end(arpack_eigs_t *ae,
                     float * sout, float * vout) {
  int i,ret=0;  
  FINTEGER n=ae->n,nev=ae->nev,ncv=ae->ncv;
  int nconv;

  logical *select=NEWA(logical,ncv);
  float *s=NEWA(float,ncv*2);
  int *perm=NEWA(int,nev);

  if(ae->info<0) {
    ret=ae->info;
    goto error;
  }

  {
    FINTEGER ierr;
    logical rvec=1;
    float sigma;
    const char *bmat="I",*which="LM";
    float tol=0;
  
    sseupd_(&rvec,"All",select, s,
            ae->v,&n, &sigma, bmat, &n, which, &nev, 
            &tol, ae->resid, &ncv, ae->v, &n, 
            ae->iparam, ae->ipntr, ae->workd, ae->workl, &ae->lworkl,
            &ierr);

    if(ierr!=0) {
      printf("arpack_eigs_end: sseupd_ error: %d\n",ierr);
      ret=ierr;
      goto error;
    }
    ret=nconv=ae->iparam[4];
    assert(nconv<=nev);

  }

  /* order v by s */

  fvec_sort_index(s,nconv,perm); 
  
  if(vout) 
    for(i=0;i<nconv;i++) 
      memcpy(vout+n*(long)i, ae->v+n*(long)(nconv-1-perm[i]), sizeof(float)*n);

  if(sout) 
    for(i=0;i<nconv;i++) 
      sout[i]=s[nconv-1-perm[i]];

  

 error: 
  free(select); 
  free(perm); 
  free(s);

  free(ae->resid); 
  free(ae->workl);
  free(ae->workd);
  free(ae->iparam);
  free(ae->ipntr);
  free(ae->v);

  free(ae);

  return ret;
}



#else

arpack_eigs_t *arpack_eigs_begin(int n,int nev) {
  fprintf(stderr,"Error: Yael not compiled with Arpack!");
  abort();
} 

int arpack_eigs_step(arpack_eigs_t *ae,
                     float **x, float **y) {
  fprintf(stderr,"Error: Yael not compiled with Arpack!");
  abort();
}

int arpack_eigs_end(arpack_eigs_t *ae,
                    float * sout, float * vout) {
  fprintf(stderr,"Error: Yael not compiled with Arpack!");
  abort();
}



#endif

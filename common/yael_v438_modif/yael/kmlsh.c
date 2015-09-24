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
#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "vector.h"
#include "nn.h"
#include "sorting.h"
#include "kmeans.h"

#include "kmlsh.h"


/*---------- Utils ----------*/

static inline int get_maxincell (const int * boundaries, int nclust)
{
  int c, nincell, maxnidx = 0;
  for (c = 0 ; c < nclust ; c++) {
    nincell = boundaries[c+1] - boundaries[c];
    if (nincell > maxnidx)
      maxnidx = nincell;
  }
  return maxnidx;
}



/*---------- k-NN list handling ----------*/

nnlist_t * nnlist_new_noalloc (int n, int k)
{
  nnlist_t * l = (nnlist_t *) malloc (sizeof (nnlist_t));
  l->n = n;
  l->k = k;
  return l;
}


nnlist_t * nnlist_new (int n, int k)
{
  nnlist_t * l = (nnlist_t *) malloc (sizeof (nnlist_t));
  l->n = n;
  l->k = k;
  l->idx = ivec_new_set (n * (long) k, -1);
  l->dis = fvec_new_set (n * (long) k, 1e30);
  return l;
}


void nnlist_delete (nnlist_t * l)
{
  free (l->dis);
  free (l->idx);
  free (l);
}


/* A few function to store and update the list of NN */
static inline void nnlist_add (nnlist_t * l, int lno, int elidx, float eldis)
{
  int * idx = l->idx + (long) lno * l->k;
  float * dis = l->dis + (long) lno * l->k;

  /* Check if distance is lowest than greatest value */
  if (eldis > dis[0])
    return;

  /* check if idx can be found somewhere. If so, skip it */
  if (ivec_index (idx, l->k, elidx) != -1)
    return;

  idx[0] = elidx;
  dis[0] = eldis;

  /* The first element is the keeper, and should be the max value. */
  long posmax = fvec_arg_max (dis, l->k);
  if (posmax == 0)
    return;

  idx[0] = idx[posmax];
  dis[0] = dis[posmax];
  idx[posmax] = elidx;
  dis[posmax] = eldis;
}


void nnlist_addn (nnlist_t * l, int lno, int n, int * idx, float * dis)
{
  int i; 
  for (i = 0 ; i < n ; i++)
    nnlist_add (l, lno, idx[i], dis[i]);
}


/*---------- k-means LSH structure handling ----------*/


kmlsh_t * kmlsh_new (int nhash, int nclust, int d)
{
  int i;
  kmlsh_t * lsh = (kmlsh_t *) malloc (sizeof (kmlsh_t));

  lsh->d = d;
  lsh->nhash = nhash;
  lsh->nclust = nclust;
  lsh->centroids = (float **) malloc (nhash * sizeof (float *));
  
  /* Make a single malloc (so lazy) */
  for (i = 0 ; i < nhash ; i++)
    lsh->centroids[i] = fvec_new (d * nclust);

  return lsh;
}


void kmlsh_delete (kmlsh_t * lsh)
{
  int i;
  for (i = 0 ; i < lsh->nhash ; i++)
    free (lsh->centroids[i]);
  free (lsh->centroids);
  free (lsh);
}



/* Learn several k-means using different sampling strategies on the learning vectors */
/* n is the number of vectors possibly used as input of k-means, while 
   nlearn is the number of vectors actually used for the k-means.       */
/* WARNING: this function is called "bvec" for the moment to avoid name conflict 
   with Yael's function. The name will change when the function will be integrated in yael */
void kmlsh_learn_xvec (kmlsh_t * lsh, int n, int nlearn, const void * v, 
		       int flags, int vec_type)
{
  int h;

  /* k-means parameters */
  int nb_iter_max = KMLSH_NB_ITER_MAX;
  int d = lsh->d;
  int nclust = lsh->nclust;
  int nt = flags & KMLSH_NT;
  int kmeans_flags = nt | KMEANS_INIT_RANDOM;
  int verbose = !(flags & KMLSH_QUIET);

  /* Weird way to do polymorphism in C*/
  const float * vf = (float *) v;
  const unsigned char * vb = (unsigned char *) v;

  if (nlearn == 0) {
    nlearn = lsh->nclust * 50; 
    if (nlearn > n) nlearn = (n * 3) / 4;
  }
  assert (nlearn <= n);
  if (verbose)
    fprintf (stderr, "n=%d  nlearn=%d  nclust=%d\n", n, nlearn, lsh->nclust);
  else kmeans_flags |= KMEANS_QUIET;

  float * vlearn = fvec_new (nlearn * d);
  for (h = 0 ; h < lsh->nhash ; h++) {
    /* Construct a first subset of vectors for learning with target size */
    int * perm = ivec_new_random_idx (n, nlearn);

    if (vec_type == KMLSH_VECTYPE_FVEC)
      fvec_cpy_subvectors (vf, perm, d, nlearn, vlearn);
    else if (vec_type == KMLSH_VECTYPE_BVEC)
      b2fvec_cpy_subvectors (vb, perm, d, nlearn, vlearn);
    else assert (0);
    free (perm);

    /* perform the k-means based on the selected vectors */
    if (verbose)
      fprintf (stderr, "Learn K-means using %d vectors\n", nlearn);
    kmeans (d, nlearn, nclust, nb_iter_max, vlearn, kmeans_flags, 0, 1, 
	    lsh->centroids[h], NULL, NULL, NULL);

    /* Optionnally, write the intermediate quantization indexes */
    if (flags & KMLSH_WRITE_INTER_NHASH) {
      char * stmp = (char *) malloc (256);
      snprintf (stmp, 256, "tmp/nhash-intermediate_%d.kmlsh", h+1);
      kmlsh_write (stmp, lsh);
      free (stmp);
    }
  }
  free (vlearn);
}
 

/* Same as kmlsh_learn, but also create the structure */
kmlsh_t * kmlsh_new_learn_bvec (int nhash, int nclust, int d, int n, int nlearn, 
				const unsigned char * v, int flags)
{
  kmlsh_t * lsh = kmlsh_new (nhash, nclust, d);
  kmlsh_learn_xvec (lsh, n, nlearn, v, flags, KMLSH_VECTYPE_BVEC);
  return lsh;
}


/* Same as kmlsh_learn, but also create the structure */
kmlsh_t * kmlsh_new_learn_fvec (int nhash, int nclust, int d, int n, int nlearn, 
				const float * v, int flags)
{
  kmlsh_t * lsh = kmlsh_new (nhash, nclust, d);
  kmlsh_learn_xvec (lsh, n, nlearn, v, flags, KMLSH_VECTYPE_FVEC);
  return lsh;
}


/* Quantize the descriptors and order them by cell. */
void kmeans_cohash_xvec (const kmlsh_t * lsh, int h, const void * v, int n, 
			 int * perm, int * boundaries, int flags, int vec_type)
{
  long i, j;
  int * idx = ivec_new (n);    /* To store index id */
  float * dis = fvec_new (n);  /* to store (unused) distances to k-NN */
  int nt = flags & KMLSH_NT;

  /* Kmeans config */
  int d = lsh->d;
  int nclust = lsh->nclust;

  /* Weird way to do polymorphism in C*/
  float * vf = (float *) v;
  const unsigned char * vb = (unsigned char *) v;
  float * vbuf = NULL;

  if (vec_type != KMLSH_VECTYPE_FVEC)
    vbuf = fvec_new (KMLSH_BLOCK_SIZE * d);   /* buffer to store the vectors cast to float */


  /* assign all the vectors using this space partitioning */
  for (i = 0 ; i < n ; i += KMLSH_BLOCK_SIZE) {
    long ninblock = KMLSH_BLOCK_SIZE;
    if (i + ninblock > n)
      ninblock = n - i;

    if (vec_type == KMLSH_VECTYPE_BVEC)
      for (j = 0 ; j < ninblock * d; j++)
        vbuf[j] = vb[i*d+j];
    else if (vec_type == KMLSH_VECTYPE_FVEC)
      vbuf = vf + i * d;
    knn_full_thread (2, ninblock, nclust, d, 1, lsh->centroids[h], 
		     vbuf, NULL, idx + i, dis + i, nt);
    if (! (flags & KMLSH_QUIET))
      fprintf (stdout, "\rQuantize %d descriptors. %6.2f%%", n, 100.0 * (i + ninblock) / (float) n);
  } 
  if (! (flags & KMLSH_QUIET))
    printf ("\n");
  if (vec_type != KMLSH_VECTYPE_FVEC)
    free (vbuf);

  int * histoidx = ivec_new_histogram (nclust + 1, idx, n);
  int a = histoidx[0], b, c; 
  histoidx[0] = 0;

  for (c = 1 ; c <= nclust ; c++) {
    b = histoidx[c];
    histoidx[c] = histoidx[c - 1] + a;
    a = b;
  }
  
  int * nocc = ivec_new_0 (nclust);
  for (i = 0 ; i < n ; i++) {
    int cell = idx[i];
    int pos = histoidx[cell] + nocc[cell]++;
    assert (nocc[cell] <= histoidx[cell+1] - histoidx[cell]);
    perm[pos] = i;
  }

  for (i = 0 ; i < nclust ; i++)
    assert (nocc[i] == histoidx[i+1] - histoidx[i]);

  ivec_cpy (boundaries, histoidx, nclust + 1);
  
  free (idx);
  free (dis);
  free (histoidx);
  free (nocc);
}


void kmeans_cohash_bvec (const kmlsh_t * lsh, int h, const unsigned char * v, int n, 
			 int * perm, int * boundaries, int flags)
{
  kmeans_cohash_xvec (lsh, h, (void *) v, n, perm, boundaries, flags, KMLSH_VECTYPE_BVEC);
}


void kmeans_cohash_fvec (const kmlsh_t * lsh, int h, const float * v, int n, 
			 int * perm, int * boundaries, int flags)
{
  kmeans_cohash_xvec (lsh, h, (void *) v, n, perm, boundaries, flags, KMLSH_VECTYPE_FVEC);
}


kmlsh_idx_t * kmlsh_idx_new (const kmlsh_t * lsh, int n)
{
  kmlsh_idx_t * lshidx = (kmlsh_idx_t *) malloc (sizeof (struct kmlsh_idx_s));

  lshidx->nhash = lsh->nhash;
  lshidx->n = n;
  lshidx->nclust = lsh->nclust;
  lshidx->perm = ivec_new (lsh->nhash * n);
  lshidx->boundaries = ivec_new (lsh->nhash * (lsh->nclust+1));
  return lshidx;
}


void kmlsh_idx_delete (kmlsh_idx_t * lshidx)
{
  free (lshidx->perm);
  free (lshidx->boundaries);
  free (lshidx);
}


kmlsh_idx_t * kmlsh_idx_new_compile_xvec (const kmlsh_t * lsh, const void * v, int n, int flags, int vec_type)
{
  int h;
  kmlsh_idx_t * lshidx = kmlsh_idx_new (lsh, n);
  int nt = flags & KMLSH_NT;

  for (h = 0 ; h < lsh->nhash ; h++) {
    int * perm = lshidx->perm + h * n;
    int * boundaries = lshidx->boundaries + h * (lsh->nclust + 1);
    
    kmeans_cohash_xvec (lsh, h, v, n, perm, boundaries, nt, vec_type);
    
    /* Optionnally, write the intermediate quantization indexes */
    if (flags & KMLSH_WRITE_INTER_NHASH) {
      char * stmp = (char *) malloc (256);
      snprintf (stmp, 256, "tmp/nhash-intermediate_%d.kmlshidx", h+1);
      kmlsh_idx_write (stmp, lshidx);
      free (stmp);
    }
  }

  return lshidx;
}

kmlsh_idx_t * kmlsh_idx_new_compile_bvec (const kmlsh_t * lsh, const unsigned char * v, int n, int flags)
{
  return kmlsh_idx_new_compile_xvec (lsh, (void *) v, n, flags, KMLSH_VECTYPE_BVEC);
}

kmlsh_idx_t * kmlsh_idx_new_compile_fvec (const kmlsh_t * lsh, const float * v, int n, int flags)
{
  return kmlsh_idx_new_compile_xvec (lsh, (void *) v, n, flags, KMLSH_VECTYPE_FVEC);
}

/* Return the number of vectors assigned to cell c for hash function h */
int kmlsh_idx_get_nvec (const kmlsh_idx_t * lshidx, int h, int c)
{
  int nvecincell = lshidx->boundaries[h * (lshidx->nclust+1) + c + 1] 
    - lshidx->boundaries[h * (lshidx->nclust+1) + c];
  return nvecincell;
}

/* Maximum number of vectors in the cell */
int kmlsh_idx_get_maxincell (const kmlsh_idx_t * lshidx, int h)
{
  int maxincell = get_maxincell (lshidx->boundaries + h * (lshidx->nclust+1), lshidx->nclust);
  return maxincell;
}



/* Return a pointer to the idxs of vectors in cell c for hash function h.
   Do not allocate any memory */
int * kmlsh_idx_get_vecids (const kmlsh_idx_t * lshidx, int h, int c)
{
  return lshidx->perm + lshidx->n * h + lshidx->boundaries[h * (lshidx->nclust + 1) + c];
}


nnlist_t * kmlsh_match_xvec (const kmlsh_t * lsh,
			const kmlsh_idx_t * lshidx_b, const void * vb, int nb,
			const kmlsh_idx_t * lshidx_q, const void * vq, int nq,
			int k, int nt, int vec_type)
{
  int h;
  long i, c;

  int d = lsh->d;
  int nclust = lsh->nclust;

  /* Weird way to do polymorphism in C*/
  const float * vfb = (float *) vb;
  const float * vfq = (float *) vq;
  const unsigned char * vbb = (unsigned char *) vb;
  const unsigned char * vbq = (unsigned char *) vq;

  /* Structure to store the list of NN hypothesis */
  nnlist_t * nnlist = nnlist_new (nq, k);

  for (h = 0 ; h < lsh->nhash ; h++) {
    /* Partition the base/query vectors based on cell co-locality */
    long maxnidx_b = kmlsh_idx_get_maxincell (lshidx_b, h);
    long maxnidx_q = kmlsh_idx_get_maxincell (lshidx_q, h);

    /* for each cluster, group the vectors and compute the exact knn-graph
       within the quantization cell. This is used to update the full graph */
    float * vbuf_b = fvec_new (maxnidx_b * d);
    float * vbuf_q = fvec_new (maxnidx_q * d);

    printf ("max in histo: base=%ld  query=%ld  (nclust=%d, nb=%d, nq=%d)\n", 
	    maxnidx_b, maxnidx_q, nclust, nb, nq);

    /* temporary arrays to store the NN variables for each cell */
    int * idxtmp = ivec_new (maxnidx_q * k);
    float * distmp = fvec_new (maxnidx_q * k);

    printf ("Group the vectors based on cell %d / %d, then exact search\n", h + 1, lsh->nhash);

    for (c = 0 ; c < nclust ; c++) {
      long nvecincell_b = kmlsh_idx_get_nvec (lshidx_b, h, c);
      long nvecincell_q = kmlsh_idx_get_nvec (lshidx_q, h, c);
      if (nvecincell_q == 0 || nvecincell_b == 0)
        continue;
       
      assert (nvecincell_b <= maxnidx_b);
      assert (nvecincell_q <= maxnidx_q);

      int * vidx_b = kmlsh_idx_get_vecids (lshidx_b, h, c); 
      int * vidx_q = kmlsh_idx_get_vecids (lshidx_q, h, c); 

      if (vec_type == KMLSH_VECTYPE_FVEC) {
        fvec_cpy_subvectors (vfb, vidx_b, d, nvecincell_b, vbuf_b);
        fvec_cpy_subvectors (vfq, vidx_q, d, nvecincell_q, vbuf_q);
      } else if (vec_type == KMLSH_VECTYPE_BVEC) {
        b2fvec_cpy_subvectors (vbb, vidx_b, d, nvecincell_b, vbuf_b);
        b2fvec_cpy_subvectors (vbq, vidx_q, d, nvecincell_q, vbuf_q);
      } else assert (0);

      /* Call the exact kNN-graph function that is applied for each cell */
      int k2 = (k < nvecincell_b ? k : nvecincell_b);
      knn_full_thread (2, nvecincell_q, nvecincell_b, d, k2, vbuf_b, vbuf_q, 
		       NULL, idxtmp, distmp, nt);  


      /* translate output of knn-graph to absolute vector indexes */
      for (i = 0 ; i < k2 * nvecincell_q ; i++) 
        idxtmp[i] = vidx_b[idxtmp[i]];

      /* update the list of NN from previous hash functions */
      for (i = 0 ; i < nvecincell_q ; i++) 
        nnlist_addn (nnlist, vidx_q[i], k2, idxtmp + i * k2, distmp + i * k2);
    }

    free (vbuf_b);
    free (vbuf_q);
    free (idxtmp);
    free (distmp);
  }
  return nnlist;
}


nnlist_t * kmlsh_match_bvec (const kmlsh_t * lsh,
            const kmlsh_idx_t * lshidx_b, const unsigned char * vb, int nb,
            const kmlsh_idx_t * lshidx_q, const unsigned char * vq, int nq,
            int k, int nt)
{
  return kmlsh_match_xvec (lsh, lshidx_b, (void *) vb, nb, lshidx_q, 
			   (void *) vq, nq, k, nt, KMLSH_VECTYPE_BVEC);
}


nnlist_t * kmlsh_match_fvec (const kmlsh_t * lsh,
            const kmlsh_idx_t * lshidx_b, const float * vb, int nb,
            const kmlsh_idx_t * lshidx_q, const float * vq, int nq,
            int k, int nt)
{
  return kmlsh_match_xvec (lsh, lshidx_b, (void *) vb, nb, lshidx_q, 
			   (void *) vq, nq, k, nt, KMLSH_VECTYPE_FVEC);
}


nnlist_t * kmlsh_ann_xvec (const void * vb, int nb,
		      const void * vq, int nq,
		      int d, int k, int nhash, int nt, int vec_type)
{
  /* pre-defined parameters */
  int nclust = (int) sqrt (nb);
  int nlearn = nclust * 100;
  if (nlearn > nb) {
    fprintf (stderr, "Warning: should use more vector in learning set\n");
    nlearn = nb;
  }

  kmlsh_t * lsh = kmlsh_new (nhash, nclust, d);
  kmlsh_learn_xvec (lsh, nb, nlearn, vb, nt, vec_type);

  /* compute the hash values for database vectors and queries */
  kmlsh_idx_t * lshidx_b = kmlsh_idx_new_compile_xvec (lsh, vb, nb, nt, vec_type);
  kmlsh_idx_t * lshidx_q = kmlsh_idx_new_compile_xvec (lsh, vq, nq, nt, vec_type);

  
  /* Perform the matching */
  nnlist_t * nnlist = kmlsh_match_xvec (lsh, lshidx_b, vb, nb, lshidx_q, vq, nq, k, nt, vec_type);

  kmlsh_delete (lsh);
  
  return nnlist;
}


nnlist_t * kmlsh_ann_bvec (const unsigned char * vb, int nb,
              const unsigned char * vq, int nq,
              int d, int k, int nhash, int nt)
{
  return kmlsh_ann_xvec ( (void *) vb, nb, (void *) vq, nq, d, k, nhash, nt, KMLSH_VECTYPE_BVEC);
}


nnlist_t * kmlsh_ann_fvec (const float * vb, int nb,
              const float * vq, int nq,
              int d, int k, int nhash, int nt)
{
  return kmlsh_ann_xvec ( (void *) vb, nb, (void *) vq, nq, d, k, nhash, nt, KMLSH_VECTYPE_FVEC);
}

/*--------------------------------------------------------------*/
/* Various Input/Output functions                               */

#define KMLSH_WRITE_ERROR(test) {if (!test) {  \
  fprintf (stderr, "# kmlsh_write: I/O error with file %s\n", filename); \
  exit (1); }}								

#define KMLSH_READ_ERROR(test) {if (!test) {  \
  fprintf (stderr, "# kmlsh_read: I/O error with file %s\n", filename); \
  exit (1); }}
  
#define KMLSH_IDX_WRITE_ERROR(test) {if (!test) {  \
  fprintf (stderr, "# kmlsh_idx_write: I/O error with file %s\n", filename); \
  exit (1); }}								

#define KMLSH_IDX_READ_ERROR(test) {if (!test) {  \
  fprintf (stderr, "# kmlsh_idx_read: I/O error with file %s\n", filename); \
  exit (1); }}
  



void kmlsh_write (const char * filename, const kmlsh_t * lsh)
{
  int h, ret;
  FILE * f = fopen (filename, "w");
  KMLSH_WRITE_ERROR (f);

  for (h = 0 ; h < lsh->nhash ; h++) {
    ret = fvec_fwrite(f, lsh->centroids[h], lsh->nclust * lsh->d);
    KMLSH_WRITE_ERROR (ret == 0);
  }
  fclose (f);
}


void kmlsh_read (const char * filename, const kmlsh_t * lsh)
{
  int h, ret;
  FILE * f = fopen (filename, "r");
  KMLSH_READ_ERROR (f);

  for (h = 0 ; h < lsh->nhash ; h++) {
    lsh->centroids[h] = fvec_new (lsh->d * lsh->nclust);
    ret = fvec_fread (f, lsh->centroids[h], lsh->d * lsh->nclust);
    KMLSH_READ_ERROR (ret == lsh->nclust * lsh->d);
  }
  fclose (f);
}


void kmlsh_idx_write (const char * filename, const kmlsh_idx_t * lshidx)
{
  long h, ret;
  FILE * f = fopen (filename, "w");
  KMLSH_IDX_WRITE_ERROR (f);


  for (h = 0 ; h < lshidx->nhash ; h++) {
    ret = ivec_fwrite (f, lshidx->perm + h * lshidx->n, lshidx->n);
    KMLSH_IDX_WRITE_ERROR (ret == 0);

    ret = ivec_fwrite (f, lshidx->boundaries + h * (lshidx->nclust + 1), lshidx->nclust + 1);
    KMLSH_IDX_WRITE_ERROR (ret == lshidx->nclust + 1);
  }
  fclose (f);
}


void kmlsh_idx_read (const char * filename, kmlsh_idx_t * lshidx)
{
  long h, ret;
  FILE * f = fopen (filename, "r");
  KMLSH_IDX_READ_ERROR (f);

  for (h = 0 ; h < lshidx->nhash ; h++) {
    ret = ivec_fread (f, lshidx->perm + h * lshidx->n, lshidx->n);
    KMLSH_IDX_READ_ERROR (ret == lshidx->n);

    ret = ivec_fread (f, lshidx->boundaries + h * (lshidx->nclust + 1), lshidx->nclust + 1);
    KMLSH_IDX_READ_ERROR (ret == lshidx->nclust + 1);
  }

  fclose (f);
}

#undef KMLSH_WRITE_ERROR
#undef KMLSH_READ_ERROR  
#undef KMLSH_IDX_WRITE_ERROR
#undef KMLSH_IDX_READ_ERROR

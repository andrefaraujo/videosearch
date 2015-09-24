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

#ifndef __kmlsh_h
#define __kmlsh_h

/*---------------------------------------------------------------------------*/
/*! @addtogroup kmlsh
 *  @{
 */


/*! @defgroup kmlsh
 * K-means LSH is an implementation of the technique described in the following paper
 * "Locality sensitive hashing: a comparison of hash function types and querying mechanisms",
 * by L. Pauleve, H. Jegou and L. Amsaleg, Pattern Recognition Letters, August 2010
 * 
 * Only the regular LSH (multiple hash functions, no multi-probe, no query adaptive)
 * is provided. This implementation is not intended for a classical query/database 
 * scenario, although it could be used for (with relatively low efficienty). 
 * Instead, it is optimized towards batch processing of large amounts of queries. 
 */


/*-----------------------------------------------------------------*/
/*! A structure to handle the list of KNN                          */
struct nnlist_s {
  long n;        /* number of points */
  long k;        /* number of nearest neighbors */
  int * idx;     /* indices of the NN */
  float * dis;   /* corresponding distances */
};

typedef struct nnlist_s nnlist_t;

/* allocate n lists of length k */
nnlist_t * nnlist_new (int n, int k);

/* same but don't allocate the index/dis lists */
nnlist_t * nnlist_new_noalloc (int n, int k);

/* free the k-NN list structure */
void nnlist_delete (nnlist_t * l);

/* add n elements to the list */
void nnlist_addn (nnlist_t * l, int lno, int n, int * idx, float * dis);


/*-----------------------------------------------------------------*/
/* A km-LSH structure defining multiple k-means                    */

#define KMLSH_NT                    0x000000ff
#define KMLSH_QUIET                 0x00010000
#define KMLSH_WRITE_INTER_NHASH     0x00020000   /* write all intermediate versions 
						    of all functions and idx structures */

#define KMLSH_BLOCK_SIZE     256
#define KMLSH_NB_ITER_MAX    8

#define KMLSH_VECTYPE_FVEC   0
#define KMLSH_VECTYPE_BVEC   1


/*! The structure that contains the parameters of the KM-LSH */
struct kmlsh_s {
  int nhash;           /* number of hash functions */
  int d;               /* vector dimensionality */
  int nclust;          /* number of cluster per cell */
  float ** centroids;  /* all centroids */
};

typedef struct kmlsh_s kmlsh_t;


/*! A structure containing the pre-processed data 
  (tables of quantized indexes) for a set of vectors  */
struct kmlsh_idx_s {
  int nhash;
  int n;               /* number of vectors stored */
  int nclust; 
  int * perm;          /* the vector ids, ordered by quantization index (all hash tables) */
  int * boundaries;  
};

typedef struct kmlsh_idx_s kmlsh_idx_t;


/* alloc the kmlsh_t structure */
kmlsh_t * kmlsh_new (int nhash, int nclust, int d);

/* free the kmlsh_t structure */
void kmlsh_delete (kmlsh_t * lsh);


/* Learn several k-means using different sampling strategies on the learning vectors */
/* n is the number of vectors possibly used as input of k-means, while 
   nlearn is the number of vectors actually used for the k-means (typically n/2). */
void kmlsh_learn_xvec (kmlsh_t * lsh, int n, int nlearn, const void * v, 
		       int flags, int vec_type);

/* Same as kmlsh_learn, but also create the structure */
kmlsh_t * kmlsh_new_learn_bvec (int nhash, int nclust, int d, int n, int nlearn, 
				const unsigned char * v, int flags);

kmlsh_t * kmlsh_new_learn_fvec (int nhash, int nclust, int d, int n, int nlearn, 
				const float * v, int flags);

/*! A function that performs the match assuming that the codes are pre-computed */
nnlist_t * kmlsh_match_xvec (const kmlsh_t * lsh,
			const kmlsh_idx_t * lshidx_b, const void * vb, int nb,
			const kmlsh_idx_t * lshidx_q, const void * vq, int nq,
			int k, int nt, int vec_type);

nnlist_t * kmlsh_match_bvec (const kmlsh_t * lsh,
            const kmlsh_idx_t * lshidx_b, const unsigned char * vb, int nb,
            const kmlsh_idx_t * lshidx_q, const unsigned char * vq, int nq,
            int k, int nt);

nnlist_t * kmlsh_match_fvec (const kmlsh_t * lsh,
            const kmlsh_idx_t * lshidx_b, const float * vb, int nb,
            const kmlsh_idx_t * lshidx_q, const float * vq, int nq,
            int k, int nt);

/* Approximate search with pre-defined parameters. 
   The parameter nhash controls the trade-off quality/efficiency/memory (number of hash functions).
   flags is mainly use to set the number of processor cores */
nnlist_t * kmlsh_ann_xvec (const void * vb, int nb,
		      const void * vq, int nq,
		      int d, int k, int nhash, int flags, int vec_type);

nnlist_t * kmlsh_ann_bvec (const unsigned char * vb, int nb,
              const unsigned char * vq, int nq,
              int d, int k, int nhash, int flags);

nnlist_t * kmlsh_ann_fvec (const float * vb, int nb,
              const float * vq, int nq,
              int d, int k, int nhash, int flags);


/* Alloc/Free the index associated with a KM-LSH structure */
kmlsh_idx_t * kmlsh_idx_new (const kmlsh_t * lsh, int n);
void kmlsh_idx_delete (kmlsh_idx_t * lshidx);


/* compute the metadata associated with index of a KM-LSH structure */
kmlsh_idx_t * kmlsh_idx_new_compile_xvec (const kmlsh_t * lsh, const void * v, int n, int flags, int vec_type);

kmlsh_idx_t * kmlsh_idx_new_compile_bvec (const kmlsh_t * lsh, const unsigned char * v, int n, int flags);

kmlsh_idx_t * kmlsh_idx_new_compile_fvec (const kmlsh_t * lsh, const float * v, int n, int flags);

/* Return the number of vectors assigned to cell c for hash function h */
int kmlsh_idx_get_nvec (const kmlsh_idx_t * lshidx, int h, int c);

/* Maximum number of vectors in the cell */
int kmlsh_idx_get_maxincell (const kmlsh_idx_t * lshidx, int h);

/* Return a pointer to the idxs of vectors in cell c for hash function h.
   Do not allocate any memory (do not free the vector, free the structure instead) */
int * kmlsh_idx_get_vecids (const kmlsh_idx_t * lshidx, int h, int c);


/* Quantize the descriptors and order them by cell. */
void kmeans_cohash_xvec (const kmlsh_t * lsh, int h, const void * v, int n, 
			 int * perm, int * boundaries, int flags, int vec_type);

void kmeans_cohash_bvec (const kmlsh_t * lsh, int h, const unsigned char * v, int n, 
			 int * perm, int * boundaries, int flags);

void kmeans_cohash_fvec (const kmlsh_t * lsh, int h, const float * v, int n, 
			 int * perm, int * boundaries, int flags);


/*-----------------------------------------------------------------*/
/* I/O                                                             */

/* write the kmlsh_t structure on disk.  */
void kmlsh_write (const char * filename, const kmlsh_t * lsh);

/* read the kmlsh_t structure from disk (must be allocated with kmlsh_new) */
void kmlsh_read (const char * filename, const kmlsh_t * lsh);

/* write the kmlsh_idx_t index codes on disk. */
void kmlsh_idx_write (const char * filename, const kmlsh_idx_t * lshidx);

/* write the kmlsh_idx_t index codes on disk */
void kmlsh_idx_read (const char * filename, kmlsh_idx_t * lshidx);



/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif

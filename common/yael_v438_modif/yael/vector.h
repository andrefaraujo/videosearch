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

#ifndef __vector_h
#define __vector_h

#include <stdio.h>

/*-------------- Basic math operations ----------------*/


/*! generate a random sample, mean 0 variance 1 */
double gaussrand ();


/*---------------------------------------------------------------------------*/
/*! @addtogroup vector
 *  @{
 */



/*! @defgroup vector 
 * Vectors are represented as C arrays of basic elements. Functions
 * operating on them are prefixed with:
 *
 * ivec_: basic type is int
 *
 * fvec_: basic type is float
 *
 * Vector sizes are passed explicitly, as long int's to allow for
 * large arrays on 64 bit machines. Vectors can be free'd with free().
 *
 *
 * Arrays of vectors are stored contiguously in memory. An array of n
 * float vectors of dimension d is 
 * 
 *   float *fv
 *
 * The i'th element of vector j of vector array vf, where 0 <= i < d
 * and 0 <= j < n is
 * 
 *   vf[ j * d + i ]
 * 
 * It can also be seen as a column-major matrix of size d, n.
 *
 */ 


/*! Alloc a new aligned vector of floating point values -- to be
 *  de-allocated with free. Some operations may be faster if input
 *  arrays are allocated with this function (data is suitably
 *  aligned). */
float * fvec_new (long n);

/*! Alloc an int array -- to be de-allocated with free. */
int *ivec_new (long n);

/*! Alloc an byte array -- to be de-allocated with free. */
unsigned char *bvec_new (long n);

/*! Alloc a long array -- to be de-allocated with free */
long long * lvec_new (long n);

/*! Alloc an int array -- to be de-allocated with free. */
double * dvec_new (long n);

/*! create a vector initialized with 0's. */
float * fvec_new_0 (long n);
double * dvec_new_0 (long n);

/*! create a vector initialized with 0's. */
int *ivec_new_0 (long n);
unsigned char *bvec_new_0 (long n);

/*! create a vector initialized with 0's. */
long long *lvec_new_0 (long n);

/*! create a vector initialized with NaN (to trace errors) */
float *fvec_new_nan (long n);

/*!  create a vector initialized with a specified value. */
float *fvec_new_set (long n, float val);

/*!  create a vector initialized with a specified value. */
int *ivec_new_set (long n, int val);

/*!  create a vector initialized with uniformly drawn samples in [0,1) */
float *fvec_new_rand (long n);

/*!  create a vector initialized with gaussian samples */
float *fvec_new_randn (long n);


/*!  same as fvec_randn, with seed for thread-safety */
void fvec_randn_r (float * v, long n, unsigned int seed); 

/*!  same as fvec_new_rand, with seed for thread-safety */
float *fvec_new_rand_r (long n, unsigned int seed);

/*!  same as fvec_new_randn, with seed for thread-safety */
float *fvec_new_randn_r (long n, unsigned int seed);

/*!  new vector [a,a+1,...b-1] */
int * ivec_new_range (long a, long b);

/*!  new vector initialized with another vector */
int * ivec_new_cpy (const int * v, long n);

/*!  new vector initialized with another vector */
float * fvec_new_cpy (const float * v, long n);

/*!  random permutation of 0..n-1 */ 
int *ivec_new_random_perm (int n);

/*!  select k random indexes among n (without repetition) */ 
int *ivec_new_random_idx  (int n, int k);

/*!  same as ivec_new_random_perm, thread-safe, with a random seed */ 
int *ivec_new_random_perm_r (int n, unsigned int seed);

/*!  same as ivec_new_random_idx, thread-safe with a random seed  */ 
int *ivec_new_random_idx_r  (int n, int k, unsigned int seed);

/*! resize a vector (realloc). Usage: v = fvec_resize (v, n). */
float * fvec_resize (float * v, long n);

/*! resize a vector (realloc). Usage: v = fvec_resize (v, n). */
int * ivec_resize (int * v, long n);

/*!  count occurrences
   @param k is the range of the values that may be encountered (assuming start at 0). Values outside the range trigger an assertion!
   @param v is the vector of values to be histrogramized, of length n
*/
int * ivec_new_histogram (int k, const int * v, long n);

/*!  same as ivec_new_histogram, but values falling out of range are clipped (counted in the nearest bin) */
int * ivec_new_histogram_clip (int k, int * v, long n);


/*! count occurrences: maps [vmin,vmax) to 0..k-1
   @param vmin   min val of range
   @param vmax   max val of range
   @param k      nb of bins 
   @param v is the vector of values to be histrogramized, of length n
*/
int * fvec_new_histogram_clip (float vmin,float vmax, int k, float *v, long n);


/*!  compute a hash value for the vector */
int ivec_hash (const int * v, long n);

/*!  all occurences of a value by another in a vector */
void ivec_replace (int * v, long n, int val, int replace_val);

/*!  count occurences of a value in the vector */
long ivec_count_occurrences (const int * v, long n, int val);
long fvec_count_occurrences (const float * v, long n, float val);

/*!  count the number of values below a threshold */
long fvec_count_lt(const float * v, long n, float val);
long ivec_count_lt(const int * v, long n, int val);

/*!  count number of values above a threshold */
long fvec_count_gt(const float * v, long n, float val);
long ivec_count_gt(const int * v, long n, int val);

/*!  count number of values in a range (min <= x < max) */
long fvec_count_inrange(const float * v, long n, float vmin, float vmax);
long ivec_count_inrange(const int * v, long n, int vmin, int vmax);

/*! count the number of nan values */
long fvec_count_nan (const float * v, long n);
long fvec_count_nonfinite (const float * v, long n);

long fvec_count_0 (const float *val, long n);

/*---------------------------------------------------------------------------*/
/* Input/Output functions                                                    */
/* I/O of a single vector is supported only if it is smaller than 2^31       */ 
/*---------------------------------------------------------------------------*/


/*!  Read the number of vectors in a file and their dimension 
  (vectors of same size). Output the number of bytes of the file. */
long fvecs_fsize (const char * fname, int *d_out, int *n_out);
long ivecs_fsize (const char * fname, int *d_out, int *n_out);
long bvecs_fsize (const char * fname, int *d_out, int *n_out);
long lvecs_fsize (const char * fname, int *d_out, int *n_out);


/*!  write a vector into an open file */
int ivec_fwrite(FILE *f, const int *v, int d);
int fvec_fwrite(FILE *f, const float *v, int d);

/*!  write a vector without the dimension header */
int ivec_fwrite_raw (FILE *f, const int *v, long d);
int fvec_fwrite_raw (FILE *f, const float *v, long d);
int bvec_fwrite_raw (FILE *f, const unsigned char *v, long d);

int ivec_write_raw (const char * fname, const int *v, long d);
int fvec_write_raw (const char * fname, const float *v, long d);
int bvec_write_raw (const char * fname, const unsigned char *v, long d);

/*!  write a set of vectors into an open file */
int ivecs_fwrite(FILE *f, int d, int n, const int *v);
int fvecs_fwrite (FILE *fo, int d, int n, const float *vf);

/*!  several integer vectors of identical length into an file */
int ivecs_write (const char *fname, int d, int n, const int *v);
int ivecs_write_txt (const char * fname, int d, int n, const int *v);
int fvecs_write (const char *fname, int d, int n, const float *vf);
int fvecs_write_txt (const char * fname, int d, int n, const float *vf);

int bvecs_write (const char *fname, int d, int n, const unsigned char *v);

/*!  load float vectors from file.
 *
 * Returns nb of vectors read, or <0 on error */
int fvecs_new_read (const char *fname, int *d_out, float **vf);

int fvecs_new_fread_max (FILE *f, int *d_out, float **vf, long nmax);


/*! mmap vectors from a file. 
 * 
 * The returned memory area is read-only. 
 *
 * WARNING, the i'th element of vector j of vector array vf, where
 * 0 <= i < d and 0 <= j < n is
 * 
 *   vf[ j * (d + 1) + i ]
 * 
 * (mind the d+1)
 * the file remains open and there is no deallocation function (yet)
 */
int fvecs_new_mmap (const char *fname, int *d_out, float **vf);

int ivecs_new_mmap (const char *fname, int *d_out, int **vi);




/* The behavior of bvecs_new_read in not completely consistent 
   with the one of fvecs_new_read (can not read stream)         */
int bvecs_new_read (const char *fname, int *d_out, unsigned char **v_out);

int lvecs_new_read (const char *fname, int *d_out, long long **v_out);

/*! load a file of byte vectors, and convert them to float on-the-fly */
int b2fvecs_new_read (const char *fname, int *d_out, float **v_out);

/*! reads sparse vectors and return them as dense. d must be known */
int fvecs_new_read_sparse (const char *fname, int d, float **vf_out);

/*! read siftgeo, return as bvecs + metadata as fvecs */
int bvecs_new_from_siftgeo(const char *fname, 
			   int *d_v_out, unsigned char **v_out,
			   int *d_meta_out, float **meta_out); 

/*!  load float vector without allocating memory 
 *
 * Fills n*d array with as much vectors read from fname as possible.
 * Returns nb of vectors read, or <0 on error.
 */
int fvecs_read (const char *fname, int d, int n, float *v);

/*! read some vector from a bvec file and put them into a fvec vector */
int b2fvecs_read (const char *fname, int d, int n, float *v);


/*!  read a vector from a text file (one line per vector) */
int fvecs_read_txt (const char *fname, int d, int n, float *v);

/*!  read a single vector from a file

 * Fill a with a single float vector from fname offset o_f into file a
 * Returns <0 on error
 */
int fvec_read (const char *fname, int d, float *a, int o_f);

/*!  load float vectors from an open file. Return the dimension */
int fvec_fread (FILE * f, float * v, int d_alloc);

/* Read raw vectors from disk */

int fvec_fread_raw (FILE *f, float * v, long n);
int ivec_fread_raw (FILE *f, int * v, long n);
int bvec_fread_raw (FILE *f, unsigned char * v, long n);


float *fvec_new_fread_raw(FILE * f, long n);
int * ivec_new_fread_raw(FILE * f, long d);
unsigned char *bvec_new_fread_raw(FILE * f, long n);

float * fvec_new_read_raw(const char * fname, long d);
int * ivec_new_read_raw(const char * fname, long d);
unsigned char *bvec_new_read_raw(const char * fname, long d);

/*! load a set of n vectors from an open file. 
  Return the number of vectors that have been read. */
long fvecs_fread (FILE * f, float * v, long n, int d_alloc);

long ivecs_fread (FILE * f, int * v, long n, int d_alloc);

long bvecs_fread (FILE * f, unsigned char * v, long n, int d_alloc);

long lvecs_fread (FILE * f, long long * v, long n, int d_alloc);

long b2fvecs_fread (FILE * f, float * v, long n);


/*!  read and allocate a an integer vector file */
int * ivec_new_read(const char *fname, int *d_out);

/*!  read an integer vector file from an open file and return the dimension */
int ivec_fread (FILE *f, int * v, int d_alloc);

/*!  read a byte vector file from an open file and return the dimension */
int bvec_fread (FILE *f, unsigned char * v, int d_alloc);
int b2fvec_fread (FILE * f, float * v);

/*!  read an long vector file from an open file and return the dimension */
int lvec_fread (FILE *f, long long * v, int d_alloc);


/*!  read several integer vectors from an ivec file. Return number read */
int ivecs_new_read (const char *fname, int *d_out, int **vi);

/*! load a few of ivecs */
int ivecs_new_fread_max (FILE *f, int *d_out, int **vi, long nmax);


/*!  display a float vector */
void fvec_print (const float * v, int n);
void fvec_fprintf (FILE * f, const float *v, int n, const char *fmt);

/*!  display an integer vector */
void ivec_print (const int * v, int n);
void ivec_fprintf (FILE * f, const int *v, int n, const char *fmt);

/*!  find first index of val (return -1 if not found) */
long ivec_index (const int * v, long n,int val);


/*---------------------------------------------------------------------------*/
/* Cast vector type                                                          */
/*---------------------------------------------------------------------------*/

/*! cast a vector of int into a (new) vector of floats */
float * ivec2fvec (const int * v, long n);

/*! cast a vector of int into a (new) vector of floats */
float * bvec2fvec (const unsigned char * v, long n);

/*! cast a vector of int into a vector of floats. No internal allocation. */
void bvectofvec (const unsigned char * v, float * vb, long n);

void fvectodvec (const float *a, double *b, long n); 


/*---------------------------------------------------------------------------*/
/* Vector manipulation and elementary operations                             */
/*---------------------------------------------------------------------------*/

/*!  Set all the components of the vector v to 0 */
void fvec_0 (float * v, long n);
void ivec_0 (int * v, long n);

void fvec_nan (float * v, long n);

/*! Set all the components of the vector to random values in [0,1[ */
void fvec_rand (float * v, long n);

/*! Set all the components of the vector to gaussian values */
void fvec_randn (float * v, long n);


/*!  are all values 0? */
int fvec_all_0 (const float * v, long n);
int ivec_all_0 (const int * v, long n);

/*! are all vals >= 0? */
int fvec_all_ge0 (const float * v, long n);
int ivec_all_ge0 (const int * v, long n);

/*!  are all values finite? */
int fvec_all_finite (const float * v, long n);

/*!  Set all the components of the vector v to the value val */
void fvec_set (float * v, long n, float val);
void ivec_set (int * v, long n, int val);

/*!  copy the vector from v2 to v1 */
void ivec_cpy (int * vdest, const int * vsource, long n);
void fvec_cpy (float * vdest, const float * vsource, long n);
void bvec_cpy (unsigned char * vdest, const unsigned char * vsource, long n);

/*!  Increment or decrement a vector by a scalar value */
void fvec_incr (float * v, long n, double scal);
void fvec_decr (float * v, long n, double scal);
void ivec_incr (int * v, long n, int scal);
void ivec_decr (int * v, long n, int scal);


/*!  Multiply or divide a vector by a scalar */
void fvec_mul_by (float * v, long n, double scal);
void fvec_div_by (float * v, long n, double scal);
/* resciprocal division */
void fvec_rdiv_by (float * v, long n, double scal);

/*!  Add or subtract two vectors. The result is stored in v1. */
void fvec_add (float * v1, const float * v2, long n);
void fvec_sub (float * v1, const float * v2, long n);

/*! v1 := v2-v1 */
void fvec_rev_sub (float * v1, const float * v2, long n);

/*! v1 := v1 + v2 * scal */
void fvec_add_mul (float * v1, const float * v2, long n, double scal);

/*!  Component-wise multiplication or division of two vectors (result in v1) */
void fvec_mul (float * v1, const float * v2, long n);
void fvec_div (float * v1, const float * v2, long n);

/*!  Normalize the vector for the given Minkowski norm. 
  The function return the norm of the original vector. 
  If the vector is all 0, it will be filled with NaNs. 
  This case can be identified when the return value is 0. 
  Infinty norm can be obtained with norm=-1 */
double fvec_normalize (float * v, long n, double norm);

/*!  This function normalize a set of n d-dimensional vectors. 
  It returns the number of vectors whose norms was 0 (for which 
  the normalization has put some NaN values). */
int fvecs_normalize (float * v, long n, long d, double norm);

void fvec_round (float * v, long n);
void fvec_sqrt (float * v, long n);
void fvec_sqr (float * v, long n);
void fvec_exp (float * v, long n);
void fvec_log (float * v, long n);
void fvec_neg (float * v, long n);

/*!  signed square-root: y = sign(x)*sqrt(abs(x)) */
void fvec_ssqrt (float * v, long n);

/*! signed power: y = sign(x) * pow(abs(x), scal) */
void fvec_spow (float * v, long n, double scal);

/*! 2-stage normalization (like Lowe's SIFT normalization) */
void fvec_normalize_2stage(float * v, long n, double scal);

void ivec_add (int * v1, const int * v2, long n);
void ivec_sub (int * v1, const int * v2, long n);
void ivec_mul_by (int * v1,long n, int scal);
void ivec_mod_by (int * v1,long n, int scal);
void ivec_add_scalar (int * v, long n, int scal);
void fvec_add_scalar (float * v, long n, float scal);


/*! Replace the "Not a number" values by a given value */
long fvec_purge_nans(float * v, long n, float replace_value);
int fvec_purge_nonfinite(float * v, long n, float replace_value);

/*!  Shrink the vector, removing "Not a number" and inf values. 
  Returns new size */
long fvec_shrink_nonfinite(float * v, long n);

/*!  find 1st occurrence of a non-finite element */
long fvec_index_nonfinite (float * v, long n);

/*! revert order of vector */
void fvec_revert(float *v, long n); 

/*! swap two vectors */
void fvec_swap(float *v1, float *v2, long n); 


/*---------------------------------------------------------------------------*/
/* Vector measures and statistics                                            */
/*---------------------------------------------------------------------------*/

/*!  compute the sum of the value of the vector */
double fvec_sum (const float * v, long n);
long long ivec_sum (const int * v, long n);

/*! cumulative sum */
void fvec_cumsum(float * v, long n);
void ivec_cumsum(int *v, long n); 

/*! opposite of cumsum: v[i] := v[i] - v[i-1] */
void fvec_cumdiff(float *v, long n);
void ivec_cumdiff(int *v, long n);

/*!  compute the sum of the product of the vector */
double fvec_product (const float * v, long n);
long long ivec_product (const int * v, long n);

/*!  sum of squared components */
double fvec_sum_sqr (const float * v, long n);
long long ivec_sum_sqr (const int * v, long n);

/*!  compute the sum of the value of the vector */
double fvec_mean (const float * v, long n);
long long ivec_mean (const int * v, long n);

/*!  compute the norm of a given vector (norm=-1 => infinty norm) */
double fvec_norm (const float * v, long n, double norm);

/*!  compute squared norm 2 */
double fvec_norm2sqr (const float * v, long n);

/*!  count the number of non-zeros elements */
long fvec_nz (const float * v, long n);
long ivec_nz (const int * v, long n);

/*!  compute the positions of the non-null positions.
  return the number of non-zeros positions. */
int fvec_find (const float *v, int n, int ** nzpos_out);
int ivec_find (const int *v, int n, int ** nzpos_out);

/*!  perform a random permutation on the elements of the vector */
void ivec_shuffle (int *v, long n);

/*!  entropy of the probability mass function represented by the vector */
double fvec_entropy (const float *pmf, int n);

/*!  entropy of a binary variable */
double binary_entropy (double p);

double ivec_unbalanced_factor(const int *hist, long n);

/*---------------------------------------------------------------------------*/
/* Distances and similarities                                                */
/*---------------------------------------------------------------------------*/

/*!  Return the Hamming distance (i.e., the number of different elements) */
long ivec_distance_hamming (const int * v1, const int * v2, long n);

/*!  Return the L2 distance between vectors */
double fvec_distance_L2 (const float * v1, const float * v2, long n);

/*!  Return the L1 distance between vectors */
double fvec_distance_L1 (const float * v1, const float * v2, long n);

/*!  Return the square L2 distance between vectors */
double fvec_distance_L2sqr (const float * v1, const float * v2, long n);

/*!  inner product between two vectors */
double fvec_inner_product (const float * v1, const float * v2, long n);

/*---------------------------------------------------------------------------
 * Sparse vector handling
 *
 * sparse vectors are represented with: int *idx, float *val, int nz.
 * 
 * Element idx[i] of vector is val[i], for i=0..nz-1
 *
 * for i=1..nz-1,  idx[i-1]<idx[i]
 *---------------------------------------------------------------------------*/


/*!  convert a vector to a sparse vector. 
  Return the number of non-zeros positions  */
int fvec_to_spfvec (float * v, int n, int ** idx_out, float ** v_out);
int ivec_to_spivec (int * v, int n, int ** idx_out, int ** v_out);


/*!  convert a sparse vector into a full vector */
float * spfvec_to_fvec (int * idx, float * v, int nz, int n);
int * spivec_to_ivec (int * idx, int * v, int nz, int n);

/*!  inner product between two sparse vectors */
float spfvec_inner_product (int *idx1, float *val1, int nz1, 
			    int *idx2, float *val2, int nz2);


/*---------------------------------------------------------------------------*/
/* Elaborate vector manipulations                                            */
/*---------------------------------------------------------------------------*/

/*! on output,
 *
 * sl_out[0] =               v[      0] + ... + v[sl[0]-1]
 *
 * sl_out[i] = sl_out[i-1] + v[sl[i-1]] + ... + v[sl[i]-1] for 0<i<n
 */
void ivec_accumulate_slices(const int *v,int *sl,int n); 


/*!  mapping operator: dest[i]:=src[map[i]] for i=0..n-1 */
void fvec_map (const float *src,const int *map,int n,float *dest);

/*!  mapping operator: dest[i]:=src[map[i]] for i=0..n-1 */
void ivec_map (const int *src,const int *map,int n,int *dest);

/*!  inverse mapping operator: dest[imap[i]]:=src[i] for i=0..n-1 */
void fvec_imap(const float *src,const int *imap,int n,float *dest);

/*! 
 * for i=0..n-1, do 
 *    accu[assign[i]] += a[i]
 */ 
void fvec_splat_add(const float *a,int n,
                    const int *assign,float *accu); 

/*! 
 * for i=0..n-1, do 
 *    accu[i] += a[assign[i]]
 */ 
void fvec_isplat_add(const float *a,int n,
                     const int *assign,float *accu); 

/*! return input vector duplicated n times, with a value added each time 
 * 
 * @param nrepeat   nb of times to repeat input vector
 * @param inc       inc*i is added to all elements of i^th repeated vector
 *
 * @return          vector of size n*nrepeat
 */
int* ivec_repeat_with_inc(const int *a,int n,
                          int nrepeat, int inc);

/*! Copy the set of nout vectors in v (seen as a set of vectors), indexed by idx, in vout
 * @param v         the set of input vectors
 * @param idx       the indexes of the vectors to be copied
 * @param d         vectors' dimensionality
 * @param nout      number of vectors copied
 * @param vout      output vector (must be allocated externally with size nout)
 */  
void fvec_cpy_subvectors (const float * v, int * idx, int d, int nout, float * vout);

void b2fvec_cpy_subvectors (const unsigned char * v, int * idx, int d, int nout, float * vout);

/*! simple type conversion */
void ivec_to_fvec(const int *v, float *f, long n); 


/*! @} */
#endif

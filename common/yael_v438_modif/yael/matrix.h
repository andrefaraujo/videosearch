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

#ifndef __matrix_h
#define __matrix_h

/*---------------------------------------------------------------------------*/
/*! @addtogroup matrix
 *  @{
 */


/*! @defgroup matrix
 * Matrix functions
 *
 * All matrices are stored in column-major order (like Fortran and Matlab) and
 * indexed from 0 (like C, unlike Fortran). The declaration:
 *
 *     a(m, n) 
 * 
 * means that element a(i, j) is accessed with a[ i * m + j ] where
 *
 *     0 <= i < m and 0 <= j < n
 *
 * WARNING some matrix functions assume row-major storage! (noted with RM) 
 */



/*---------------------------------------------------------------------------*/
/* Standard operations                                                       */
/*---------------------------------------------------------------------------*/

/*! Allocate a new nrow x ncol matrix */
float *fmat_new (int nrow, int ncol);
float *fmat_new_0 (int nrow, int ncol);



/*!  Matrix multiplication. 
 * 
 * This function maps to the BLAS sgemm, assuming all matrices are packed
 * 
 * @param left(m,k)   left operand
 * @param right(k,n)  right perand
 * @param result(m,n) result matrix
 * @param m           nb of rows of left matrix and of result
 * @param n           nb of columns of right matrix and of result
 * @param k           nb of columns of left matrix and nb of rows of right matrix
 * @param transp transp[0] (resp. transp[1]) should be set to 'N' if
 *               the left (resp. right) matrix is in column-major (Fortran) 
 *               order and to 'T' if it is row-major order (C order). 
 *               The result is always in column-major order
 */

void fmat_mul_full(const float *left, const float *right,
                   int m, int n, int k,
                   const char *transp,
                   float *result);

/*! same as fmat_mul_full, matrices may be non-packed (yes, this is close to sgemm) */
void fmat_mul_full_nonpacked(const float *left, const float *right,
                             int m, int n, int k,
                             const char *transp,
                             int ld_left, int ld_right, 
                             float *result,
                             int ld_result);



/*!  same as fmat_mul, allocates result
 */
float* fmat_new_mul_full(const float *left, const float *right,
                         int m, int n, int k,
                         const char *transp);

/*! same as fmat_mul_full, all in standard order */
void fmat_mul (const float *left, const float *right, int m, int n, int k, float *mout);

/*! same as fmat_mul_full, left(k,m) transposed */
void fmat_mul_tl (const float *left, const float *right, int m, int n, int k, float *mout);

/*! same as fmat_mul_full, right(n,k) transposed */
void fmat_mul_tr (const float *left, const float *right, int m, int n, int k, float *mout);

/*! same as fmat_mul_full, left(k,m) and right(n,k) transposed */
void fmat_mul_tlr (const float *left, const float *right, int m, int n, int k, float *mout);


float* fmat_new_mul (const float *left, const float *right, int m, int n, int k);
float* fmat_new_mul_tl (const float *left, const float *right, int m, int n, int k);
float* fmat_new_mul_tr (const float *left, const float *right, int m, int n, int k);
float* fmat_new_mul_tlr (const float *left, const float *right, int m, int n, int k);



/*! solve the linear squares system a*x = b with n unknowns and m equations. 
 *
 * @param m      number of unkowns
 * @param n      number of equations
 * @param a(m,n) transposed matrix of the system 
 * @param b(n)   right-hand side of the equation 
 * @param x(m)   solution 
 * 
 * @return 0 if ok, else an error code (see sgels doc)
 */
int fmat_solve_ls_t(int m, int n, const float *a, const float *b, float *x);


/*! display the matrix in matlab-parsable format */
void fmat_print (const float *a, int nrow, int ncol);


/*! same as fmat_print but matrix is in row-major order */
void fmat_print_tranposed (const float *a, int nrow, int ncol);



/*---------------------------------------------------------------------------*/
/* Matrix manipulation functions                                             */
/*---------------------------------------------------------------------------*/

/*! Extract a submatrix.
 *  
 * @param a        the matrix (at least nrow by ncol)
 * @param nrow     nb of rows of input matrix
 * @param nrow_out nb of rows of output matrix
 * @param ncol     nb of columns of output matrix 
 * @return         the extracted submatrix 
 */

float *fmat_get_submatrix (const float *a, int nrow, 
                           int nrow_out,
                           int ncol);

int *imat_get_submatrix (const int *a, int nrow, 
                           int nrow_out,
                           int ncol);


/*! return the submatrix defined by a list of columns  */
float *fmat_new_get_columns (const float *a, int nrow, int ncolout, const int *cols);

void fmat_get_columns (const float *a, int d, int ncolout, const int *cols, float *out);

/*! return the matrix defined by 

out(i, j) = a(rows[i], cols[j])

*/
void fmat_get_rows_cols(const float *a, int d, 
                          int n_row, const int *rows, 
                          int n_col, const int *cols, 
                          float *out);

void fmat_shuffle_columns(float *a, int nrow, int ncol);

/*! produce a vector by taking a particular row of a matrix */
float *fmat_new_get_row (const float *a, int nrow, int ncol, int row);

/*!  produce a matrix composed of the rows indicated by the vector rows */
float *fmat_new_get_rows (const float *a, int d, int n, int nrowout, const int *rows);

/*! per-column sum of matrix elements. Output is a vector of length ncol */
void fmat_sum_columns (const float * a, int nrow, int ncol, float * sums);
float *fmat_new_sum_columns (const float *a, int nrow, int ncol);

/*! per-row sum of matrix elements */
void fmat_sum_rows (const float * a, int nrow, int ncol, float * sums);
float *fmat_new_sum_rows (const float *a, int nrow, int ncol);


float *fmat_new_vstack(const float *a,int da,
                       const float *b,int db,
                       int n);


/*! Matrix transposition
 * 
 * @param a         the matrix (nrow by ncol
 * @param ncol      number of columns of original matrix
 * @param nrow      number of rows of original matrix
 * @return          transposed copy of the matrix (and void for the inplace version)
 */

float *fmat_new_transp (const float *a, int ncol, int nrow);

/*! RM 
 * a is ncol-by-nrow
 * accu is k-by-k
 *
 * for i=0..ncol-1,j=0..nrow-1, do 
 *    accu(row_assign[i],col_assign[j]) += a(i,j)
 *
 */ 
void fmat_splat_separable(const float *a,int nrow,int ncol,
                          const int *row_assign,const int *col_assign,
                          int k,
                          float *accu); 

/*! Accumulate columns depending on assignment
 *
 * a(nrow, ncol)
 * accu(nrow, ?) 
 *
 * for i=0..ncol-1 do 
 *    accu(:, assign[i]) += a(:, i)
 *
 */ 
void fmat_splat_separable_1D(const float *a,int nrow,int ncol,
                             const int *assign,
                             float *accu); 




int *imat_joint_histogram(int n,int k,int *row_assign,int *col_assign);


/*! removes 0-filled columns of a matrix. Returns new number of columns */
int fmat_remove_0_columns(float *a, int d, int n); 

/*! replaces each column with colmn * (norm2sqr of column) ^ pw */
void fmat_normalize_columns_l2sqr_pow(float *a, int d, int n, float pw); 

/*! replaces each column with colmn * (norm2sqr of column + eps) ^ pw */
void fmat_normalize_columns_l2sqr_pow_robust(float *a, int d, int n, float pw, float eps);

/*---------------------------------------------------------------------------*/
/* Special matrices                                                          */
/*---------------------------------------------------------------------------*/


/*! RM  produce a new matrix of size nrow x ncol, filled with gaussian values */
float * fmat_new_rand_gauss (int nrow, int ncol);

/*! produce a random orthogonal basis matrix of size d*d */
float *random_orthogonal_basis (int d);

/*! Construct a Hadamard matrix of dimension d using the Sylvester construction.
   d should be a power of 2 */
float * hadamard (int d);


/*---------------------------------------------------------------------------*/
/* Statistical matrix operations                                             */
/*---------------------------------------------------------------------------*/

/* compute average of v matrix columns, subtract it to v and return average */
float *fmat_center_columns(int d,int n,float *v); 

/* subtract a vector from all columns of a matrix, m_i := m_i - avg*/
void fmat_subtract_from_columns(int d,int n,float *m,const float *avg);

/* add a vector from all columns of a matrix, m_i := m_i + avg*/
void fmat_add_to_columns(int d,int n,float *m,const float *avg);

/* reverse: m_i := avg - m_i */
void fmat_rev_subtract_from_columns(int d,int n,float *m,const float *avg);










/*! Compute covariance of a set of vectors
 * 
 * @param v(d,n)  vectors to compute covariance
 * @param avg(d)  on output, average vector (can be NULL)
 * @param assume_centered  assumes the data is centered (avg not used)
 * 
 * @return (d,d)  covariance matrix
 */
float *fmat_new_covariance (int d, int n, const float *v,
                        float *avg,int assume_centered);


/*! Perform the Principal Component Analysis of a set of vectors
 *
 * @param v(d,n)  vectors to perform the PCA on. The vectors are assumed to be centered already!
 * @param singvals(d)  corresponding singular values (may be NULL)
 *
 * @return (d,d) matrix of eigenvectors (column-stored). 
 */
float *fmat_new_pca(int d,int n,const float *v,
                    float *singvals); 

/*! same as fmat_pca, but return only a few vectors
 *
 * @param v(d,n)  vectors to perform the PCA on. The vectors are assumed to be centered already!
 * @param singvals(nev)  corresponding singular values (may be NULL)
 *
 * @return (d,nev) matrix of eigenvectors. To transform a
 *                 vector a low-dimension space, multiply by the d2<nev first lines of the matrix
 */

float *fmat_new_pca_part(int d,int n,int nev,
                         const float *v,float *singvals); 


/*! Compute SVD decomposition of a matrix: a = u * diag(singvals) * v'
 *
 * @param d            nb of rows of matrix a
 * @param n            nb of columns of matrix a
 * @param ns           nb of singular values to compute 
 * @param a(d,n)       matrix to compute singular vals for
 * @param singvals(ns) output singular values (may be NULL)
 * @param u(d,ns)      left orthogonal matrix (may be NULL)
 * @param v(n,ns)      right orthogonal matrix (may be NULL)
 */
int fmat_svd_partial(int d, int n, int ns, const float *a,
                     float *singvals, float *u, float *v); 

/*! with additionnal options
 */
int fmat_svd_partial_full(int n,int m,int nev,const float *a,int a_transposed,
                          float *s,float *vout,float *uout,int nt);


/*! Compute the PCA eigenvalues and eigenvectors from covariance matrix */
float *fmat_new_pca_from_covariance(int d, const float *cov, float *singvals);
void fmat_pca_from_covariance(int d,const float *cov, float *singvals, float * pcamat);


/* Online version of PCA */

/* a structure to perform an online PCA */
struct pca_online_s {
  int n;            /* number of samples used for the pca. A negative values means that 
		       the normalization of mu and cov has already been performed */
  int d;
  float * mu;
  float * cov;
  float * eigvec;
  float * eigval;
};

typedef struct pca_online_s pca_online_t;

/*! Construct the online PCA structure */
pca_online_t * pca_online_new (int d);

/*! Free memory associated with the online PCA structure */
void pca_online_delete (struct pca_online_s * pca);

/*! Accumulate information for PCA for n input vectors */
void pca_online_accu (struct pca_online_s * pca, const float * v, long n);

/*! compute the mean and covariance matrix */
void pca_online_cov (struct pca_online_s * pca);

/*! Online PCA: compute the mean and the eigenvectors.
   Also compute the covariance matrix if not already done. 
   The output is stored in the structure itself  */
void pca_online_complete (struct pca_online_s * pca);

/*! Same function as pca_online_complete but for compute only the eigenvectors 
  associated with the nev largest eigenvalues */
void pca_online_complete_part (struct pca_online_s * pca, int nev);


/*! Project some vectors according to a PCA structure */
void pca_online_project (const pca_online_t * pca, const float * v, float * vo, int d, long n, int dout);



/*---------------------------------------------------------------------------*/
/*! @} */
/*---------------------------------------------------------------------------*/

#endif 

/*---------------------------------------------------------------------------*/


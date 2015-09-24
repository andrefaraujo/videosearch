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

#ifndef __eigs_h
#define __eigs_h


/*---------------------------------------------------------------------------*/
/*! @addtogroup linearalgebra
 *  @{  */

/*! Compute the eigenvalues and eigvectors of a symmetric matrix m
  @param d           dimension of the square matrix m
  @param m(d,d)      the matrix (first elements are the first row)
  @param eigval(d)   on output the eigenvalues (unsorted)
  @param eigvec(d,d) on output, eigenvector j is eigvec(:,j)

  @return            =0 for success, else an error code (see info in lapack's dsygv documentation)
  
  the vectors eigval and eigvec must be allocated externally
*/
int eigs_sym (int d, const float * m, float * eigval, float * eigvec);


/*! Solve a generalized eigenvector problem */
int geigs_sym (int d, const float * a, const float * b, float * eigval, 
	       float * eigvec);


/*! Re-ordering of the eigenvalues and eigenvectors for a given criterion
   @param criterion equal to 0 for ascending order, descending otherwise
*/
void eigs_reorder (int d, float * eigval, float * eigvec, int criterion);


/*! same as eigs_sym, but returns only part of the vectors 
 * 
 * @param nev           nb of eigenvectors/values to return
 * @param eigval(nev)   the n eigenvalues
 * @param eigvec(d,nev) eigenvector j is eigvec(:,j)
 *
 * @return =0 for success, else an error code (see info in ssaupd or ierr in sseupd from arpack's documentation)
 */
int eigs_sym_part (int d, const float * m, int nev, float * eigval, float * eigvec);


/*! thin wrapper around the partial eigenvalue Arpack function */
typedef struct arpack_eigs_t arpack_eigs_t; 


/*! begin partial eigenvalue computation -- user should have a matrix multiplication function at hand
 *
 * @param n           dimension of the square matrix
 * @param nev         nb of eigenvectors/values to return
 */
arpack_eigs_t *arpack_eigs_begin(int n,int nev); 

/*! one iteration
 * @param x    *x_out is the array that should be multiplied (size n)
 * @param y    *y_out is result of the multiplication (size n)
 * @return     >0 compute y := A * x
               0: stop iteration
               <0, error (call arpack_eigs_end for cleanup)
 */
int arpack_eigs_step(arpack_eigs_t *,
                     float **x_out, float **y_out); 

/*! result and cleanup 
 * @param sout   eigenvalues
 * @param vout   eigenvectors
 * @return       nb of filled-in eigenvals and eigenvecs (may be below nev if some did not converge)
 */
int arpack_eigs_end(arpack_eigs_t *,
                    float * sout, float * vout); 







/*! @} */
#endif

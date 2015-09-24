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

#ifndef __binheap_h
#define __binheap_h

#include <stdlib.h>

/*---------------------------------------------------------------------------*/
/*! @addtogroup binheap
 *  @{
 */


/*! @defgroup binheap
  This structure is used, in particular, to find the maxk smallest
  elements of a possibly unsized stream of values. 
*/

/*! Binary heap used as a maxheap. 
  Element (label[1],val[1]) always contains the maximum value of the binheap. 
*/
struct fbinheap_s {
  float * val;     /*!< valid values are val[1] to val[k] */
  int * label;     /*!< idem for labels */
  int k;           /*!< number of elements stored  */
  int maxk;        /*!< maximum number of elements */
};

typedef struct fbinheap_s fbinheap_t;

/*! create the maxheap structure for maxk elements (maximum)
 @param maxk maximum number of elements to be stored in the heap*/
struct fbinheap_s * fbinheap_new (int maxk);

/*! return the size of a maxheap structure 
 @param maxk the maximum number of elements that the structure will receive */
size_t fbinheap_sizeof (int maxk); 

/*! A binheap can be stored in an externally allocated memory area 
  of fbinheap_sizeof(maxk) bytes. The fbinheap_init() function is used 
  to initialize this memory area */
void fbinheap_init (fbinheap_t *bh, int maxk);

/*! free allocated memory */
void fbinheap_delete (fbinheap_t * bh);

/*! remove all the elements from the heap */
void fbinheap_reset (fbinheap_t *bh);

/*! insert an element on the heap (if the value val is small enough) */
void fbinheap_add (fbinheap_t * bh, int label, float val);

/*! remove largest value from binheap (low-level access!) */
void fbinheap_pop (fbinheap_t * bh);

/*! add n elements on the heap (the values are added only if they 
are small enough compared to the other elements)

  @param bh the maxheap structure
  @param n the number of elements to be added
  @param labels The identifiers for the values to be added
  @param v the set of vectors to be added
 */
void fbinheap_addn (fbinheap_t * bh, int n, const int * labels, const float * v);

/*! add n elements on the heap, using the set of labels starting at label0  */
void fbinheap_addn_label_range (fbinheap_t * bh, int n, int label0, const float * v);

/*! output the labels in increasing order of associated values 
  @param bh the maxheap structure
  @pram perm the array that receive the output permutation order (pre-allocated)
*/
void fbinheap_sort_labels (fbinheap_t * bh, int * perm);

/*! output the sorted values */
void fbinheap_sort_values (fbinheap_t * bh, float * v);

/*! output both sorted results: labels and corresponding values  */
void fbinheap_sort (fbinheap_t * bh, int * labels, float *v);

/*! sort by increasing labels, ouptput sorted labels & associated values */
void fbinheap_sort_per_labels (fbinheap_t * bh, int * labels, float *v);

/*! show the heap content */
void fbinheap_display (fbinheap_t * bh);


/*! @} */











/*! idem for arbitrary data. The size of the labels is provided on construction */

struct abinheap_s {
  float * val;     /*!< valid values are val[1] to val[k] */
  void * label;
  int k;           /*!< number of elements stored  */
  int maxk;        /*!< maximum number of elements */
  int labelsize;
};

typedef struct abinheap_s abinheap_t;

/*! create the maxheap structure for maxk elements (maximum)
 @param maxk maximum number of elements to be stored in the heap*/
struct abinheap_s * abinheap_new (int maxk, int labelsize);

/*! return the size of a maxheap structure 
 @param maxk the maximum number of elements that the structure will receive */
size_t abinheap_sizeof (int maxk, int labelsize); 

/*! A binheap can be stored in an externally allocated memory area 
  of abinheap_sizeof(maxk) bytes. The abinheap_init() function is used 
  to initialize this memory area */
void abinheap_init (abinheap_t *bh, int maxk, int labelsize);

/*! free allocated memory */
void abinheap_delete (abinheap_t * bh);

/*! remove all the elements from the heap */
void abinheap_reset (abinheap_t *bh);

/*! insert an element on the heap (if the value val is small enough) */
void abinheap_add (abinheap_t * bh, void *label, float val);

/*! remove largest value from binheap (low-level access!) */
void abinheap_pop (abinheap_t * bh);

/*! output both sorted results: labels and corresponding values  */
void abinheap_sort (abinheap_t * bh, void *labels, float *v);

/*! get the i'th label (0-based) */
void *abinheap_get_label (abinheap_t *bh, int i);


#endif

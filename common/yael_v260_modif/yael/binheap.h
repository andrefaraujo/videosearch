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


#endif

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

#ifndef SORTING_H_INCLUDED
#define SORTING_H_INCLUDED

/*! @addtogroup sorting
 *  @{  */

/*! @defgroup: sorting  
 *
 * Various sorting functions + a few simple array functions that can
 * be called from python efficiently */

/*! Find the maximum elements of an array.
 * 
 * @param v(n)      array to search
 * @param maxes(k)  largest values, on output: 
 *
 *    v[maxes[0]] >= v[maxes[1]] >= ... >= v[maxes[k-1]] >= v[i] for all i not in maxes.
 */
void fvec_k_max(const float *v, int n, int *maxes, int k);


/*! Find the minimum elements of an array.
 * See find_k_max. 
 */
void fvec_k_min(const float *v, int n, int *mins, int k);


/*! finds the ranks of a few values in a large set.
 *
 * Finds the ranks the values would have if the set was sorted by
 * *decreasing* order.
 *
 * @param tab(n)         unsorted table in which ranks are to be found
 * @param vals(nval)     values whose ranks are to be found 
 * @param minranks(nval) minimum rank of each value (may be NULL) 
 * @param maxranks(nval) maximum rank of each value + 1 (may be NULL) 
 * 
 * On ouput, for value 0 <= i < nval, 
 * - minranks[i]-1 is the highest index of values > vals[i]
 * - maxranks[i] is the lowest index of values < vals[i]
 * - if vals[i] is in the array, elements minranks[i] to maxranks[i]-1 have value vals[i]
 * 
 * The algorithm is in O(n * log(nval)). If nval is larger, it is 
 * more efficient to sort the array.
 */ 
void fvec_ranks_of(const float *tab,int n,
                     const float *vals,int nval,
                     int *minranks,int *maxranks);

/*! idem but ranks in increasing array */
void fvec_ranks_inc_of(const float *tab, int n,
                         const float *vals, int nval,
                         int *minranks, int *maxranks);


/*---------------------------------------------------------------------------*/
/* Simple index functions (useful to call from C)                            */
/*---------------------------------------------------------------------------*/

/*! Replace ilabels[i] with the location of ilabels[i] in the table labels.
 *
 * @param labels(nres)    
 * @param ilabels(nilabels) 
 * 
 * On output ilabels is modified: 
 * 
 *   labels[ilabels_out[i]] = ilabels[i] for 0<=i<nilabels or -1 if there is none
 */
void find_labels (const int *labels, int nres, int *ilabels, int nilabels);

/*! return the smallest value of a vector */
float fvec_min (const float *f, long n);
int ivec_min (const int * f, long n);


/*! return the largest value of a vector */
float fvec_max(const float *f, long n);
int ivec_max (const int *f, long n);

/*! return the position of the smallest element of a vector.
  First position in case of ties, n should be >0. */
int fvec_arg_min (const float *f, long n);

/*! return the position of the largest elements of a vector.
  First position in case of ties, n should be >0. */
int fvec_arg_max (const float *f, long n);


/*! computes the median of a float array. Array modified on output! */
float fvec_median (float *f, int n);

/* computes the median of a float array (without modifying this array) */
float fvec_median_const (const float *f, int n);


/*! find quantile. 
 *
 * @return value v such that q elements are <= v
 */
float fvec_quantile (float *f,int n,int q);


/*! in-place sort */
void ivec_sort (int *tab, int n);

/*! return permutation to sort an array. 
 *
 * @param tab(n)    table to sort
 * @param perm(n)   output permutation that sorts table
 * 
 * On output,  
 *
 *      tab[perm[0]] <= tab[perm[1]] <= ... <= tab[perm[n-1]]
 *
 * Is stable. 
 */
void ivec_sort_index (const int *tab, int n, int *perm);

/*! fill-in iperm so that iperm[perm[i]]=i for i=0..n-1 */
void ivec_invert_perm (const int *perm, int n, int *iperm); 


/*! in-place sort */
void fvec_sort (float *v, int n);

/*! in-place sort for several vectors */
void fvecs_sort (float *v, int d, int n);

/*! return permutation to sort an array. See ivec_sort_index. */
void fvec_sort_index (const float *tab, int n, int *perm);

/*! Apply a permutation to a vector. The permutation is 
 * typically generated using the ivec_sort_index function. In that 
 *  case the function outputs a sorted array. 
 */
void ivec_sort_by_permutation (int * v, const int * order, int n);

void fvec_sort_by_permutation (float * v, const int * order, int n);


/************ operations on sorted int arrays */

/*! count occurrences of val in sorted vector */
int ivec_sorted_count_occurrences (const int *v, int n, int val);

/*! find index of highest value <= val (=-1 if all values are > val) */
int ivec_sorted_find (const int *v, int n, int val);

/*! count the number of distinct values in the input fvector */
int ivec_sorted_count_unique (const int *v, int n);

/*! count the number of occurrences of several values */
int ivec_sorted_count_occurrences_multiple (const int *v, int n,
                                            const int *vals, int nval);


/*! merge several sorted sets
 * 
 * There are k sets. Set 0 <= i < k has size sizes[i]. Element j of
 * set i is vals[i][j], and the arbitrary associated index is
 * lists[i][j]. The set is ordered, so on input 
 *
 *    vals[i][0] <= vals[i][1] <= ... <= vals[i][sizes[i]-1]
 * 
 * On ouput, *labels_out and *vals_out contain an ordered set with all values.
 *
 * @return total number of elements (=sum(sizes[i],i=0..k-1))
 */
int merge_ordered_sets (const int **labels, const float **vals,
                        const int *sizes, int k,
                        int **labels_out, float **vals_out); 


/*! remove largest values from an array 
 *
 * Finds the smallest value m of vals, compresses array labels by
 * removing labels[i] for which vals[i] < m * ratio returns new size
 * of labels array. 
 *
 * NB that on output, vals[i] does not correspond to labels[i] any more!
 * 
 */
int compress_labels_by_disratio (int *labels, const float *vals, int n, 
				 float ratio); 


/*! @} */
#endif

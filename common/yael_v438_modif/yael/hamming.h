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

/* Hamming distances. The binary vector length should be a power of 8 */
#ifndef __hamming_h
#define __hamming_h

typedef unsigned char uint8;
typedef unsigned short uint16;
typedef unsigned int uint32;
typedef unsigned long long uint64;

typedef long long int64;


/* matching elements (those returned) */
typedef struct hammatch_s {
  int qid;        /* query id */
  int bid;        /* base id */
  uint16 score;   /* Hamming distance */
} hammatch_t;


/* Define individual Hamming distance for various sizes.
   ncodes is given in bytes, therefore the actual number of bits is 8*ncodes.
   The generic one is slow while optimization is available for specific sizes */
uint16 hamming (const uint8 *bs1, const uint8 * bs2, int ncodes);



/* Compute a set of Hamming distances between na and nb binary vectors */
void compute_hamming (uint16 * dis, const uint8 * a, const uint8 * b, 
                      int na, int nb, int ncodes);




/* Counting the number of matches or of cross-matches (with actually returning them)
   Useful to be used with function that assume pre-allocated memory                  */
void match_hamming_count (const uint8 * bs1, const uint8 * bs2, int n1, int n2, 
                          int ht, int ncodes, size_t * nptr);

void crossmatch_hamming_count (const uint8 * dbs, int n, int ht, 
                               int ncodes, size_t * nptr);


/* For 1 query signature, compute the hamming distance and report those below a given 
   threshold in a structure array */
void match_hamming_thres (const uint8 * bs1, const uint8 * bs2, 
                          int n1, int n2, int ht, int ncodes, size_t bufsize, 
                          hammatch_t ** hmptr, size_t * nptr);


/* The same but with pre-allocation (typically used with match_hamming_count) */
size_t match_hamming_thres_prealloc (const uint8 * bs1, const uint8 * bs2, 
                                     int n1, int n2, int ht, int ncodes, 
                                     int * idx, uint16 * hams);

                                   /* Compute all cross-distances between two sets of binary vectors */
void crossmatch_hamming (const uint8 * dbs, long n, int ht, int ncodes, 
                         long bufsize, hammatch_t ** hmptr, size_t * nptr);

/* alternative variant with pre-allocated external memory.
   return number of elements for safety check. 
   Typical usage is to first invoke crossmatch_hamming_count, allocate memory,
   and then invoke crossmatch_hamming_prealloc */

size_t crossmatch_hamming_prealloc (const uint8 * dbs, long n, int ht, int ncodes,  
                                    int * idx, uint16 * hams);

/* Threaded versions, when OpenMP is available */
#ifdef _OPENMP
void compute_hamming_thread (uint16 * dis, const uint8 * a, const uint8 * b, 
                             int na, int nb, int ncodes);

size_t match_hamming_thres_nt (const uint8 * bs1, const uint8 * bs2, int n1, int n2, 
                              int ht, int ncodes, int nt, int ** keys, uint16 ** ham);
#endif /* _OPENMP */


#endif /* __hamming_h */

 
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

/* This code was written by Herve Jegou. Contact: herve.jegou@inria.fr  */
/* Last change: June 1st, 2010                                          */
/* This software is governed by the CeCILL license under French law and */
/* abiding by the rules of distribution of free software.               */
/* See http://www.cecill.info/licences.en.html                          */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hamming.h"

/* If SSE4.2 is available, use the specific processor instructions */
#ifdef __SSE4_2__
#include <nmmintrin.h>
#define hamming_32(pa,pb) _mm_popcnt_u32((*((const uint32 *) (pa)) ^ *((const uint32 *) (pb))))
#define hamming_64(pa,pb) _mm_popcnt_u64((*((const uint64 *) (pa)) ^ *((const uint64 *) (pb))))
#endif

#define hamming_128(a,b)  (hamming_64((const uint64 *) (a),(const uint64 *) (b))+hamming_64(((const uint64 *) (a)) + 1, ((const uint64 *) (b)) + 1))

#define MIN(a,b) ((a)>(b) ? (b) : (a))

/* Define the Hamming distance by selecting the most appropriate function,
 using the generic version as a backup */


/* the slice size is set to avoid testing the buffer size too often */
#define HAMMATCH_SLICESIZE 16

/* For functions that compute distances by blocks */
#define HAM_BLOCKSIZE  128

/* geometric re-allocation: add a constant size plus a relative 50% of additional memory */
#define HAMMATCH_REALLOC_NEWSIZE(oldsize) (HAMMATCH_SLICESIZE+((oldsize * 5) / 4))


static uint16 uint8_nbones[256] = {
    0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
    4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};



/*-------------------------------------------------------*/
/* Elementary Hamming distance computation: unoptimized  */

uint16 hamming (const uint8 *bs1, const uint8 * bs2, int ncodes)
{
  int i;
  uint16 ham = 0;

  for (i = 0; i < ncodes ; i++) {
    ham += uint8_nbones[*bs1 ^ *bs2];
    bs1++;
    bs2++;
  }

  return ham;
}


#ifndef __SSE4_2__
#warning "SSE4.2 NOT USED FOR HAMMING DISTANCE COMPUTATION. Consider adding -msse4!"

static uint16 hamming_32 (const uint32 * bs1, const uint32 * bs2)
{
  uint16 ham = 0;
  uint32 diff = ((*bs1) ^ (*bs2));

  ham = uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff];
  return ham;
}


static uint16 hamming_64 (const uint64 * bs1, const uint64 * bs2)
{
  uint16 ham = 0;
  uint64 diff = ((*bs1) ^ (*bs2));

  ham = uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff & 255];
  diff >>= 8;
  ham += uint8_nbones[diff];

  return ham;
}

#endif


/*-------------------------------------------------------*/
/* Compute a set of Hamming distances                    */
static void compute_hamming_32 (uint16 * dis, const uint32 * a, const uint32 * b, int na, int nb)
{
  int i, j;
  const uint32 * pb = (const uint32 *) b;
  for (j = 0 ; j < nb ; j++) {
    const uint32 * pa = (const uint32 *) a;
    for (i = 0 ; i < na ; i++) {
      *dis = hamming_32 (pa, pb);
      pa++;
      dis++;
    }
    pb++;
  }
}


static void compute_hamming_64 (uint16 * dis, const uint64 * a, const uint64 * b, int na, int nb)
{
  int i, j;
  const uint64 * pb = (const uint64 *) b;
  for (j = 0 ; j < nb ; j++) {
    const uint64 * pa = (const uint64 *) a;
    for (i = 0 ; i < na ; i++) {
      *dis = hamming_64 (pa, pb);
      pa++;
      dis++;
    }
    pb++;
  }
}


static void compute_hamming_128 (uint16 * dis, const uint64 * a, const uint64 * b, int na, int nb)
{
  int i, j;
  const uint64 * pb = (const uint64 *) b;
  for (j = 0 ; j < nb ; j++) {
    const uint64 * pa = (const uint64 *) a;
    for (i = 0 ; i < na ; i++) {
      *dis = hamming_128 ((const uint64 *) pa, (const uint64 *) pb);
      pa += 2;
      dis++;
    }
    pb += 2;
  }
}


void compute_hamming (uint16 * dis, const uint8 * a, const uint8 * b,
                      int na, int nb, int ncodes)
{
  switch (ncodes) {
    case 4:  compute_hamming_32 (dis, (const uint32 *) a, (const uint32 *) b, na, nb);  return;
    case 8:  compute_hamming_64 (dis, (const uint64 *) a, (const uint64 *) b, na, nb);  return;
    case 16: compute_hamming_128 (dis, (const uint64 *) a, (const uint64 *) b, na, nb); return;
    default: fprintf (stderr, "# Warning: non-optimized version of compute_hamming\n");
  }

  int i, j;
  const uint8 * pb = b;
  for (j = 0 ; j < nb ; j++) {
    const uint8 * pa = a;
    for (i = 0 ; i < na ; i++) {
      *dis = hamming (pa, pb, ncodes);
      pa += ncodes;
      dis++;
    }
    pb += ncodes;
  }
}


/*-------------------------------------------------------*/
/* Count number of matches given a threshold            */
static void match_hamming_count_32 (const uint32 * bs1, const uint32 * bs2, int n1, int n2, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint32 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_32 (bs1, bs2) <= ht)
        posm++;
      bs2++;
    }
    bs1++;
  }
  *nptr = posm;
}


static void match_hamming_count_64 (const uint64 * bs1, const uint64 * bs2, int n1, int n2, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint64 * bs1_ = bs1;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs1 = bs1_ + i;
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_64 (bs1, bs2) <= ht)
        posm++;
      bs2 += 1;
    }
    bs1 += 1;  /* next signature */
  }
  *nptr = posm;
}


static void match_hamming_count_128 (const uint64 * bs1, const uint64 * bs2, int n1, int n2, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_128 (bs1, bs2) <= ht)
        posm++;
      bs2 += 2;
    }
    bs1  += 2;  /* next signature */
  }
  *nptr = posm;
}


void match_hamming_count (const uint8 * bs1, const uint8 * bs2, int n1, int n2, int ht, int ncodes, size_t * nptr)
{
  size_t i, j, posm = 0;

  switch (ncodes) {
    case 4:  match_hamming_count_32 ((const uint32 *) bs1, (const uint32 *) bs2, n1, n2, ht, nptr);  return;
    case 8:  match_hamming_count_64 ((const uint64 *) bs1, (const uint64 *) bs2, n1, n2, ht, nptr);  return;
    case 16: match_hamming_count_128 ((const uint64 *) bs1, (const uint64 *) bs2, n1, n2, ht, nptr); return;
    default: fprintf (stderr, "# Warning: non-optimized version of match_hamming_count\n");
  }

  const uint8 * bs2_ = bs2;
  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming (bs1, bs2, ncodes) <= ht)
        posm++;
      bs2 += ncodes;
    }
    bs1  += ncodes;  /* next signature */
  }
  *nptr = posm;
}


/* Count number of cross-matches given a threshold            */
static void crossmatch_hamming_count_32 (const uint32 * dbs, int n, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint32 * bs1 = dbs;

  for (i = 0 ; i < n ; i++) {
    const uint32 * bs2 = bs1 + 1;
    for (j = i + 1 ; j < n ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_32 (bs1, bs2) <= ht)
        posm++;
      bs2++;
    }
    bs1++;
  }
  *nptr = posm;
}


static void crossmatch_hamming_count_64 (const uint64 * dbs, int n, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint64 * bs1 = dbs;

  for (i = 0 ; i < n ; i++) {
    const uint64 * bs2 = bs1 + 1;
    for (j = i + 1 ; j < n ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_64 (bs1, bs2) <= ht)
        posm++;
      bs2++;
    }
    bs1++;
  }
  *nptr = posm;
}


static void crossmatch_hamming_count_128 (const uint64 * dbs, int n, int ht, size_t * nptr)
{
  size_t i, j, posm = 0;
  const uint64 * bs1 = dbs;

  for (i = 0 ; i < n ; i++) {
    const uint64 * bs2 = bs1 + 2;

    for (j = i + 1 ; j < n ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming_128 (bs1, bs2) <= ht)
        posm++;
      bs2 += 2;
    }
    bs1 += 2;
  }
  *nptr = posm;
}


void crossmatch_hamming_count (const uint8 * dbs, int n, int ht, int ncodes, size_t * nptr)
{
  switch (ncodes) {
    case 4:  crossmatch_hamming_count_32 ((const uint32 *) dbs, n, ht, nptr);  return;
    case 8:  crossmatch_hamming_count_64 ((const uint64 *) dbs, n, ht, nptr);  return;
    case 16: crossmatch_hamming_count_128 ((const uint64 *) dbs, n, ht, nptr); return;
    default: fprintf (stderr, "# Warning: non-optimized version of crossmatch_hamming_count\n");
  }

  size_t i, j, posm = 0;
  const uint8 * bs1 = dbs;
  for (i = 0 ; i < n ; i++) {
    const uint8 * bs2 = bs1 + ncodes;

    for (j = i + 1 ; j < n ; j++) {
      /* collect the match only if this satisfies the threshold */
      if (hamming (bs1, bs2, ncodes) <= ht)
        posm++;
      bs2 += ncodes;
    }
    bs1  += ncodes;  /* next signature */
  }

  *nptr = posm;
}


/*-------------------------------------------------------*/
/* Return all matches given a threshold                  */

/* Compute hamming distance and report those below a given threshold in a structure array */
hammatch_t * hammatch_new (int n)
{
  return (hammatch_t *) malloc (n * sizeof (hammatch_t));
}


hammatch_t * hammatch_realloc (hammatch_t * m, int n)
{
  return (hammatch_t *) realloc (m, n * sizeof (hammatch_t));
}


static void match_hamming_thres_32 (const uint32 * bs1, const uint32 * bs2, int n1, int n2, int ht,
                                    size_t bufsize, hammatch_t ** hmptr, size_t * nptr)
{
  size_t i, j, posm = 0;
  uint16 h;
  *hmptr = hammatch_new (bufsize);
  hammatch_t * hm = *hmptr;
  const uint32 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      h = hamming_32 (bs1, bs2);

      if (h <= ht) {    /* Enough space to store another match ? */
        if (posm >= bufsize) {
          bufsize = HAMMATCH_REALLOC_NEWSIZE (bufsize);
          *hmptr = hammatch_realloc (*hmptr, bufsize);
          assert (*hmptr != NULL);
          hm = (*hmptr) + posm;
        }
        hm->qid = i;
        hm->bid = j;
        hm->score = h;
        hm++;
        posm++;
      }
      bs2++;  /* next signature */
    }
    bs1++;
  }
  *nptr = posm;
}


static void match_hamming_thres_64 (const uint64 * bs1, const uint64 * bs2, int n1, int n2, int ht,
                                    size_t bufsize, hammatch_t ** hmptr, size_t * nptr)
{
  size_t i, j, posm = 0;
  uint16 h;
  *hmptr = hammatch_new (bufsize);
  hammatch_t * hm = *hmptr;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      h = hamming_64 (bs1, bs2);

      if (h <= ht) {    /* Enough space to store another match ? */
        if (posm >= bufsize) {
          bufsize = HAMMATCH_REALLOC_NEWSIZE (bufsize);
          *hmptr = hammatch_realloc (*hmptr, bufsize);
          assert (*hmptr != NULL);
          hm = (*hmptr) + posm;
        }
        hm->qid = i;
        hm->bid = j;
        hm->score = h;
        hm++;
        posm++;
      }
      bs2++;  /* next signature */
    }
    bs1++;
  }
  *nptr = posm;
}


static void match_hamming_thres_128 (const uint64 * bs1, const uint64 * bs2, int n1, int n2, int ht,
                                     size_t bufsize, hammatch_t ** hmptr, size_t * nptr)
{
  size_t i, j, posm = 0;
  uint16 h;
  *hmptr = hammatch_new (bufsize);
  hammatch_t * hm = *hmptr;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      h = hamming_128 (bs1, bs2);

      if (h <= ht) {    /* Enough space to store another match ? */
        if (posm >= bufsize) {
          bufsize = HAMMATCH_REALLOC_NEWSIZE (bufsize);
          *hmptr = hammatch_realloc (*hmptr, bufsize);
          assert (*hmptr != NULL);
          hm = (*hmptr) + posm;
        }
        hm->qid = i;
        hm->bid = j;
        hm->score = h;
        hm++;
        posm++;
      }
      bs2 += 2;  /* next signature */
    }
    bs1 += 2;
  }
  *nptr = posm;
}


void match_hamming_thres (const uint8 * bs1, const uint8 * bs2,
                          int n1, int n2, int ht, int ncodes, size_t bufsize,
                          hammatch_t ** hmptr, size_t * nptr)
{
  switch (ncodes) {
    case 4:  match_hamming_thres_32 ((const uint32 *) bs1, (const uint32 *) bs2, n1, n2, ht, bufsize, hmptr, nptr);  return;
    case 8:  match_hamming_thres_64 ((const uint64 *) bs1, (const uint64 *) bs2, n1, n2, ht, bufsize, hmptr, nptr);  return;
    case 16: match_hamming_thres_128 ((const uint64 *) bs1, (const uint64 *) bs2, n1, n2, ht, bufsize, hmptr, nptr); return;
    default: fprintf (stderr, "# Warning: non-optimized version of match_hamming_thres\n");
  }

  size_t i, j, posm = 0;
  uint16 h;
  *hmptr = hammatch_new (bufsize);
  hammatch_t * hm = *hmptr;
  const uint8 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming (bs1, bs2, ncodes);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        if (posm >= bufsize) {
          bufsize = HAMMATCH_REALLOC_NEWSIZE (bufsize);
          *hmptr = hammatch_realloc (*hmptr, bufsize);
          assert (*hmptr != NULL);
          hm = (*hmptr) + posm;
        }
        hm->qid = i;
        hm->bid = j;
        hm->score = h;
        hm++;
        posm++;
      }
      bs2 += ncodes;  /* next signature */
    }
    bs1 += ncodes;
  }

  *nptr = posm;
}


static size_t match_hamming_thres_prealloc_32 (const uint32 * bs1,
                                               const uint32 * bs2,
                                               int n1, int n2, int ht,
                                               int * idx, uint16 * hams)
{
  size_t i, j, posm = 0;
  uint16 h;
  const uint32 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming_32 (bs1, bs2);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        *idx = i; idx++;
        *idx = j; idx++;

        *hams = h;
        hams++;
        posm++;
      }
      bs2++;  /* next signature */
    }
    bs1++;
  }
  return posm;
}


static size_t match_hamming_thres_prealloc_64 (const uint64 * bs1,
                                               const uint64 * bs2,
                                               int n1, int n2, const int ht,
                                               int * idx, uint16 * hams)
{
  size_t i, j, posm = 0;
  uint16 h;
  const uint64 * bs1_ = bs1;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs1 = bs1_ + i;
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming_64 (bs1, bs2);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        *idx = i; idx++;
        *idx = j; idx++;

        *hams = h;
        hams++;
        posm++;
      }
      bs2++;  /* next signature */
    }
  }
  return posm;
}

#ifdef NOTDEFINED
/* Blocked version -> not faster, not used */
 static size_t match_hamming_thres_prealloc_64 (const uint64 * bs1,
                                               const uint64 * bs2,
                                               const int n1, const int n2, const int ht,
                                               int * idx, uint16 * hams)
{
  size_t i, j, posm = 0, bli, blj;
  uint16 h;
  const uint64 * bs1_ = bs1;
  const uint64 * bs2_ = bs2;

  for (bli = 0 ; bli < n1 ; bli += HAM_BLOCKSIZE) {
    const size_t bli_end = MIN(bli+HAM_BLOCKSIZE, n1);

    for (blj = 0 ; blj < n2 ; blj += HAM_BLOCKSIZE) {
      const size_t blj_end = MIN(blj+HAM_BLOCKSIZE, n2);

      for (i = bli ; i < bli_end ; i++) {
        bs1 = bs1_ + i;
        bs2 = bs2_ + blj;
        for (j = blj ; j < blj_end ; j++) {
          /* Here perform the real work of computing the distance */
          h = hamming_64 (bs1, bs2);

          /* collect the match only if this satisfies the threshold */
          if (h <= ht) {
            /* Enough space to store another match ? */
            *idx = i; idx++;
            *idx = j; idx++;

            *hams = h;
            hams++;
            posm++;
          }
          bs2++;  /* next signature */
        }
        bs1++;
      }
    }
  }
  return posm;
}
#endif


static size_t match_hamming_thres_prealloc_128 (const uint64 * bs1,
                                                const uint64 * bs2,
                                                int n1, int n2, int ht,
                                                int * idx, uint16 * hams)
{
  size_t i, j, posm = 0;
  uint16 h;
  const uint64 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming_128 (bs1, bs2);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        *idx = i; idx++;
        *idx = j; idx++;

        *hams = h;
        hams++;
        posm++;
      }
      bs2+=2;  /* next signature */
    }
    bs1+=2;
  }
  return posm;
}


size_t match_hamming_thres_prealloc (const uint8 * bs1, const uint8 * bs2,
                                     int n1, int n2, int ht, int ncodes,
                                     int * idx, uint16 * hams)
{
  switch (ncodes) {
    case 4:  return match_hamming_thres_prealloc_32 ((const uint32 *) bs1,
              (const uint32 *) bs2, n1, n2, ht, idx, hams);
    case 8:  return match_hamming_thres_prealloc_64 ((const uint64 *) bs1,
              (const uint64 *) bs2, n1, n2, ht, idx, hams);
    case 16: return match_hamming_thres_prealloc_128 ((const uint64 *) bs1,
              (const uint64 *) bs2, n1, n2, ht, idx, hams);
    default: fprintf (stderr, "# Warning: non-optimized version of match_hamming_thres\n");
  }

  size_t i, j, posm = 0;
  uint16 h;
  const uint8 * bs2_ = bs2;

  for (i = 0 ; i < n1 ; i++) {
    bs2 = bs2_;
    for (j = 0 ; j < n2 ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming (bs1, bs2, ncodes);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        *idx = i; idx++;
        *idx = j; idx++;

        *hams = h;
        hams++;
        posm++;
      }
      bs2 += ncodes;  /* next signature */
    }
    bs1 += ncodes;
  }

  return posm;
}


void crossmatch_hamming (const uint8 * dbs, long n, int ht, int ncodes,
                    long bufsize, hammatch_t ** hmptr, size_t * nptr)
{
  size_t i, j, posm = 0;
  uint16 h;
  *hmptr = hammatch_new (bufsize);
  hammatch_t * hm = *hmptr;
  const uint8 * bs1 = dbs;

  for (i = 0 ; i < n ; i++) {
    const uint8 * bs2 = bs1 + ncodes;

    for (j = i + 1 ; j < n ; j++) {

      /* Here perform the real work of computing the distance */
      h = hamming (bs1, bs2, ncodes);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        if (posm >= bufsize) {
          bufsize = HAMMATCH_REALLOC_NEWSIZE (bufsize);
          *hmptr = hammatch_realloc (*hmptr, bufsize);
          assert (*hmptr != NULL);
          hm = (*hmptr) + posm;
        }

        hm->qid = i;
        hm->bid = j;
        hm->score = h;
        hm++;
        posm++;
      }
      bs2 += ncodes;
    }
    bs1  += ncodes;  /* next signature */
  }

  *nptr = posm;
}


size_t crossmatch_hamming_prealloc (const uint8 * dbs, long n, int ht,
                               int ncodes, int * idx, uint16 * hams)
{
  size_t i, j, posm = 0;
  uint16 h;
  const uint8 * bs1 = dbs;

  for (i = 0 ; i < n ; i++) {
    const uint8 * bs2 = bs1 + ncodes;

    for (j = i + 1 ; j < n ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming (bs1, bs2, ncodes);

      /* collect the match only if this satisfies the threshold */
      if (h <= ht) {
        /* Enough space to store another match ? */
        *idx = i; idx++;
        *idx = j; idx++;

        *hams = h;
        hams++;
        posm++;
      }
      bs2 += ncodes;
    }
    bs1 += ncodes;  /* next signature */
  }
  return posm;
}


/*-------------------------------------------*/
/* Threaded versions, if OpenMP is available */
#ifdef _OPENMP

#define HAMBLOCK 128
#define MIN(a,b) ((a)<(b) ? (a) : (b))

void compute_hamming_thread (uint16 * dis, const uint8 * a, const uint8 * b,
                             int na, int nb, int ncodes)
{
  size_t i, j;
#pragma omp parallel shared (dis, a, b, na, nb) private (i, j)
    {
#pragma omp for
      for (j = 0 ; j < nb ; j++)
	      for (i = 0 ; i < na ; i++)
	        dis[j * na + i] = hamming (a + i * ncodes, b + j * ncodes, ncodes);
    }
}


size_t match_hamming_thres_nt (const uint8 * bs1, const uint8 * bs2, int n1, int n2,
                              int ht, int ncodes, int nt, int ** keys, uint16 ** ham)
{
  size_t bl, nmatches;
  const int nblock1 = 1 + (n1 - 1) / HAMBLOCK;
  const int nblock2 = 1 + (n2 - 1) / HAMBLOCK;
  const int nblock = nblock1 * nblock2;

  size_t * blcount = malloc ((nblock + 1) * sizeof (*blcount));
  blcount[0] = 0;

  #pragma omp parallel for private(bl)
  for (bl = 0 ; bl < nblock ; bl++) {
    size_t bl1 = bl / nblock1;
    size_t bl2 = bl % nblock1;

    size_t s1 = bl1 * HAMBLOCK;
    size_t s2 = bl2 * HAMBLOCK;
    size_t nb1 = MIN(n1 - s1, HAMBLOCK);
    size_t nb2 = MIN(n2 - s2, HAMBLOCK);

    match_hamming_count (bs1 + s1 * ncodes, bs2 + s2 * ncodes,
                         nb1, nb2, ht, ncodes, blcount + bl + 1);
  }

  /* accumulate count to determine offset */
  nmatches = 0;
  for (bl = 1 ; bl <= nblock ; bl++) {
    if (blcount[bl] > 500)
      fprintf (stderr, "bl %ld -> %ld matches  (bl-1/cum = %ld)\n", bl-1, blcount[bl], blcount[bl-1]);
    blcount[bl] = blcount[bl-1] + blcount[bl];
  }
  nmatches = blcount[nblock];
  fprintf (stderr, "nmatches = %d\n", nmatches);

  *keys = malloc (nmatches * 2 * sizeof(**keys));
  *ham = malloc (nmatches * sizeof(**ham));

 #pragma omp parallel for private(bl)
    for (bl = 0 ; bl < nblock ; bl++) {
      size_t bl1 = bl / nblock1;
      size_t bl2 = bl % nblock1;

      size_t s1 = bl1 * HAMBLOCK;
      size_t s2 = bl2 * HAMBLOCK;
      size_t nb1 = MIN(n1 - s1, HAMBLOCK);
      size_t nb2 = MIN(n2 - s2, HAMBLOCK);

      match_hamming_thres_prealloc (bs1 + s1 * ncodes, bs2 + s2 * ncodes,
                                    nb1, nb2, ht, ncodes,
                                    (int*) (*keys) + blcount[bl] * 2,
                                    (uint16*) (*ham) + blcount[bl]);
    }

  free (blcount);
  return nmatches;
}

#endif /* _OPENMP */

#undef HAM_BLOCKSIZE

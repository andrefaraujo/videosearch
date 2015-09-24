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
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "ivf.h"
#include "hamming.h"


/* geometric re-allocation add a relative 50% of additional memory */
#define IVF_REALLOC_NEWSIZE(oldsize) (4+((oldsize * 3) / 2))


int ivf_checksum (const ivf_t * ivf)
{
  return (ivf->checksum == IVFCHECKSUM);
}


ivf_t *ivf_new (int k, int elemsize, int seg_size)
{
  int i;

  if (seg_size == 0)
    seg_size = DEFAULT_SEG_SIZE;

  ivf_t *ivf = (ivf_t *) malloc (sizeof (ivf_t));
  if (!ivf)
    return NULL;

  ivf->checksum = IVFCHECKSUM;
  ivf->k = k+1;  /* k+1 to be able to use the indexing convention of matlab */
  ivf->n = 0;
  ivf->elem_size = elemsize;

  ivf->nbelems = (int *) malloc (sizeof (int) * ivf->k);
  ivf->seg_size = (int *) malloc (sizeof (int) * ivf->k);
  ivf->ids = (int **) malloc (sizeof (*ivf->ids) * ivf->k);
  if (!ivf->ids) {
    ivf_delete (ivf);
    return NULL;
  }

  ivf->adat = (unsigned char **) malloc (sizeof (*ivf->adat) * ivf->k);
  if (!ivf->ids) {
    ivf_delete (ivf);
    return NULL;
  }

  /* a minimum allocated segment size by default equal to seg_size */
  for (i = 0; i < ivf->k ; i++) {
    ivf->seg_size[i] = seg_size;
    ivf->nbelems[i] = 0;
    ivf->ids[i] = (int *) malloc (sizeof (**ivf->ids) * ivf->seg_size[i]);
    if (!ivf->ids[i]) {
      ivf_delete (ivf);
      return NULL;
    }
  }

  /* store additional info only if elem_size > 0 */
  for (i = 0; i < ivf->k ; i++) {
    ivf->adat[i] = (unsigned char *) malloc (ivf->elem_size * ivf->seg_size[i]);
    if (!ivf->adat[i]) {
      ivf_delete (ivf);
      return NULL;
    }
  }

  return ivf;
}


void ivf_delete (ivf_t * ivf)
{
  int i;

  if (ivf == NULL)
    return;

  if (ivf->ids)
    for (i = 0; i < ivf->k; i++)
      free (ivf->ids[i]);

  if (ivf->adat)
    for (i = 0; i < ivf->k; i++)
      free (ivf->adat[i]);

  free (ivf->nbelems);
  free (ivf->seg_size);
  free (ivf->ids);
  free (ivf->adat);
  free (ivf);
}


void ivf_addn (ivf_t * ivf, const int * ids, const int * keys, 
		 const unsigned char * adat, const int n)
{
  int i, w, j, elem_size = ivf->elem_size;
  for (i = 0; i < n; i++) {
    w = keys[i];
    if (! (w >= 0 && w < ivf->k)) {
      fprintf (stderr, "# Invalid key id : %d (range must be %d-%d)\n", 
	       w, 0, ivf->k-1);
      return;
    }

    /* First check if a realloc is required or not */
    if (ivf->nbelems[w] + 1 == ivf->seg_size[w]) {       /* -> YES, it is */
      ivf->seg_size[w] = IVF_REALLOC_NEWSIZE (ivf->seg_size[w]);

      ivf->ids[w] = (int *) realloc (ivf->ids[w], (int) ivf->seg_size[w]
			     * sizeof (*ivf->ids[w]));
      assert (ivf->ids[w]);

      if (elem_size == 0)
        continue;

      ivf->adat[w] = (unsigned char *)realloc (ivf->adat[w], (int) ivf->seg_size[w] * elem_size);
      assert (ivf->adat[w]);
    }
    j = ivf->nbelems[w];

    /* Store the id of the image or vector */
    ivf->ids[w][j] = ids[i];

    if (elem_size > 0)
      memcpy (ivf->adat[w] + j * elem_size, adat + i * elem_size, elem_size);

    /* update the number of elements */
    ivf->nbelems[w]++;
  }
  ivf->n += n;
}


void ivf_add (ivf_t * ivf, int key, int id, const unsigned char * val)
{
  ivf_addn (ivf, &id, &key, val, 1);
}


int ivf_get_nb_elems (const ivf_t * ivf, int key)
{
  assert (key >= 0 && key < ivf->k);
  return ivf->nbelems[key];
}


int * ivf_get_ids (const ivf_t * ivf, int key)
{
  assert (key >= 0 && key < ivf->k);
  return ivf->ids[key];
}


unsigned char * ivf_get_vals (const ivf_t * ivf, int key)
{
  assert (key >= 0 && key < ivf->k);
  return ivf->adat[key];
}


unsigned char * ivf_find_vals (const ivf_t * ivf, int * keys, int * ids, int n)
{
  unsigned char * dat = (unsigned char *) malloc (ivf->elem_size * n);
	int j = 0, i, f;
	
	while (j<n) {
	   f = 0;
	 	for (i = 0; i < ivf->nbelems[keys[j]]; i++) {
			if (ivf->ids[keys[j]][i] == ids[j])	{
				f = 1;
				memcpy(dat + j * ivf->elem_size, ivf->adat[keys[j]] + i * ivf->elem_size, ivf->elem_size);
				j++;
				if (j>=n)
					break;
				if (keys[j] != keys[j-1])
					break;
			}		
		}
		if (f!=1)  /* id not found in the given key list, wrong input */
		{
			free (dat);
			return NULL;
		}

	}

  return dat;
}


void ivf_display (const ivf_t * ivf)
{
  int i, j;
  printf ("Nb lists  %d\n", ivf->k);

  /* for each segment, display the contents */
  for (i = 0; i < ivf->k; i++) {

    if (ivf->nbelems[i] > 0)
      fprintf (stdout, "[ List %d ] %d elements (seg_size: %d)\n", i,
               ivf->nbelems[i], ivf->seg_size[i]);
    else continue;

    for (j = 0; j < ivf->nbelems[i]; j++)
      printf ("%8d / ", (int) ivf->ids[i][j]);
    printf ("\n");
  }
}


/* Count the total number of elements (descriptors) in the inverted file */
int ivf_count_nbelems (const ivf_t * ivf)
{
  int tot = 0, i;
  for (i = 0 ; i < ivf->k ; i++)
    tot += ivf->nbelems[i];
  return tot;
}


/* compute the imbalance factor */
double ivf_imbalance_factor (const ivf_t * ivf)
{
  int * hist = ivf->nbelems;
  int n = ivf->k, vw;
  double tot = 0, uf = 0;

  for (vw = 0 ; vw < n ; vw++) {
    tot += hist[vw];
    uf += hist[vw] * (double) hist[vw];
  }
  return uf * n / (tot * tot);
}

/* I/O */
#define WRITECHECK(a,n) if(fwrite(a,sizeof(*a),n,f)!=n) {perror("ivf_fwrite"); fprintf (stderr, "ivfc - LINE %d\n", __LINE__); return 0; }

#define READCHECK(a,n) if(fread(a,sizeof(*a),n,f)!=n) {perror("ivf_read"); fprintf (stderr, "ivfc - LINE %d\n", __LINE__); return NULL; }

int ivf_save (const char * fname, const ivf_t * ivf)
{
  int i;
  FILE * f = fopen (fname, "w");
  if (!f) {
    perror ("ivf_save - can't open the file");
    return 0;
  }

  WRITECHECK (&ivf->checksum, 1);
  WRITECHECK (&ivf->k, 1);
  WRITECHECK (&ivf->n, 1);
  WRITECHECK (&ivf->elem_size, 1);
  WRITECHECK (ivf->nbelems, ivf->k);
  
  for (i = 0 ; i < ivf->k ; i++) 
    WRITECHECK (ivf->ids[i], ivf->nbelems[i]);

  /* optionally write the complementary information */
  if (ivf->elem_size > 0)
    for (i = 0 ; i < ivf->k ; i++) 
      WRITECHECK (ivf->adat[i], ivf->elem_size * ivf->nbelems[i]);

  fclose (f);
  return 1;
}


ivf_t * ivf_load (const char * fname)
{
  int i;
  FILE * f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "Unable to open the file %s\n", fname);
    return NULL;
  }

  ivf_t *ivf = (ivf_t *) malloc (sizeof (ivf_t));
  if (!ivf) 
    return NULL;

  READCHECK (&ivf->checksum, 1);
  if (ivf->checksum != IVFCHECKSUM) {
    fprintf (stderr, "# ivf_fread: incorrect checksum\n");
    return NULL;
  }

  READCHECK (&(ivf->k), 1);
  READCHECK (&(ivf->n), 1);
  READCHECK (&(ivf->elem_size), 1);

  ivf->nbelems = (int *) malloc (sizeof (int) * ivf->k);
  ivf->seg_size = (int *) malloc (sizeof (int) * ivf->k);
  ivf->ids = (int **) malloc (sizeof (*ivf->ids) * ivf->k);
  if (!ivf->ids || !ivf->seg_size || !ivf->ids) {
    ivf_delete (ivf);
    return NULL;
  }

  READCHECK (ivf->nbelems, ivf->k);
  for (i = 0 ; i < ivf->k ; i++)
    ivf->seg_size[i] = ivf->nbelems[i] > DEFAULT_SEG_SIZE ? ivf->nbelems[i] : DEFAULT_SEG_SIZE;
  
  for (i = 0 ; i < ivf->k ; i++) {
    ivf->ids[i] = (int *) malloc (sizeof (**ivf->ids) * ivf->seg_size[i]);
    assert (ivf->ids[i]); /* not that good as a check */
    READCHECK (ivf->ids[i], ivf->nbelems[i]);
  }

  /* optionally read the complementary information */
  if (ivf->elem_size > 0) {
    ivf->adat = (unsigned char **) malloc (sizeof (*ivf->adat) * ivf->k);

    for (i = 0 ; i < ivf->k ; i++) {
      ivf->adat[i] = (unsigned char *) malloc (ivf->elem_size * ivf->seg_size[i]);
      assert (ivf->adat[i]);
      if (ivf->nbelems[i] > 0) 
	READCHECK (ivf->adat[i], ivf->elem_size * ivf->nbelems[i]);
    }
  }

  fclose (f);
  return ivf;
}
 
#undef WRITECHECK
#undef READCHECK

ivfmatch_t * ivfmatch_new (int n)
{
  return (ivfmatch_t *) malloc (n * sizeof (ivfmatch_t));
}

ivfmatch_t * ivfmatch_realloc (ivfmatch_t * m, int n)
{
  return (ivfmatch_t *) realloc (m, n * sizeof (ivfmatch_t));
}


/* Compute the set of matches according to Hamming distance threshold
   Parameters
   ids     may be anything (just used to identify the queries submitted
   keys    the quantization indexes
   adat    the set of binary signature associated with the query
   n       the number of input query vector
   
*/
ivfmatch_t * ivf_hequery (const ivf_t * ivf, 
                          const int * qids, const int * keys, 
                          const unsigned char * adat, const int nq,
                          int * buffer_size, int ht)
{						
  int i, j, posm = 0;
  int bufsize = *buffer_size;
  int elem_size = ivf->elem_size;

  ivfmatch_t * matches = ivfmatch_new (bufsize);
  ivfmatch_t * m = matches; /* For arithmetic pointer optimization */

  const unsigned char * qbs = adat;
  unsigned int h;

  for (i = 0 ; i < nq ; i++) {
    const int qid = qids[i];
    int listno = keys[i];
    int listlen = ivf_get_nb_elems (ivf, listno);
    int * listids = ivf->ids[listno];
    const unsigned char * dbs = ivf->adat[listno];

    for (j = 0 ; j < listlen ; j++) {
      /* Here perform the real work of computing the distance */
      h = hamming (qbs, dbs, elem_size);

      /* collect the match only if this satisfies the threshold */ 
      if (h <= ht) {

        /* Enough space to store another match ? */
        if (posm >= bufsize) {
          /*	  fprintf (stderr, "Realloc match buffer: %d -> %d\n", bufsize, IVF_REALLOC_NEWSIZE (bufsize)); */
	        bufsize = IVF_REALLOC_NEWSIZE (bufsize);
	        matches = ivfmatch_realloc (matches, bufsize);
	        assert (matches != NULL);
	        m = matches + posm;
	      }

        m->qid = qid;
        m->bid = listids[j];
        m->score = h; 
        m++;
        posm++;
      }
      dbs += elem_size;  /* next signature in inverted list */
    }
    qbs += elem_size;  /* next binary signature */
  }

  /* output the number of elements that been actually selected */
  *buffer_size = posm;
  return matches;
}


/* Collect matches */
hammatch_t ** ivf_he_collect (const ivf_t * ivf, const int * keys,
                              const unsigned char * qbs, int nq,
                              int ht, size_t * nmatches)
{
  int i, nbufinit = 512;
    
  /* Match entities and number of matches per query */
  hammatch_t ** hmlist = (hammatch_t **) malloc (sizeof(*hmlist) * nq);
  
#ifdef _OPENMP
#pragma omp parallel for private (i)
  for (i = 0 ; i < nq ; i++) {
    match_hamming_thres (qbs + i * ivf->elem_size, ivf->adat[keys[i]], 
                         1, ivf_get_nb_elems (ivf, keys[i]),  /* size of the inverted list */
                         ht, ivf->elem_size, nbufinit, hmlist+i, nmatches+i);
  }
#else
  for (i = 0 ; i < nq ; i++) {
    match_hamming_thres (qbs + i * ivf->elem_size, ivf->adat[keys[i]], 
                         1, ivf_get_nb_elems (ivf, keys[i]),  /* size of the inverted list */
                         ht, ivf->elem_size, nbufinit, hmlist+i, nmatches+i);
  }
#endif
  
  return hmlist;
}
  


ivfmatch_t * ivf_hequeryw (const ivf_t * ivf, 
                            const int * qids, const int * keys,
                            const unsigned char * qbs, int nq,
                            int ht, size_t * totmatches,
                            const float * score_map_, const float * list_w_)
{
  size_t i, j;  
  
  /* Match entities to count number of matches per query */
  size_t * nmatches = (size_t *) malloc (sizeof(*nmatches) * nq);
  hammatch_t ** hmlist = ivf_he_collect (ivf, keys, qbs, nq, ht, nmatches);
                                            
  /* compute the cumulative number of matches */
  size_t * cumnmatches = (size_t *) malloc (sizeof (*cumnmatches) * (nq+1));
  cumnmatches[0] = 0;
  for (i = 0 ; i < nq ; i++)
    cumnmatches[i+1] = nmatches[i] + cumnmatches[i];
  *totmatches = cumnmatches[nq];
  
  /* Populate the output structure */
  ivfmatch_t * matches = ivfmatch_new (*totmatches);
  
  /* if score_map is undefined, just returns the Hamming distances
     listweight should be NULL in this case to avoid an unexpected behavior */
  float * score_map = (float *) score_map_;
  if (score_map == NULL) {
    score_map = (float *) malloc ((ht+1) * sizeof (*score_map));
    for (i = 0 ; i <= ht ; i++) 
        score_map[i] = i;
  }

  /* list_w is typically  usedfor idf. Set to 1 by default. */
  float * list_w = (float *) list_w_;
  if (list_w == NULL) {
    list_w = (float *) malloc (ivf->k * sizeof (*list_w));
    for (i = 0 ; i < ivf->k ; i++)
      list_w[i] = 1;
  }

  for (i = 0 ; i < nq ; i++) {
    const int * listids = ivf->ids[keys[i]];
    
    ivfmatch_t * m = matches + cumnmatches[i];
    hammatch_t * hm = hmlist[i];
    
    for (j = 0 ; j < nmatches[i] ; j++) {
      m->qid = qids[i];
      m->bid = listids[hm->bid];
      m->score = score_map[hm->score] * list_w[keys[i]];
      m++;
      hm++;
    }
  }
  
  for (i = 0 ; i < nq ; i++)
    free (hmlist[i]);
  free (hmlist);
  
  free (cumnmatches);
  free (nmatches);
  
  if (score_map_ == NULL)
    free (score_map);
  if (list_w == NULL)
    free (list_w);
  return matches;
}


/* Collect cross-matches with Hamming distance */
hammatch_t ** ivf_he_collect_crossmatches (const ivf_t * ivf, int ht, size_t * nmatches)
{
  int i, nbufinit = 512;
  
  /* Match entities and number of matches per query */
  hammatch_t ** hmlist = (hammatch_t **) malloc (sizeof(*hmlist) * ivf->k);
  
#ifdef _OPENMP
#pragma omp parallel for private (i)
  for (i = 0 ; i < ivf->k ; i++) {
    crossmatch_hamming (ivf->adat[i], ivf_get_nb_elems (ivf, i), 
                        ht, ivf->elem_size, nbufinit, hmlist+i, nmatches+i);
    
    hammatch_t *m = hmlist[i];
    const int * listids = ivf->ids[i];
    long j, n = nmatches[i];
    
    for (j = 0 ; j < n ; j++) {
      m->qid = listids[m->qid];
      m->bid = listids[m->bid];
      m++;
    }
    
  }
#else
  for (i = 0 ; i < ivf->k ; i++) {
    crossmatch_hamming (ivf->adat[i], ivf_get_nb_elems (ivf, i), 
                        ht, ivf->elem_size, nbufinit, hmlist+i, nmatches+i);
    hammatch_t *m = hmlist[i];
    const int * listids = ivf->ids[i];
    int j;
    for (j = 0 ; j < nmatches[i] ; j++) {
        m->qid = listids[m->qid];
        m->bid = listids[m->bid];
        m++;
    }
  }
#endif

  return hmlist;
}


/* Collect cross-matches with Hamming distance */
void ivf_he_crossmatches_prealloc (const ivf_t * ivf, int ht, 
                                   int * idx, uint16 * hams, 
                                   size_t * cumnmatches)
{
  long i;
#ifdef _OPENMP
#pragma omp parallel for private (i)
  for (i = 0 ; i < ivf->k ; i++) {
    crossmatch_hamming_prealloc (ivf->adat[i], ivf_get_nb_elems (ivf, i), 
                                 ht,  ivf->elem_size,
                                 idx + 2 * cumnmatches[i], 
                                 hams + cumnmatches[i]);
    
    long n = cumnmatches[i+1] - cumnmatches[i];
    int * m = idx + 2 * cumnmatches[i];
    
    const int * listids = ivf->ids[i];
    
    long j;
    for (j = 0 ; j < n ; j++) {
      *m = listids[*m]; m++;
      *m = listids[*m]; m++;
    }    
  }
#else
  for (i = 0 ; i < ivf->k ; i++) {
    int nout = crossmatch_hamming_prealloc (ivf->adat[i], ivf_get_nb_elems (ivf, i), 
                                            ht, ivf->elem_size,
                                            idx + 2 * cumnmatches[i], 
                                            hams + cumnmatches[i]);
    
    long n = cumnmatches[i+1] - cumnmatches[i];
    assert (nout == n);
    
    int * m = idx + 2 * cumnmatches[i];
    
    const int * listids = ivf->ids[i];
    
    long j;
    for (j = 0 ; j < n ; j++) {
      *m = listids[*m]; m++;
      *m = listids[*m]; m++;
    }    
  }
#endif
}


/* Couny cross-matches with Hamming distance */
void ivf_he_count_crossmatches (const ivf_t * ivf, int ht, size_t * nmatches)
{
  long i;
  
#ifdef _OPENMP
#pragma omp parallel for private (i)
  for (i = 0 ; i < ivf->k ; i++) {
    crossmatch_hamming_count (ivf->adat[i], ivf_get_nb_elems (ivf, i),  
                              ht, ivf->elem_size, nmatches+i);   
  }
#else
  for (i = 0 ; i < ivf->k ; i++) {
    crossmatch_hamming_count (ivf->adat[i], ivf_get_nb_elems (ivf, i), 
                              ht, ivf->elem_size, nmatches+i);   
  }
#endif
}








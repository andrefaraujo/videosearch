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


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "binheap.h"
#include "sorting.h"


size_t fbinheap_sizeof (int maxk) 
{
  return sizeof (fbinheap_t) + maxk * (sizeof (float) + sizeof (int));
}


void fbinheap_init (fbinheap_t *bh, int maxk) 
{
  bh->k = 0;
  bh->maxk = maxk;

  /* weirdly index the nodes from 1 (for the root) */
  bh->val = (float *) ((char *) bh + sizeof (fbinheap_t)) - 1;
  bh->label = (int *) ((char *) bh + sizeof (fbinheap_t)
		       + maxk * sizeof (*bh->val)) - 1;
} 


fbinheap_t * fbinheap_new (int maxk)
{
  int bhsize = fbinheap_sizeof (maxk);
  fbinheap_t * bh = (fbinheap_t *) malloc (bhsize);
  fbinheap_init (bh,maxk);
  return bh;
}

void fbinheap_reset (fbinheap_t *bh)
{
  bh->k = 0;
}

void fbinheap_delete (fbinheap_t * bh)
{
  free (bh);
}


void fbinheap_pop (fbinheap_t * bh)
{
  assert (bh->k > 0);

  float val = bh->val[bh->k];
  int i = 1, i1, i2;

  while (1) {
    i1 = i << 1;
    i2 = i1 + 1;

    if (i1 > bh->k)
      break;

    if (i2 == bh->k + 1 || bh->val[i1] > bh->val[i2]) {
      if (val > bh->val[i1])
        break;
      bh->val[i] = bh->val[i1];
      bh->label[i] = bh->label[i1];
      i = i1;
    } 
    else {
      if (val > bh->val[i2])
        break;
      
      bh->val[i] = bh->val[i2];
      bh->label[i] = bh->label[i2];
      i = i2;
    }
  }
  
  bh->val[i] = bh->val[bh->k];
  bh->label[i] = bh->label[bh->k];
  bh->k--;
}


static void fbinheap_push (fbinheap_t * bh, int label, float val)
{
  assert (bh->k < bh->maxk);

  int i = ++bh->k, i_father;

  while (i > 1) {
    i_father = i >> 1;
    if (bh->val[i_father] >= val)  /* the heap structure is ok */
      break; 

    bh->val[i] = bh->val[i_father];
    bh->label[i] = bh->label[i_father];

    i = i_father;
  }
  bh->val[i] = val;
  bh->label[i] = label;
}


void fbinheap_add (fbinheap_t * bh, int label, float val)
{
  if (bh->k < bh->maxk) {
    fbinheap_push (bh, label, val);
    return;
  }

  if (val < bh->val[1]) {  
    fbinheap_pop (bh);
    fbinheap_push (bh, label, val);
  }
}


void fbinheap_addn (fbinheap_t * bh, int n, const int * label, const float * v)
{
  int i;

  for (i = 0 ; i < n && bh->k < bh->maxk; i++)
    if(!isnan(v[i]))      
      fbinheap_push (bh, label[i], v[i]);

  float lim=bh->val[1];
  for ( ; i < n; i++) {
    if(v[i]<lim) {
      fbinheap_pop (bh);
      fbinheap_push (bh, label[i], v[i]);
      lim=bh->val[1];
    }      
  }
}


void fbinheap_addn_label_range (fbinheap_t * bh, int n, int label0, const float * v)
{
  int i;

  for (i = 0 ; i < n && bh->k < bh->maxk; i++)
    if(!isnan(v[i]))
      fbinheap_push (bh, label0+i, v[i]);

  float lim=bh->val[1];
  for ( ; i < n; i++) { /* optimized loop for common case */
    if(v[i]<lim) {
      fbinheap_pop (bh);
      fbinheap_push (bh, label0+i, v[i]);
      lim=bh->val[1];
    }      
  }

}


static int cmp_floats (const void * v1, const void * v2)
{
  if (*(float *) v2 == *(float *) v1)
    return 0;
  
  return (*(float *) v1 > *(float *) v2) ? 1 : -1;
}


void fbinheap_sort_labels (fbinheap_t * bh, int * perm)
{
  int i;
  fvec_sort_index (bh->val + 1, bh->k, perm);
  for (i = 0 ; i < bh->k ; i++)
    perm[i] = bh->label[perm[i]+1];
}


void fbinheap_sort_values (fbinheap_t * bh, float * v)
{
  memcpy (v, bh->val + 1, bh->k * sizeof (*v));
  qsort (v, bh->k, sizeof (*v), &cmp_floats);
}


void fbinheap_sort_per_labels (fbinheap_t * bh, int * perm, float *v) {
  ivec_sort_index(bh->label+1, bh->k, perm);
  int i;
  if(v) {
    for(i=0;i<bh->k;i++)  {
      int heappos = perm[i] + 1;
      perm[i] = bh->label[heappos];
      v[i] = bh->val[heappos];
    }
  } else {
    for(i=0;i<bh->k;i++)  
      perm[i] = bh->label[perm[i] + 1];    
  }     
}



void fbinheap_sort (fbinheap_t * bh, int * perm, float *v)
{
  int i, heappos;
  /* TODO use binheap structure to get in the correct order */
  fvec_sort_index (bh->val + 1, bh->k, perm);
  for (i = 0 ; i < bh->k ; i++) {
    heappos = perm[i] + 1;
    perm[i] = bh->label[heappos];
    v[i] = bh->val[heappos];
  }
}


void fbinheap_display (fbinheap_t * bh)
{
  int i;
  printf ("[nel = %d / maxel = %d] ", bh->k, bh->maxk);

  for (i = 1 ; i <= bh->k ; i++)
    printf ("%d %.6f / ", bh->label[i], bh->val[i]);
  printf ("\n");
}

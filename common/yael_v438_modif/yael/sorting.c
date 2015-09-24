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

#include <assert.h>
#include <string.h>
#include <math.h>

#include "sorting.h"
#include "machinedeps.h"
#include "binheap.h"
#include "vector.h"

#define NEWA(type,n) (type*)malloc(sizeof(type)*(n))


static int compare_for_k_min (const void *v1, const void *v2)
{
  return *(*(float **) v1) > *(*(float **) v2) ? 1 : -1;
}



/*--------------------------------------------------------------------------
  The following function are related to the Hoare selection algorithm (also know as quickselect). 
  This is a "lazy" version of the qsort algorithm. 
  It is used to find some quantile of the values in a table 
*/

#define PERM(i) (*perm[i])
#define SWAPFPTR(i,j) {const float* tmp = perm[i]; perm[i] = perm[j]; perm[j] = tmp; }

/* order perm[i0..i1-1] such that *perm[i] <= *perm[j]
   forall i0<=i<q and q<=j<i1  */
static void hoare_selectp (const float **perm, int i0, int i1, int q)
{
  float pivot = PERM(i0);
  int j0, j1, lim;
  assert (i1 - i0 > 1 && q > i0 && q < i1);

  for (j0 = i0, j1 = i1 ; 1 ; ) {
    while (j0 < j1 - 1) {
      j0++;
      if (PERM(j0) > pivot)
        goto endseginf;
    }
    lim = j1;
    break;
  endseginf:
    while (j1 - 1 > j0) {
      j1--;
      if (PERM(j1) <= pivot)
        goto endsegsup;
    }
    lim = j0;
    break;
  endsegsup:
    SWAPFPTR (j0, j1);
  }
  assert (lim > i0);
  if (lim == i1) {
    SWAPFPTR (i0, i1 - 1);
    lim = i1 - 1;
  }
  if (lim == q)
    return; 
  else if (q < lim)
    hoare_selectp (perm, i0, lim, q);
  else
    hoare_selectp (perm, lim, i1, q);
}

#undef PERM
#undef SWAPFPTR


/*  The same Hoare algorithm, but which that modifies its input */
#define PERM(i) f[i]
#define SWAPFLOAT(i,j) {float tmp = f[i]; f[i]=f[j]; f[j]=tmp; }

static void hoare_select_f (float *f, int i0, int i1, int q)
{
  float pivot = PERM(i0);
  int j0, j1, lim;
  assert (i1 - i0 > 1 && q > i0 && q < i1);

  for (j0 = i0, j1 = i1 ; 1 ; ) {
    while (j0 < j1 - 1) {
      j0++;
      if (PERM(j0) > pivot)
        goto endseginf;
    }
    lim = j1;
    break;
  endseginf:
    while (j1 - 1 > j0) {
      j1--;
      if (PERM(j1) <= pivot)
        goto endsegsup;
    }
    lim = j0;
    break;
  endsegsup:
    SWAPFLOAT (j0, j1);
  }
  assert (lim > i0);
  if (lim == i1) {
    SWAPFLOAT (i0, i1 - 1);
    lim = i1 - 1;
  }
  if (lim == q)
    return;                     /* mission accomplished */
  if (q < lim)
    hoare_select_f (f, i0, lim, q);
  else
    hoare_select_f (f, lim, i1, q);
}

#undef PERM
#undef SWAPFLOAT




/*--- Find k biggest elements of array val[n] and order them ---*/

/* version based on the Hoare algorithm (also called qselect)*/
static void fvec_k_max_hoare (const float *val, int n, int *idx, int k)
{

  const float **idx_ptr = NEWA (const float *, n);

  int i;
  for (i = 0; i < n; i++)
    idx_ptr[i] = val + i;

  if (k < n)
    hoare_selectp (idx_ptr, 0, n, n - k);

  /* sort upper part of array */
  qsort (idx_ptr + n - k, k, sizeof (*idx_ptr), compare_for_k_min);

  for (i = 0; i < k; i++)
    idx[i] = idx_ptr[n - 1 - i] - val;  /* Pointer arithmetic */

  free (idx_ptr);
}


/* maxheap version */
static void fvec_k_max_maxheap (const float *val, int n,
                                       int *idx, int k)
{
  fbinheap_t *mh = fbinheap_new (k);
  int i;

  for (i = 0; i < n; i++)
    fbinheap_add (mh, i, -val[i]);       /* -val because we want maxes instead of mins */

  fbinheap_sort_labels (mh,idx);
  fbinheap_delete (mh);
}


void fvec_k_max (const float *val, int n, int *idx, int k)
{
  assert (k <= n);

  if (n == 0 || k == 0)
    return;

  /* TODO: find out where the limit really is */
  if (n > 4 * k)
    fvec_k_max_maxheap (val, n, idx, k);
  else
    fvec_k_max_hoare (val, n, idx, k);

  /* The 2 algorithms are not strictly equivalent because they order
     same-values differently:

     - hoare orders them arbitrarily (but it is reproducible). TODO:
     compare values and pointers to make it deterministic

     - maxheap orders lexicographically by (val[i],i).
   */
}




/*--- Idem for smallest ---*/

/* Hoare version */
static void fvec_k_min_hoare (const float *val, int n, int *idx, int k)
{

  const float **idx_ptr = NEWA (const float *, n);

  int i;
  for (i = 0; i < n; i++)
    idx_ptr[i] = val + i;

  if (k < n)
    hoare_selectp (idx_ptr, 0, n, k);

  /* sort lower part of array */
  qsort (idx_ptr, k, sizeof (*idx_ptr), compare_for_k_min);

  for (i = 0; i < k; i++)
    idx[i] = idx_ptr[i] - val; 

  free (idx_ptr);
}


/* maxheap version */
static void fvec_k_min_maxheap (const float *val, int n,
                                     int *idx, int k)
{
  fbinheap_t *mh = fbinheap_new (k);
  int i;

  for (i = 0; i < n; i++)
    fbinheap_add (mh, i, val[i]);    /* -val because we want maxes instead of mins */

  fbinheap_sort_labels (mh, idx);
  fbinheap_delete (mh);
}


void fvec_k_min (const float *val, int n, int *idx, int k)
{
  assert (k <= n);

  if (n == 0 || k == 0)
    return;

  if (k == 1) {
    *idx = fvec_arg_min (val, n);
    return; 
  }

  /* TODO: find out where the limit really is */
  if (n > 3 * k)
    fvec_k_min_maxheap (val, n, idx, k); 
  else
    fvec_k_min_hoare (val, n, idx, k);
}



/*********************************************************************
 * Simple functions 
 *********************************************************************/

void find_labels (const int *labels, int nres, int *ilabels, int nilabels)
{
  int aux[nilabels];
  int i, j, left = nilabels;
  memcpy (aux, ilabels, sizeof (int) * nilabels);
  memset (ilabels, 0xff, sizeof (int) * nilabels);
  for (i = 0; i < nres; i++) {
    for (j = 0; j < nilabels; j++)
      if (labels[i] == aux[j]) {
        if (ilabels[j] < 0) {
          ilabels[j] = i;
          left--;
          if (!left)
            return;
        }
      }
  }
}

#ifdef HAVE_TLS

static __thread const float * tab_to_sort_f;

static int compare_for_sort_index_f (const void *v1, const void *v2)
{
#elif defined(HAVE_QSORT_R)
static int compare_for_sort_index_f (void *thunk, const void *v1, const void *v2)
{
  const float *tab_to_sort_f=thunk;
#else 
#error "please provide some kind of thread-local storage"
#endif
  
  float dt = tab_to_sort_f[*(int *)v1] - tab_to_sort_f[*(int *)v2];
  if (dt) 
    return dt>0 ? 1 : -1;
  return *(int *)v1 - *(int *)v2;
}



void fvec_sort_index(const float *tab,int n,int *perm) {
  int i;

  for (i = 0 ; i < n ; i++) 
    perm[i] = i;

#ifdef HAVE_TLS
  tab_to_sort_f = tab;
  qsort (perm, n, sizeof(int), compare_for_sort_index_f);
#elif defined(HAVE_QSORT_R)
  qsort_r (perm, n, sizeof(int), (void*)tab, compare_for_sort_index_f);
#endif
}



#ifdef HAVE_TLS

static __thread const int * tab_to_sort;

static int compare_for_sort_index (const void *v1, const void *v2)
{
#elif defined(HAVE_QSORT_R)
static int compare_for_sort_index (void *thunk, const void *v1, const void *v2)
{
  const int *tab_to_sort = thunk;
#else 
#error "please provide some kind of thread-local storage"
#endif
  
  int dt = tab_to_sort[*(int *)v1] - tab_to_sort[*(int *)v2];
  if (dt) 
    return dt;
  return *(int *)v1 - *(int *)v2;
}


void ivec_sort_index (const int *tab, int n, int *perm) 
{
  int i;

  for (i = 0 ; i < n ; i++) 
    perm[i] = i;

#ifdef HAVE_TLS
  tab_to_sort = tab;
  qsort (perm, n, sizeof(int), compare_for_sort_index);
#elif defined(HAVE_QSORT_R)
  qsort_r (perm, n, sizeof(int), (void*)tab, compare_for_sort_index);
#endif
}


void ivec_invert_perm (const int *perm, int n, int *iperm) 
{
  int i;
  for (i = 0 ; i < n ; i++) 
    iperm[perm[i]] = i;
}


static int compare_for_ivec_sort (const void *v1, const void *v2)
{
  int dt = *(int *)v1 - *(int *)v2;
  if (dt) 
    return dt;
  return v1 - v2;
}

void ivec_sort(int *tab, int n) 
{
  qsort (tab, n, sizeof(int), compare_for_ivec_sort);
}
 
 static int compare_for_fvec_sort (const void *v1, const void *v2)
{
  float dt = *(float *)v1 - *(float *)v2;
  return dt>0 ? 1 : dt<0 ? -1 : v1 - v2;
}


void fvec_sort(float *tab, int n) 
{
  qsort (tab, n, sizeof(int), compare_for_fvec_sort);
}


void fvecs_sort (float * v, int d, int n)
{
  int i;
  for (i = 0 ; i < n ; i++)
    fvec_sort (v + i * d, d);
}


/* sort according to the input permutation */
void ivec_sort_by_permutation (int * v, const int * order, int n)
{
  int i;
  int * o = malloc (n * sizeof (*o));

  for (i = 0 ; i < n ; i++)
    o[order[i]] = i;

  for (i = 0 ; i < n ; i++)
    while (o[i] != i) {
      int newpos = o[i];

      int a = v[i];
      v[i] = v[newpos];
      v[newpos] = a;

      int b = o[newpos];
      o[newpos] = o[i];
      o[i] = b;
    }

  free (o);
}


void fvec_sort_by_permutation (float * v, const int * order, int n)
{
  int i;
  float * o = malloc (n * sizeof (*o));

  for (i = 0 ; i < n ; i++)
    o[order[i]] = i;

  for (i = 0 ; i < n ; i++)
    while (o[i] != i) {
      int newpos = o[i];

      float a = v[i];
      v[i] = v[newpos];
      v[newpos] = a;

      float b = o[newpos];
      o[newpos] = o[i];
      o[i] = b;
    }

  free (o);
} 


/*********************************************************************
 * Get ranks of a few values
 *********************************************************************/

/* Histogram with irregular bins */
typedef struct {
  float val;
  int no;                       /* index in the vals table */
  int n_gt;                     /* nb of values greater than val */
  int n_eq;                     /* nb of values equal to val */

  int copy_from_no;             /* handle duplicate values */
} stop_t;

static int cmp_stops_for_increasing (const void *v1, const void *v2)
{
  float f1=((stop_t *) v1)->val;
  float f2=((stop_t *) v2)->val;
  /* NaN considered smallest */
  return isnan(f2) || f2 < f1 ? 1 : -1;
}


static stop_t *make_stops (const float *tab, int n,
                           const float *vals, int nval, int *nstop_out)
{
  int i;
  stop_t *stops = NEWA (stop_t, nval + 1);

  for (i = 0; i < nval; i++) {
    stop_t *stop = stops + i;
    stop->val = vals[i];
    stop->no = i;
    stop->n_gt = stop->n_eq = 0;
  }
  {                             /* lower sentinel to count values below the lowest val and NaNs */
    stop_t *stop = stops + nval;
    stop->val = -1e30;          /* -1.0/0.0 does not work ?? */
    stop->no = -1;
    stop->n_gt = stop->n_eq = 0;
  }
  qsort (stops, nval + 1, sizeof (stop_t), &cmp_stops_for_increasing);
  for(i=1;i<nval+1;i++) assert(stops[i-1].val<=stops[i].val);
  assert (stops[0].no == -1 || !isfinite(stops[0].val));

  /* handle duplicate vals */
  int nstop = nval + 1;
  i = 1;
  while (i < nstop) {
    if (stops[i].val == stops[i - 1].val) {
      int copied_no = stops[i].no;
      /* shift array */
      memmove (stops + i, stops + i + 1, (nstop - i - 1) * sizeof (stop_t));
      nstop--;
      stops[nstop].no = copied_no;
      stops[nstop].copy_from_no = stops[i - 1].no;
    } else i++;
  }

  for (i = 0; i < n; i++) {
    float v = tab[i];
    /* bisection, property: stops[t0].val <= v < stops[t1].val
     * NaNs will be counted in stops[0].n_gt */
    int t0 = 0, t1 = nstop;
    while (t0 + 1 != t1) {
      int med = (t0 + t1) / 2;
      if (stops[med].val <= v)
        t0 = med;
      else
        t1 = med;
    }
    if (stops[t0].val < v)
      stops[t0].n_gt++;
    else  /* NaN's and other weird values here... */
      stops[t0].n_eq++;
  }

  *nstop_out = nstop;

  return stops;
}

void fvec_ranks_of (const float *tab, int n,
                      const float *vals, int nval,
                      int *minranks, int *maxranks)
{

  int i, nstop;

  stop_t *stops = make_stops (tab, n, vals, nval, &nstop);

  /* sum up ranks */
  int rank = 0;
  for (i = nstop - 1; i > 0; i--) {
    stop_t *stop = stops + i;
    int no = stop->no;
    rank += stop->n_gt;
    if (minranks)
      minranks[no] = rank;
    rank += stop->n_eq;
    if (maxranks)
      maxranks[no] = rank;
  }

  /* handle duplicates */
  for (i = nstop; i <= nval; i++) {
    int src = stops[i].copy_from_no;
    int dst = stops[i].no;
    if (minranks)
      minranks[dst] = minranks[src];
    if (maxranks)
      maxranks[dst] = maxranks[src];
  }

  free (stops);
}


void fvec_ranks_inc_of (const float *tab, int n,
                          const float *vals, int nval,
                          int *minranks, int *maxranks)
{
  int i, nstop;

  stop_t *stops = make_stops (tab, n, vals, nval, &nstop);

  /* sum up ranks. We don't add stops[0].n_eq because it contains NaNs
     (to be ignored) */
  int rank = stops[0].n_gt;

  for (i = 1; i < nstop; i++) {
    stop_t *stop = stops + i;
    int no = stop->no;
    if (minranks)
      minranks[no] = rank;
    rank += stop->n_eq;
    if (maxranks)
      maxranks[no] = rank;
    rank += stop->n_gt;
  }

  /* handle duplicates */
  for (i = nstop; i <= nval; i++) {
    int src = stops[i].copy_from_no;
    int dst = stops[i].no;
    if (minranks)
      minranks[dst] = minranks[src];
    if (maxranks)
      maxranks[dst] = maxranks[src];
  }

  free (stops);
}

int ivec_sorted_find (const int *v,int n,int val) 
{
  if(n==0 || v[0]>val) 
    return -1;
  
  int i0 = 0;
  int i1 = n;

  /* find an occurrence */
  while (i1 > i0 + 1) {
    int i2 = (i0 + i1) / 2;
    if (val < v[i2]) 
      i1 = i2;
    else
      i0 = i2;
  }
  return i0;
}



/*! count occurrences of val in sorted vector */
int ivec_sorted_count_occurrences (const int *v, int n, int val) 
{
  int i0 = ivec_sorted_find (v, n, val);
  
  if (i0<0 || v[i0]!=val) 
    return 0;
  
  /* explore surroundings */
  int i1 = i0 + 1; 
  while (i1 < n && v[i1] == val) 
    i1++;
  i0--;
  while (i0 >= 0 && v[i0] == val) 
    i0--;
  
  return i1 - i0 - 1;
}


int ivec_sorted_count_occurrences_multiple(const int *v,int n,
                                           const int *vals,int nval) 
{
  int i;
  int accu=0;

  for (i = 0 ; i < nval ; i++) 
    accu += ivec_sorted_count_occurrences (v, n, vals[i]);

  return accu;
}


/*! count unique occurrences  */
int ivec_sorted_count_unique (const int *v, int n) 
{
  int i = 0;
  int count = 0;
  
  while (i < n) {
    count++;
    int v0 = v[i];
    while (i < n && v[i] == v0) 
      i++;
  }
  
  return count;
}



float fvec_median (float *f, int n)
{
  if(n == 0) 
    return 0.0 / 0.0; 

  if (n == 1) 
    return f[0];

  int halfn = n / 2;
  int j;

  hoare_select_f (f, 0, n, halfn);

  float min_upper = f[halfn];
  for (j = halfn + 1; j < n; j++)
    if (f[j] < min_upper)
      min_upper = f[j];

  if (n % 2 == 1)
    return min_upper;
  else {
    float max_lower = f[0];
    for (j = 1; j < halfn; j++)
      if (f[j] > max_lower)
        max_lower = f[j];
    return 0.5 * (min_upper + max_lower);
  }
}


float fvec_median_const (const float *f, int n) 
{
  float *f2 = fvec_new (n); 
  fvec_cpy (f2, f, n);
  float med = fvec_median (f2, n);
  free (f2);
  return med;
}


float fvec_min (const float *f, long n) 
{
  assert (n > 0);
  float m = f[0];
  long i;
  for (i = 1 ; i < n ; i++) 
    if (f[i] < m) 
      m = f[i];
  return m;
}


int ivec_min (const int *f, long n) 
{
  assert (n > 0);
  int m = f[0];
  long i;
  for (i = 1 ; i < n ; i++) 
    if (f[i] < m) 
      m = f[i];
  return m;
}


float fvec_max (const float *f, long n) 
{
  assert (n > 0);
  float m = f[0];
  long i;
  for (i = 1 ; i < n ; i++) 
    if (f[i] > m) 
      m = f[i];
  return m;
}


int ivec_max (const int *f, long n) 
{
  assert (n > 0);
  int m = f[0];
  long i;
  for (i = 1 ; i < n ; i++) 
    if (f[i] > m) 
      m = f[i];
  return m;
}


int fvec_arg_max (const float *f, long n) 
{
  assert (n > 0);
  float m = f[0];
  long i,i0 = 0;
  for (i = 1 ; i < n ; i++) 
    if (f[i] > m) {
      m = f[i]; 
      i0 = i; 
    }
  return i0;
}


int fvec_arg_min (const float *f, long n) 
{
  assert (n > 0);
  float m = f[0];
  long i, i0 = 0;
  for (i = 1 ; i < n ; i++) 
    if(f[i] < m) {
      m = f[i]; 
      i0 = i; 
    }
  return i0;
}


float fvec_quantile (float *f, int n, int q)
{
  assert (n>0);
  if (n==1) return f[0];
  if (q==0) return fvec_min(f, n); /* any value */
  if (q>=n) return fvec_max(f, n);
  
  hoare_select_f (f, 0, n, q);

  float min_upper = f[q];
  int j;
  for (j = q + 1; j < n; j++)
    if (f[j] < min_upper)
      min_upper = f[j];
  
  return min_upper;
}

/*********************************************************************
 * Ordered set merge
 *********************************************************************/


int merge_ordered_sets (const int **labels,const float **vals,
			const int *sizes,int k,
			int **labels_out,float **vals_out) {
  int i,j;
  int n_out = ivec_sum (sizes, k);

  int *all_labels = ivec_new (n_out);
  float *all_vals = fvec_new (n_out);

  /* Maxheap:
   * * maxheap label = index of table in 0..k-1
   * * maxheap val = - (label from labels table)
   *  
   * If maxheap val does not fit in a float (if label>2**24), it
   * triggers an assertion. Time to implement a maxheap with int
   * values...
   */
  fbinheap_t *mh = fbinheap_new(k);

  /* current index on table k */ 
  int indices[k];

  for ( i = 0 ; i < k ; i++) {
    if (sizes[i] == 0) 
      continue;
    indices[i] = 0;
    int label = labels[i][0];
    float mh_val = -label;
    assert ((int)(-mh_val) == label || !"lost precision in int->float conversion");
    fbinheap_add (mh, i, mh_val);
  }
  
  int all_i = 0;
  while (mh->k>0) {    

    /* smallest available label */    
    i = mh->label[1];       /* index of table */
    j = (int)(-mh->val[1]); /* label */

    /* I don't dare compiling with -DNDEBUG */    
    /* assert(j==labels[i][indices[i]]); */

    all_labels[all_i] = j;
    all_vals[all_i] = vals[i][indices[i]];
    all_i++;

    /* remove handled label */
    fbinheap_pop (mh);
    
    indices[i]++;
    if (indices[i] < sizes[i]) { /* push next label from this table */
      int label = labels[i][indices[i]];
      float mh_val = -label;
      assert ((int)(-mh_val) == label || !"lost precision in int->float conversion");
      fbinheap_add (mh, i, mh_val);
    }
  }
  fbinheap_delete (mh);  
  assert (all_i == n_out);

  *labels_out = all_labels;
  *vals_out = all_vals;
  return n_out;
}



int compress_labels_by_disratio (int *labels, const float *vals, int n, float ratio) {
  if (!n) 
    return 0;
  int i;
  float min = vals[0];
  for (i = 1 ; i < n ; i++)
    if (vals[i] < min) 
      min = vals[i];

  float thresh = min * ratio;
  int j = 0;
  for (i = 0 ; i < n ; i++) 
    if (vals[i] < thresh)
      labels[j++] = labels[i];

  /* clear out the rest */
  for (i = j ; i < n ; i++) 
    labels[i]=-1;
        
  return j;
}

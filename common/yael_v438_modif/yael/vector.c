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

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <errno.h>

/* for mmap */
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>


#include "vector.h"


#ifdef __linux__
#include <malloc.h>
#else
static void *memalign(size_t ignored,size_t nbytes) {
  return malloc(nbytes); 
} 
#endif


/*-------------------------------------------------------------*/
/* Allocation                                                  */
/*-------------------------------------------------------------*/



float *fvec_new (long n)
{
  float *ret = (float *) memalign (16, sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "fvec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


double *dvec_new (long n)
{
  double *ret = (double *) memalign (16, sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "dvec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


 int *ivec_new (long n)
{
  int *ret = (int *) malloc (sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "ivec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}

unsigned char * bvec_new (long n)
{
  unsigned char *ret = (unsigned char *) memalign (16, sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "bvec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


long long * lvec_new (long n)
{
  long long *ret = (long long *) malloc (sizeof (*ret) * n);
  if (!ret) {
    fprintf (stderr, "ivec_new %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


float *fvec_new_0 (long n)
{
  float *ret = (float *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "fvec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


double * dvec_new_0 (long n)
{
  double *ret = (double *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "fvec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


long long *lvec_new_0 (long n)
{
  long long *ret = (long long*) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "fvec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}


float *fvec_new_nan (long n)
{
  float *ret = fvec_new(n);
  fvec_nan(ret,n);
  return ret;
}

/*-------------------------------------------------------------*/
/* Random                                                      */
/*-------------------------------------------------------------*/

/* Generate Gaussian random value, mean 0, variance 1 (from Python source) */


static double drand_r(unsigned int *seed) {
  return rand_r(seed)/((double)RAND_MAX + 1.0);
}

#define NV_MAGICCONST  1.71552776992141

static double gaussrand_r (unsigned int *seed)
{
  double z;
  while (1) {
    float u1, u2, zz;
    u1 = drand_r (seed);
    u2 = drand_r (seed);
    z = NV_MAGICCONST * (u1 - .5) / u2;
    zz = z * z / 4.0;
    if (zz < -log (u2))
      break;
  }
  return z;
}

double gaussrand () {
  unsigned int seed=lrand48();
  return gaussrand_r(&seed);
}


void fvec_rand (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = drand48();
}

void fvec_rand_r (float * v, long n, unsigned int seed)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = drand_r(&seed);
}


void fvec_randn (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = gaussrand();
}

void fvec_randn_r (float * v, long n, unsigned int seed)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = gaussrand_r(&seed);
}


float *fvec_new_rand (long n) 
{
  float * f = fvec_new (n);
  fvec_rand (f, n);
  return f;
}

float * fvec_new_randn (long n)
{
  float * f = fvec_new (n);
  fvec_randn(f,n);
  return f;
}

float * fvec_new_rand_r (long n, unsigned int seed)
{
  float * f = fvec_new (n);
  fvec_rand_r(f,n,seed);
  return f;
}

float * fvec_new_randn_r (long n, unsigned int seed)
{
  float * f = fvec_new (n);
  fvec_randn_r(f,n,seed);
  return f;
}


int * ivec_new_random_idx (int n, int k)
{
  return ivec_new_random_idx_r (n, k, lrand48());
}

int * ivec_new_random_idx_r (int n, int k, unsigned int seed)
{
  int *idx = ivec_new (n);
  int i;

  for (i = 0; i < n; i++)
    idx[i] = i;

  for (i = 0; i < k ; i++) {
    int j = i +  rand_r(&seed) % (n - i);
    /* swap i and j */
    int p = idx[i];
    idx[i] = idx[j];
    idx[j] = p;
  }

  return idx;
}

int *ivec_new_random_perm (int n)
{
  return ivec_new_random_idx (n, n - 1);
}

int *ivec_new_random_perm_r (int n, unsigned int seed)
{
  return ivec_new_random_idx_r (n, n - 1, seed);
}


void ivec_shuffle (int * v, long n)
{
  int i;

  for (i = 0; i < n - 1; i++) {
    int j = i + random () % (n - i);
    /* swap i and j */
    int p = v[i];
    v[i] = v[j];
    v[j] = p;
  }
}

/*-------------------------------------------------------------*/
/* Allocate & initialize                                       */
/*-------------------------------------------------------------*/


int *ivec_new_0 (long n)
{
  int *ret = (int *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "ivec_new_0 %ld : out of memory\n", n);
    abort();
  }
  return ret;
}

unsigned char *bvec_new_0 (long n)
{
  unsigned char *ret = bvec_new (n);
  if (!ret) {
    fprintf (stderr, "ivec_new_0 %ld : out of memory\n", n);
    abort();
  }
  memset (ret, 0, n);
  return ret;
}



float *fvec_new_set (long n, float val)
{
  int i;
  float *ret = (float *) calloc (sizeof (*ret), n);
  if (!ret) {
    fprintf (stderr, "fvec_new_set %ld : out of memory\n", n);
    abort();
  }

  for (i = 0 ; i < n ; i++)
    ret[i] = val;

  return ret;
}


int *ivec_new_set (long n, int val)
{
  int i;
  int *ret = ivec_new(n);

  for (i = 0 ; i < n ; i++)
    ret[i] = val;

  return ret;
}


int * ivec_new_range (long a, long b)
{
  int i;
  int *ret = ivec_new(b - a);

  for (i = a ; i < b ; i++)
    ret[i - a] = i;

  return ret;
}


int * ivec_new_cpy (const int * v, long n)
{
  int *ret = ivec_new(n);  
  memcpy (ret, v, n * sizeof (*ret));
  return ret;
}


float * fvec_new_cpy (const float * v, long n) {
  float *ret = fvec_new(n);
  memcpy (ret, v, n * sizeof (*ret));
  return ret;  
}


/*-------------------------------------------------------------*/
/* resize                                                      */
/*-------------------------------------------------------------*/


float * fvec_resize (float * v, long n)
{
  float * v2 = realloc (v, n * sizeof (*v));
  return v2;
}


int * ivec_resize (int * v, long n)
{
  int * v2 = realloc (v, n * sizeof (*v));
  return v2;
}



/*-------------------------------------------------------------*/
/* statistics                                                  */
/*-------------------------------------------------------------*/


int *ivec_new_histogram (int k, const int *v, long n)
{
  long i;
  int *h = ivec_new_0 (k);

  for (i = 0; i < n; i++) {
    assert (v[i] >= 0 && v[i] < k);
    h[v[i]]++;
  }

  return h;
}


int * ivec_new_histogram_clip (int k, int * v, long n) {
  long i;
  int *h = ivec_new_0 (k);

  for (i = 0; i < n; i++) {
    int val=v[i];
    if(val>=k) val=k-1;
    if(val<0) val=0;
    h[val]++;
  }

  return h;
}

int * fvec_new_histogram_clip (float vmin,float vmax, int k, float *v, long n) {
  long i;
  int *h = ivec_new_0 (k);
  float istep=k/(vmax-vmin);

  for (i = 0; i < n; i++) {
    float vf=(v[i]-vmin)*istep;
    if(vf>=k)           h[k-1]++;
    else if(vf<0)       h[0]++;
    else if(!isnan(vf)) h[(int)(floor(vf))]++;
  }

  return h;
}

long ivec_count_occurrences(const int * v, long n, int val) {
  long count=0;
  while(n--) if(v[n]==val) count++;
  return count;
}

long fvec_count_occurrences(const float * v, long n, float val) {
  long count=0;
  while(n--) if(v[n]==val) count++;
  return count;
}


long fvec_count_lt(const float * v, long n, float val) {
  long count=0;
  while(n--) if(v[n]<val) count++;
  return count;
}

long ivec_count_lt(const int * v, long n, int val) {
  long count=0;
  while(n--) if(v[n]<val) count++;
  return count;
}

long fvec_count_gt(const float * v, long n, float val) {
  long count=0;
  while(n--) if(v[n]>val) count++;
  return count;
}

long ivec_count_gt(const int * v, long n, int val) {
  long count=0;
  while(n--) if(v[n]>val) count++;
  return count;
}

long fvec_count_inrange(const float * v, long n, float vmin, float vmax) {
  long count=0;
  while(n--) if(v[n]>=vmin && v[n]<vmax) count++;
  return count;
}

long ivec_count_inrange(const int * v, long n, int vmin, int vmax) {
  long count=0;
  while(n--) if(v[n]>=vmin && v[n]<vmax) count++;
  return count;
}

long fvec_count_nan (const float * v, long n)
{
  long i, nnan = 0;
  for (i = 0 ; i < n ; i++)
    if (isnan(v[i]))
      nnan++;

  return nnan;
}

long fvec_count_nonfinite (const float * v, long n)
{
  long i, nnf = 0;
  for (i = 0 ; i < n ; i++)
    if (!isfinite(v[i]))
      nnf++;

  return nnf;
}

long fvec_count_0 (const float *val, long n)
{
  int n0 = 0;
  while (n--)
    if (val[n] == 0)
      n0++;
  return n0;
}

/*-------------------------------------------------------------*/
/* mappings                                                    */
/*-------------------------------------------------------------*/


int *ivec_repeat_with_inc(const int *a,int n,
                          int nrepeat, int inc) {
  int *ret=ivec_new(nrepeat*n); 
  int i;
  for(i=0;i<nrepeat;i++) {
    ivec_cpy(ret+i*n, a, n); 
    ivec_add_scalar(ret+i*n, n, i*inc);
  }
  return ret;
}


void fvec_splat_add(const float *a,int n,
                    const int *assign,float *accu) {
  int i;
  for(i=0;i<n;i++) 
    accu[assign[i]] += a[i];
}


void fvec_isplat_add(const float *a,int n,
                     const int *assign,float *accu) {
  int i;
  for(i=0;i<n;i++) 
    accu[i] += a[assign[i]];
  
}

void fvec_map(const float *src,const int *map,int n,float *dest) {
  int i;
  for(i=0;i<n;i++) 
    dest[i]=src[map[i]];
}

void ivec_map (const int *src,const int *map,int n,int *dest) {
  int i;
  for(i=0;i<n;i++) 
    dest[i]=src[map[i]];
}

void fvec_imap(const float *src,const int *imap,int n,float *dest) {
  int i;
  for(i=0;i<n;i++) 
    dest[imap[i]]=src[i];
}


int ivec_hash(const int * v, long n) {
  unsigned int *k=(unsigned int*)v;
  unsigned int b    = 378551;
  unsigned int a    = 63689;
  unsigned int hash = 0;
  int i;
  for(i = 0; i < n; i++) {
    hash = hash * a + k[i];
    a    = a * b;
  }
  return hash;
}

void ivec_replace(int * v, long n,int val,int replace_val) {
  while(n--) if(v[n]==val) v[n]=replace_val;  
}


void ivec_accumulate_slices(const int *v,int *sl,int n) {
  int i;
  int accu=0,j=0;
  for(i=0;i<n;i++) {
    while(j<sl[i]) 
      accu+=v[j++];
    sl[i]=accu;
  }
}



/*---------------------------------------------------------------------------*/
/* Input/Output functions                                                    */
/*                                                                           */
/* To avoid repeating too much code, many functions that are the same        */
/* for ivec, fvec, bvec, etc. are implemented with prefix xvec_ and          */
/* take as 1st argument unitsize, ther size of one element of the type       */
/* at hand                                                                   */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/* Input functions                                                           */
/*---------------------------------------------------------------------------*/

static long xvecs_fsize(long unitsize, const char * fname, int *d_out, int *n_out)
{
  int d, ret; 
  long nbytes;

  *d_out = -1;
  *n_out = -1;

  FILE * f = fopen (fname, "r");
  
  if(!f) {
    fprintf(stderr, "xvecs_fsize %s: %s\n", fname, strerror(errno));
    return -1;
  }
  /* read the dimension from the first vector */
  ret = fread (&d, sizeof (d), 1, f);
  if (ret == 0) { /* empty file */
    *n_out = 0;
    return ret;
  }
  
  fseek (f, 0, SEEK_END);
  nbytes = ftell (f);
  fclose (f);
  
  if(nbytes % (unitsize * d + 4) != 0) {
    fprintf(stderr, "xvecs_size %s: weird file size %ld for vectors of dimension %d\n", fname, nbytes, d);
    return -1;
  }

  *d_out = d;
  *n_out = nbytes / (unitsize * d + 4);
  return nbytes;
}


long fvecs_fsize (const char * fname, int *d_out, int *n_out) {
  return xvecs_fsize(sizeof(float), fname, d_out, n_out);
}

long ivecs_fsize (const char * fname, int *d_out, int *n_out)
{
  return xvecs_fsize (sizeof(int), fname, d_out, n_out);
}

long bvecs_fsize (const char * fname, int *d_out, int *n_out)
{
  return xvecs_fsize (sizeof(unsigned char), fname, d_out, n_out);
}

long lvecs_fsize (const char * fname, int *d_out, int *n_out)
{
  return xvecs_fsize (sizeof(long long), fname, d_out, n_out);
}




int fvecs_new_read (const char *fname, int *d_out, float **vf_out)
{

  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_new_read: could not open %s\n", fname);
    perror ("");
    *d_out=-1;
    vf_out=NULL;
    return -1;
  }

  int n=fvecs_new_fread_max(f,d_out,vf_out,-1);

  fclose (f);
  
  return n;
}


int fvecs_new_fread_max (FILE *f, int *d_out, float **vf_out, long nmax)
{
  long n;
  int d = -1;

  float *vf=NULL;
  long nalloc=0;
  
  for(n=0; nmax<0 || n<nmax ; n++) {
    int new_d;

    if (fread (&new_d, sizeof (int), 1, f) != 1) {
      if (feof (f))
        break;
      else {
        perror ("fvecs_new_read error 1");
        goto error;
      }
    }

    if (n == 0)
      d = new_d;
    else if (new_d != d) {
      fprintf (stderr, "fvecs_new_read: non-uniform vectors sizes: new_d=%d d=%d\n",new_d,d);
      goto error;
    }

    if((n+1)*d>nalloc) {
      do {
        nalloc=nalloc==0 ? 1024 : nalloc*2;
      } while((n+1)*d>nalloc);
      vf=realloc(vf,nalloc*sizeof(float));
      assert(vf || !"out of memory!");
    }

    if (fread (vf+n*d, sizeof (float), d, f) != d) {
      fprintf (stderr, "fvecs_new_read error 3");
      goto error;
    }

  }
  vf=realloc(vf,n*d*sizeof(float));

  *vf_out = vf;
  *d_out = d;
  return n;
 error: 
  free(vf); 
  *vf_out = NULL;
  *d_out = -1;
  return -1;
}

int bvecs_new_from_siftgeo(const char *fname, 
			    int *d_out, unsigned char **v_out,
			    int *d_meta_out, float **meta_out) {

  FILE * f = fopen(fname, "r"); 
  if(!f) {
    fprintf(stderr, "could not open %s: %s\n", fname, strerror(errno));
    return -1;
  }

#define READANDCHECK(a,n) if(fread(a, sizeof(*(a)), n, f) != n) {fprintf(stderr, "weird format in %s\n", fname); goto err; }
  
  int n = 0, na = 0;
  int d = -1;
  unsigned char * v = NULL; 
  float * meta = NULL;

  for(;;) {
    float buf[9];
    int header_read = fread(buf, 1, sizeof(buf), f);

    if(header_read == 0 && feof(f)) break;
    else if(header_read != sizeof(buf)) {
      fprintf(stderr, "bvecs_new_from_siftgeo: error in point header in %s\n", fname); 
      goto err;
    }

    int d2; 
    READANDCHECK(&d2, 1)
    if(n == 0) {
      if(!(d2 >= 0 && d2 < 100000)) {
        fprintf(stderr, "bvecs_new_from_siftgeo: weird descriptor dim %d in %s\n", d2, fname); 
        goto err;        
      }        
      d = d2;
    }
    else if(d2 != d) {
      fprintf(stderr, "bvecs_new_from_siftgeo: weird dim in %s (expect %d found %d)\n", fname, d, d2); 
      goto err;
    } 
      
    if(n >= na) {
      na = na == 0 ? 512 : na*3/2;
      if(meta_out) meta = realloc(meta, sizeof(float) * 9 * na);
      v = realloc(v, d * sizeof(float) * na); 
    }
    
    if(meta_out) memcpy(meta + 9 * n, buf, sizeof(buf));
    
    READANDCHECK(v + d * n, d);
    
    n++;
  }
#undef READANDCHECK

  fclose(f); 

  *d_out = d; 
  *v_out = v; 
  if(meta_out) *meta_out = meta;
  if(d_meta_out) *d_meta_out = 9;
  return n;

  err: 
  fclose(f); 
  free(v); 
  free(meta); 
  *d_out = -1; 
  *v_out = NULL; 
  if(meta_out) *meta_out = NULL;
  if(d_meta_out) *d_meta_out = -1;
  return -1;
}




int fvecs_new_read_sparse (const char *fname, int d, float **vf_out) {
  float *vf=NULL;
  long n=0,na=0;
  float *vals=fvec_new(d);
  int *idx=ivec_new(d);
  
  FILE *f = fopen (fname, "r");
#define E(msg) {                                                \
  fprintf (stderr, "fvecs_new_read_sparse %s: " msg , fname);   \
  perror ("");                                                  \
  free(vf); free(vals); free(idx);                              \
  return -1;                                                    \
}
  if (!f) E("");
  
  while(!feof(f)) {
    int nz,ret,nz2;
    ret=fread(&nz,sizeof(int),1,f);
    if(ret!=1) {
      if(feof(f)) break;
      E("err 1");
    }
    if(fread(idx,sizeof(int),nz,f)!=nz) E("err 2");
    if(fread(&nz2,sizeof(int),1,f)!=1) E("err 3");
    if(nz!=nz2) E("err 4");
    if(fread(vals,sizeof(float),nz,f)!=nz) E("err 5");
    
    if(n>=na) {
      na=(na+1)*3/2;
      vf=realloc(vf,na*sizeof(float)*d);
    }
    
    float *dense=spfvec_to_fvec (idx,vals,nz,d);
    memcpy(vf+n*d,dense,sizeof(float)*d);
    free(dense);
    
    n++;       
  }
#undef E
  free(vals);
  free(idx);
  fclose(f);
  *vf_out=vf;
  return n;
}


int fvecs_new_mmap (const char *fname, int *d_out, float **vf) {
  int fd = open (fname, O_RDONLY);
  if (fd < 0) {
    fprintf(stderr,"[fi]vecs_new_mmap: could not open %s",fname);
    perror ("");
    *d_out=-1; *vf=NULL;
    return -1;
  }
  /* map whole file */
  struct stat sb;
  fstat (fd, &sb);
  void *mm = mmap (NULL, sb.st_size, PROT_READ, MAP_SHARED, fd, 0);
  if (mm == MAP_FAILED) {
    fprintf(stderr,"[fi]vecs_new_mmap: could not mmap %s",fname);
    perror ("");
    close(fd);
    *d_out=-1; *vf=NULL;
    return -1;
  }
  int d=(*(int*)mm);
  *d_out=d;
  *vf=(void*)(((int*)mm));
  long unitsz=sizeof(int)+sizeof(float)*d;
  assert(sb.st_size % unitsz == 0);
  return sb.st_size / unitsz;
}

int ivecs_new_mmap (const char *fname, int *d_out, int **vf) {
  return fvecs_new_mmap(fname,d_out,(float**)vf);
}




int fvecs_read (const char *fname, int d, int n, float *a)
{
  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_read: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  long i;
  for (i = 0; i < n; i++) {
    int new_d;

    if (fread (&new_d, sizeof (int), 1, f) != 1) {
      if (feof (f))
        break;
      else {
        perror ("fvecs_read error 1");
        fclose(f);
        return -1;
      }
    }

    if (new_d != d) {
      fprintf (stderr, "fvecs_read error 2: unexpected vector dimension\n");
      fclose(f);
      return -1;
    }

    if (fread (a + d * (long) i, sizeof (float), d, f) != d) {
      fprintf (stderr, "fvecs_read error 3\n");
      fclose(f);
      return -1;
    }
  }
  fclose (f);

  return i;
}


int b2fvecs_read (const char *fname, int d, int n, float *v)
{
  int n_new; 
  int d_new;
  bvecs_fsize (fname, &d_new, &n_new);
  assert (d_new == d);
  assert (n <= n_new);

  FILE * f = fopen (fname, "r");
  assert (f || "b2fvecs_read: Unable to open the file");
  b2fvecs_fread (f, v, n);
  fclose (f);
  return n;
}


int fvecs_read_txt (const char *fname, int d, int n, float *v)
{
  long i, ret;
  FILE * f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "fvecs_read_txt: could not open %s\n", fname);
    perror ("");
    return -1;
  }

  for (i = 0; i < n * d ; i++) {
    ret = fscanf (f, "%f", v + i);
    if (ret != 1) {
      if (feof (f))
	break;
      else {
        perror ("fvecs_read_txt error 1");
        /* TODO free mem */
        return -1;
      }
    }
  }
  fclose (f);

  return i / d;
}




static int xvec_fread (long unit_size, FILE * f, void * v, int d_alloc)
{
  int d;
  int ret = fread (&d, sizeof (int), 1, f);

  if (feof (f))
    return 0;

  if (ret != 1) {
    perror ("# xvec_fread error 1");
    return -1;
  }

  if (d < 0 || d > d_alloc) {
    fprintf(stderr, "xvec_fread: weird vector size (expect %d found %d)\n", d_alloc, d);
    return -1;
  }

  ret = fread (v, unit_size, d, f);
  if (ret != d) {
    perror ("# xvec_fread error 2");
    return -1;
  }

  return d;
}


int fvec_fread (FILE * f, float * v, int d_alloc)
{
  return xvec_fread(sizeof(float), f, v, d_alloc);
}

int ivec_fread (FILE * f, int * v, int d_alloc)
{
  return xvec_fread(sizeof(int), f, v, d_alloc);
}

int bvec_fread (FILE * f, unsigned char * v, int d_alloc)
{
  return xvec_fread(sizeof(unsigned char), f, v, d_alloc);
}

int lvec_fread (FILE * f, long long * v, int d_alloc)
{
  return xvec_fread(sizeof(long long), f, v, d_alloc);
}

int fvec_read (const char *fname, int d, float *a, int o_f) {
  FILE * f = fopen(fname, "r");

  if(!f) {
    fprintf(stderr, "cannot open %s: %s\n", fname, strerror(errno));
    return -1;
  }
  
  int ret = fseek(f, (4 * d + 4) * o_f, SEEK_SET);
  if(ret < 0) {
    fprintf(stderr, "seek error in %s: %s\n", fname, strerror(errno));
    fclose(f);
    return -1;
  }

  ret = fvec_fread(f, a, d);

  fclose(f);

  return ret;
}





static long xvecs_fread (long unit_size, FILE * f, void * v, long n, int d_alloc)
{
  long i = 0, d = -1, ret;
  for (i = 0 ; i < n ; i++) {
    if (feof (f))
      break;

    ret = xvec_fread (unit_size, f, (char *) v + unit_size * i * d, d_alloc);

    if (ret == 0)  /* eof */
      break;

    if (ret == -1) /* err */
      return 0;

    if (i == 0) {
      d = ret;
      if(d != d_alloc) 
	fprintf(stderr, "xvecs_fread: warn allocated d = %d, found d = %ld\n",
		d_alloc, d);
    }
      

    if (d != ret) {
      perror ("# xvecs_fread: dimension of the vectors is not consistent\n");
      return 0;
    }
  }
  return i;
}


long fvecs_fread (FILE * f, float * v, long n, int d_alloc)
{
  return xvecs_fread(sizeof(float), f, v, n, d_alloc);
}

long ivecs_fread (FILE * f, int * v, long n, int d_alloc)
{
  return xvecs_fread(sizeof(int), f, v, n, d_alloc);
}

long bvecs_fread (FILE * f, unsigned char * v, long n, int d_alloc)
{
  return xvecs_fread(sizeof(unsigned char), f, v, n, d_alloc);
}

long lvecs_fread (FILE * f, long long * v, long n, int d_alloc)
{
  return xvecs_fread(sizeof(long long), f, v, n, d_alloc);
}



/* The behavior of xvecs_new_read in not completely consistent 
   with the one of fvecs_new_read (can not read stream)         */

static int xvecs_new_read (long unitsize, const char *fname, int *d_out, void **v_out)
{
  int n;
  int d;
  xvecs_fsize (unitsize, fname, &d, &n);
  unsigned char * v = malloc(unitsize * n * d);
  FILE * f = fopen (fname, "r");
  if(!f) {
    fprintf(stderr, "xvecs_new_fread: cannot open %s: %s\n", 
	    fname, strerror(errno));
    *d_out = -1; 
    *v_out = NULL;
    return -1;
  }
  xvecs_fread (unitsize, f, v, n, d);
  fclose (f);
  *d_out = d; 
  *v_out = v;
  return n;
}

int bvecs_new_read (const char *fname, int *d_out, unsigned char **v_out) {
  return xvecs_new_read(sizeof(unsigned char), fname, d_out, (void**)v_out);
}


int lvecs_new_read (const char *fname, int *d_out, long long **v_out) {
  return xvecs_new_read(sizeof(long long), fname, d_out, (void**)v_out);
}


int b2fvecs_new_read (const char *fname, int *d_out, float **v_out)
{
  int n;
  int d;
  bvecs_fsize (fname, &d, &n);
  float * v = fvec_new ((long) n * (long) d);
  FILE * f = fopen (fname, "r");
  assert (f || "bvecs_new_read: Unable to open the file");
  b2fvecs_fread (f, v, n);
  fclose (f);
  
  *v_out = v;
  *d_out = d;
  
  return n;
}


int fvec_fread_raw (FILE *f, float * v, long n) 
{
  long ret = fread (v, sizeof (*v), n, f);
  if (ret != n) {
    free(v);
    perror ("# fvec_fread_raw error");
    return -1;
  }
  return n;
}


int ivec_fread_raw (FILE *f, int * v, long n) 
{
  long ret = fread (v, sizeof (*v), n, f);
  if (ret != n) {
    free(v);
    perror ("# ivec_fread_raw error");
    return -1;
  }
  return n;
}


int bvec_fread_raw (FILE *f, unsigned char * v, long n) 
{
  long ret = fread (v, sizeof (*v), n, f);
  if (ret != n) {
    free(v);
    perror ("# bvec_fread error");
    return -1;
  }
  return n;
}


float * fvec_new_fread_raw(FILE * f, long d) 
{
  float *v = fvec_new(d);

  long ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    free(v);
    perror ("# fvec_fread error 2");
    return NULL;
  }
  return v;
}


int * ivec_new_fread_raw(FILE * f, long d) 
{
  int * v = ivec_new(d);

  long ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    free(v);
    perror ("# fvec_fread error 2");
    return NULL;
  }
  return v;
}

unsigned char *bvec_new_fread_raw(FILE * f, long d) 
{
  unsigned char *v = bvec_new(d);

  long ret = fread (v, sizeof (*v), d, f);
  if (ret != d) {
    free(v);
    perror ("# fvec_fread error 2");
    return NULL;
  }
  return v;
}


float * fvec_new_read_raw(const char * fname, long d) 
{
  FILE * f = fopen (fname, "r");
  float * v = fvec_new_fread_raw (f, d);
  fclose (f);
  return (v);
}

int * ivec_new_read_raw(const char * fname, long d) 
{
  FILE * f = fopen (fname, "r");
  int * v = ivec_new_fread_raw (f, d);
  fclose (f);
  return (v);
}


unsigned char *bvec_new_read_raw(const char * fname, long d) 
{
  FILE * f = fopen (fname, "r");
  unsigned char * v = bvec_new_fread_raw (f, d);
  fclose (f);
  return (v);
}



long b2fvecs_fread (FILE * f, float * v, long n)
{
  long i = 0, d = -1, ret;
  for (i = 0 ; i < n ; i++) {
    if (feof (f))
      break;

    ret = b2fvec_fread (f, v + i * d);
    if (ret == 0)  /* eof */
      break;

    if (ret == -1)
      return 0;

    if (i == 0)
      d = ret;

    if (d != ret) {
      perror ("# b2fvecs_fread: dimension of the vectors is not consistent\n");
      return 0;
    }
  }
  return i;
}


int ivecs_write_txt (const char * fname, int d, int n, const int *v)
{ 
  int i, j, ret = 0;
  FILE * fo = fopen (fname, "w");
  if (!fo) {
    perror ("ivecs_write_txt: cannot open file");
    return -1;
  }

  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < d ; j++)
      fprintf (fo, "%d ", v[i * d + j]);
    ret += fprintf (fo, "\n");
  }

  return ret;
}


int ivecs_new_fread_max (FILE *f, int *d_out, int **vi_out, long nmax)
{
  
  assert(sizeof(int)==sizeof(float)); /* that's all what matters */

  return fvecs_new_fread_max(f,d_out,(float**)vi_out,nmax);
}


int ivecs_new_read (const char *fname, int *d_out, int **vi_out)
{

  FILE *f = fopen (fname, "r");
  if (!f) {
    fprintf (stderr, "ivecs_new_read: could not open %s\n", fname);
    perror ("");
    *d_out=-1; *vi_out=NULL;
    return -1;
  }
    
  int n=ivecs_new_fread_max(f,d_out,vi_out,-1);
 
  fclose(f);
  
  return n;
}


int *ivec_new_read(const char *fname, int *d_out) {
  int d;
  long n;
  int *vi;
  n = ivecs_new_read(fname,&d,&vi);
  if (n<0) 
    return NULL;
  assert (n==1);
  if (d_out) 
    *d_out=d;
  return vi;
}

int b2fvec_fread (FILE * f, float * v)
{
  int d, j;
  int ret = fread (&d, sizeof (int), 1, f);
  if (feof (f))
    return 0;

  if (ret != 1) {
    perror ("# bvec_fread error 1");
    return -1;
  }

  unsigned char * vb = (unsigned char *) malloc (sizeof (*vb) * d);

  ret = fread (vb, sizeof (*vb), d, f);
  if (ret != d) {
    perror ("# bvec_fread error 2");
    return -1;
  }
  for (j = 0 ; j < d ; j++)
    v[j] = vb[j];
  free (vb);
  return d;
}


/*---------------------------------------------------------------------------*/
/* Output functions                                                          */
/*---------------------------------------------------------------------------*/

int fvec_fwrite (FILE *fo, const float *v, int d) 
{
  int ret;
  ret = fwrite (&d, sizeof (int), 1, fo);
  if (ret != 1) {
    perror ("fvec_fwrite: write error 1");
    return -1;
  }
  ret = fwrite (v, sizeof (float), d, fo);
  if (ret != d) {
    perror ("fvec_fwrite: write error 2");
    return -1;
  }  
  return 0;
}


int fvec_fwrite_raw(FILE *fo, const float *v, long d) 
{
  long ret = fwrite (v, sizeof (*v), d, fo);
  if (ret != d) {
    perror ("fvec_fwrite_raw: write error 2");
    return -1;
  }  
  return 0;
}


int ivec_fwrite_raw(FILE *fo, const int *v, long d) 
{
  long ret = fwrite (v, sizeof (*v), d, fo);
  if (ret != d) {
    perror ("ivec_fwrite_raw: write error 2");
    return -1;
  }  
  return 0;
}


int bvec_fwrite_raw(FILE *fo, const unsigned char *v, long d) 
{
  long ret = fwrite (v, sizeof (*v), d, fo);
  if (ret != d) {
    perror ("ivec_fwrite_raw: write error 2");
    return -1;
  }  
  return 0;
}


int ivec_write_raw (const char *fname, const int *v, long d)
{
  FILE * f = fopen (fname, "w");
  assert (f);
  int ret = ivec_fwrite_raw (f, v, d);
  fclose (f);
  return ret;
}

int fvec_write_raw (const char *fname, const float *v, long d)
{
  FILE * f = fopen (fname, "w");
  assert (f);
  int ret = fvec_fwrite_raw (f, v, d);
  fclose (f);
  return ret;
}

int bvec_write_raw (const char *fname, const unsigned char *v, long d)
{
  FILE * f = fopen (fname, "w");
  assert (f);
  int ret = bvec_fwrite_raw (f, v, d);
  fclose (f);
  return ret;
}

int fvecs_fwrite (FILE *fo, int d, int n, const float *vf)
{
  int i;
  /* write down the vectors as fvecs */
  for (i = 0; i < n; i++) {
    if(fvec_fwrite(fo, vf+i*d, d)<0) 
      return i;
  }
  return n;
}



int fvecs_write (const char *fname, int d, int n, const float *vf)
{
  FILE *fo = fopen (fname, "w");
  if (!fo) {
    fprintf(stderr, "fvecs_write: cannot open %s for writing", fname);
    perror ("");
    return -1;
  }

  int ret = fvecs_fwrite (fo, d, n, vf);

  fclose (fo);
  return ret;
}


int fvecs_write_txt (const char * fname, int d, int n, const float *vf)
{ 
  int i, j, ret = 0;
  FILE * fo = fopen (fname, "w");
  if (!fo) {
    perror ("fvecs_write_txt: cannot open file");
    return -1;
  }

  for (i = 0 ; i < n ; i++) {
    for (j = 0 ; j < d ; j++)
      fprintf (fo, "%f ", vf[i * d + j]);
    ret += fprintf (fo, "\n");
  }

  return ret;
}


int ivec_fwrite (FILE *f, const int *v, int d)
{
  int ret = fwrite (&d, sizeof (d), 1, f);
  if (ret != 1) {
    perror ("ivec_fwrite: write error 1");
    return -1;
  }

  ret = fwrite (v, sizeof (*v), d, f);
  if (ret != d) {
    perror ("ivec_fwrite: write error 2");
    return -2;
  }
  return 0;
}




int ivecs_fwrite(FILE *f, int d, int n, const int *v)
{
  int i;
  for (i = 0 ; i < n ; i++) {    
    ivec_fwrite (f, v, d);
    v+=d;
  }
  return n;
}


int ivecs_write (const char *fname, int d, int n, const int *v)
{
  int ret = 0;
  FILE *f = fopen (fname, "w");
  if (!f) {
    perror ("ivecs_write");
    return -1;
  }

  ret = ivecs_fwrite (f, d, n, v);

  fclose (f);
  return ret;
}


int bvec_fwrite (FILE *f, const unsigned char *v, int d)
{
  int ret = fwrite (&d, sizeof (d), 1, f);
  if (ret != 1) {
    perror ("bvec_fwrite: write error 1");
    return -1;
  }

  ret = fwrite (v, sizeof (*v), d, f);
  if (ret != d) {
    perror ("bvec_fwrite: write error 2");
    return -2;
  }
  return 0;
}




int bvecs_fwrite(FILE *f, int d, int n, const unsigned char *v)
{
  int i;
  for (i = 0 ; i < n ; i++) {    
    bvec_fwrite (f, v, d);
    v+=d;
  }
  return n;
}


int bvecs_write (const char *fname, int d, int n, const unsigned char *v)
{
  int ret = 0;
  FILE *f = fopen (fname, "w");
  if (!f) {
    perror ("bvecs_write");
    return -1;
  }

  ret = bvecs_fwrite (f, d, n, v);

  fclose (f);
  return ret;
}



void fvec_print (const float * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%g ", v[i]);
  printf ("]\n");
}


void ivec_print (const int * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%d ", v[i]);
  printf ("]\n");
}


void bvec_print (const unsigned char * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%d ", (int) v[i]);
  printf ("]\n");
}


void lvec_print (const long long * v, int n)
{
  int i;
  printf ("[");
  for (i = 0 ; i < n ; i++)
    printf ("%lld ", (long long) v[i]);
  printf ("]\n");
}


void fvec_fprintf (FILE * f, const float *v, int n, const char *fmt)
{
  int i;
  if (fmt == NULL)
    fmt = "%f ";

  for (i = 0; i < n; i++)
    fprintf (f, fmt, v[i]);
}


void ivec_fprintf (FILE * f, const int *v, int n, const char *fmt)
{
  int i;
  if (fmt == NULL)
    fmt = "%d ";

  for (i = 0; i < n; i++)
    fprintf (f, fmt, v[i]);
}


/*---------------------------------------------------------------------------*/
/* Vector cast                                                              */
/*---------------------------------------------------------------------------*/


float * ivec2fvec (const int * v, long n)
{
  long i;
  float * vo = fvec_new (n);
  for (i = 0 ; i < n ; i++)
    vo[i] = (float) v[i];

  return vo;
}


float * bvec2fvec (const unsigned char * v, long n)
{
  long i;
  float * vo = fvec_new (n);
  for (i = 0 ; i < n ; i++)
    vo[i] = (float) v[i];

  return vo;
}


void bvectofvec (const unsigned char * v, float * vb, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    vb[i] = (float) v[i];
}

void fvectodvec (const float *a, double *b, long n) {
  long i; 
  for (i = 0 ; i < n ; i++)
    b[i] = a[i];  
}


/*---------------------------------------------------------------------------*/
/* Elementary operations                                                     */
/*---------------------------------------------------------------------------*/

void fvec_0(float * v, long n)
{
  memset (v, 0, n * sizeof (*v));
}

void fvec_nan (float * v, long n) {
  memset (v, -1, n * sizeof (*v));
} 


void ivec_0(int * v, long n)
{
  memset (v, 0, n * sizeof (*v));
}

int fvec_all_0 (const float * v, long n) {
  while(n--) if(v[n]!=0) return 0;
  return 1;
}

int fvec_all_ge0 (const float * v, long n) {
  while(n--) if(v[n]<0) return 0;
  return 1;
}

int ivec_all_ge0 (const int * v, long n) {
  while(n--) if(v[n]<0) return 0;
  return 1;
}

int fvec_all_finite (const float * v, long n) {
  while(n--) if(!isfinite(v[n])) return 0;
  return 1;
}

int ivec_all_0 (const int * v, long n) {
  while(n--) if(v[n]!=0) return 0;
  return 1;
}

void fvec_set (float * v, long n, float val)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = val;
}


void ivec_set (int * v, long n, int val)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = val;
}


void ivec_cpy (int * vdest, const int * vsource, long n)
{
  memmove (vdest, vsource, n * sizeof (*vdest));
}


void fvec_cpy (float * vdest, const float * vsource, long n)
{
  memmove (vdest, vsource, n * sizeof (*vdest));
}


void bvec_cpy (unsigned char * vdest, const unsigned char * vsource, long n)
{
  memmove (vdest, vsource, n * sizeof (*vdest));
}


void fvec_incr (float * v, long n, double scal)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] += scal;
}


void fvec_decr (float * v, long n, double scal)
{
  fvec_incr(v, n, -scal);  
}


void ivec_incr (int * v, long n, int scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] += scal;
}
void ivec_decr (int * v, long n, int scal) {
  ivec_incr(v, n, -scal);  
}


void fvec_mul_by (float * v, long n, double scal)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] *= scal;
}

void fvec_div_by (float * v, long n, double scal)
{
  fvec_mul_by(v, n, 1. / scal);
}

void fvec_rdiv_by (float * v, long n, double scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] = scal / v[i];
}


void fvec_add (float * v1, const float * v2, long n)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] += v2[i];
}

void fvec_add_mul (float * v1, const float * v2, long n,double scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] += v2[i]*scal;
}

void fvec_sub (float * v1, const float * v2, long n)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] -= v2[i];
}

void fvec_rev_sub (float * v1, const float * v2, long n) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] = v2[i] - v1[i];
}

void fvec_mul (float * v1, const float * v2, long n)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] *= v2[i];
}


void fvec_div (float * v1, const float * v2, long n)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] /= v2[i];
}



void fvec_round (float * v, long n)
{
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] = (int)(v[i]+0.5);
}

void fvec_sqrt (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = sqrt (v[i]);
}


void fvec_sqr (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] =  v[i] * v[i];
}


void fvec_exp (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = exp (v[i]);
}


void fvec_log (float * v, long n)
{
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = log (v[i]);
}

void fvec_neg (float * v, long n) {
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = -v[i];
}

void fvec_ssqrt (float * v, long n) {
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = v[i]>0 ? sqrt(v[i]) : -sqrt(-v[i]);

}

void fvec_spow (float * v, long n, double scal) {
  long i;
  for (i = 0 ; i < n ; i++)
    v[i] = v[i]>0 ? pow(v[i],scal) : -pow(-v[i],scal);
}

void ivec_add (int * v1, const int * v2, long n) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] += v2[i];
}

void ivec_sub (int * v1, const int * v2, long n){
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v1[i] -= v2[i];
}

void ivec_mul_by (int * v,long n, int scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] *= scal;
}

void ivec_mod_by (int * v,long n, int scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] = v[i] % scal;
}


void ivec_add_scalar (int * v, long n, int scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] += scal;
}

void fvec_add_scalar (float * v, long n, float scal) {
  long i = 0;
  for (i = 0 ; i < n ; i++)
    v[i] += scal;
}


/*---------------------------------------------------------------------------*/
/* NaN's and non-finite elements                                             */
/*---------------------------------------------------------------------------*/



long fvec_purge_nans(float * v, long n, float replace_value) {
  long i, count=0;
  
  for(i=0;i<n;i++) if(isnan(v[i])) {
    count++;
    v[i]=replace_value;
  }
  
  return count;
}

int fvec_purge_nonfinite(float * v, long n, float replace_value) {
  long i, count=0;
  
  for(i=0;i<n;i++) if(!isfinite(v[i])) {
    count++;
    v[i]=replace_value;
  }
  
  return count;
}

long fvec_shrink_nonfinite(float * v, long n) {
  long i,j=0;
  
  for(i=0;i<n;i++) if(isfinite(v[i])) 
    v[j++]=v[i];
  return j;  
}

long fvec_index_nonfinite(float * v, long n) {
  long i;
  for(i=0;i<n;i++) if(!isfinite(v[i])) return i;
  return -1;  
}

void fvec_revert(float *v, long n) {
  long i; 
  for(i = 0; i * 2 < n; i++) {
    float tmp = v[i]; 
    v[i] = v[n - 1 - i];    
    v[n - 1 - i] = tmp; 
  }
}

void fvec_swap(float *v1, float *v2, long n) {
  long i; 
  for(i = 0; i < n; i++) {
    float tmp = v1[i]; 
    v1[i] = v2[i];    
    v2[i] = tmp; 
  }
}
  



/*---------------------------------------------------------------------------*/
/* Vector measures and statistics                                            */
/*---------------------------------------------------------------------------*/

double fvec_normalize (float * v, long n, double norm)
{
  if(norm==0) 
    return 0;

  double nr = fvec_norm (v, n, norm);

  /*  if(nr!=0)*/
  fvec_mul_by (v, n, 1. / nr);
  return nr;
}

void fvec_normalize_2stage(float * v, long n, double scal) {
  double nr = fvec_normalize (v, n, 2);
  
  if(nr == 0) return; 

  int renorm = 0;
  int i;

  for(i=0; i<n; i++) 
    if(v[i] > scal) {
      v[i] = scal; 
      renorm = 1; 
    }
  
  if(renorm) 
    fvec_normalize (v, n, 2);
    
}



int fvecs_normalize (float * v, long n, long d, double norm)
{
  int i, nNaN = 0;
  if (norm == 0) 
    return 0;

  for (i = 0 ; i < n ; i++) {
    double nr = fvec_norm (v + i * d, d, norm);
    if (nr == 0)
      nNaN++;

    fvec_mul_by (v + i * d, d, 1 / nr);
  }

  return nNaN;
}

double fvec_sum (const float * v, long n)
{
  long i;
  double s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s;
}

void fvec_cumsum(float * v, long n) {
  long i;
  double s = 0;
  for (i = 0 ; i < n ; i++) {
    s += v[i];
    v[i] = s;
  }
}

void ivec_cumsum(int *v, long n) {
  long i;
  int s = 0;
  for (i = 0 ; i < n ; i++) {
    s += v[i];
    v[i] = s;
  }
}


void fvec_cumdiff(float *v, long n) {
  long i;
  float prev = 0;
  for (i = 0 ; i < n ; i++) {
    float pprev = v[i];
    v[i] -= prev;
    prev = pprev;
  }
}

void ivec_cumdiff(int *v, long n) {
  long i;
  int prev = 0;
  for (i = 0 ; i < n ; i++) {
    int pprev = v[i];
    v[i] -= prev;
    prev = pprev;
  }
}


double fvec_mean (const float * v, long n)
{
  long i;
  double s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s / n;
}


double fvec_sum_sqr (const float * v, long n)
{
  return fvec_norm2sqr(v, n);
}


double fvec_product (const float * v, long n) {
  double p=1.0;
  while(n--) p*=v[n];
  return p;
}


long long ivec_product (const int * v, long n) {
  long long p=1;
  while(n--) p*=v[n];
  return p;  
}


long long ivec_sum (const int * v, long n)
{
  long i;
  long long s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s;
}


long long ivec_mean (const int * v, long n)
{
  long i;
  long long s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i];

  return s / n;
}


long long ivec_sum_sqr (const int * v, long n)
{
  long i;
  long long s = 0;
  for (i = 0 ; i < n ; i++)
    s += v[i]*(long long)v[i];

  return s;
}


double fvec_norm (const float * v, long n, double norm)
{
  if(norm==0) return n;

  long i;
  double s = 0;

  if(norm==1) {
    for (i = 0 ; i < n ; i++)
      s += fabs(v[i]);
    return s;
  }

  if(norm==2) {
    for (i = 0 ; i < n ; i++)
      s += v[i]*v[i];
    
    return sqrt(s);
  } 

  if(norm==-1) {
    for (i = 0 ; i < n ; i++)
      if(fabs(v[i])>s) s=fabs(v[i]);
    return s;
  }

  for (i = 0 ; i < n ; i++) {
    s += pow (v[i], norm);
  }


  return pow (s, 1 / norm);
}


double fvec_norm2sqr (const float * v, long n) {
  double s=0;
  long i;
  for (i = 0 ; i < n ; i++)
    s += v[i] * v[i];
  return s;
}

long fvec_nz (const float * v, long n)
{
  long i, nz = 0;
  for (i = 0 ; i < n ; i++)
    if (v[i] != 0)
      nz++;

  return nz;
}

long ivec_nz (const int * v, long n)
{
  long i, nz = 0;
  for (i = 0 ; i < n ; i++)
    if (v[i] != 0)
      nz++;

  return nz;
}


int fvec_find (const float *v, int n, int ** nzpos_out)
{
  int nz = fvec_nz (v, n);
  int * nzpos = ivec_new (nz);
  int i, ii = 0;

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      nzpos[ii] = i;
      ii++;
    }

  *nzpos_out = nzpos;
  return nz;
}


int ivec_find (const int *v, int n, int ** nzpos_out)
{
  int nz = ivec_nz (v, n);
  int * nzpos = ivec_new (nz);
  int i, ii = 0;

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      nzpos[ii] = i;
      ii++;
    }

  *nzpos_out = nzpos;
  return nz;
}





double fvec_entropy (const float *pmf, int n)
{
  int i;
  double minusent = 0;

  for (i = 0 ; i < n ; i++)
    if (pmf[i] > 0)
      minusent += pmf[i] * log (pmf[i]);

  return - minusent / log(2);
}

/*! @brief entropy of a binary variable */
double binary_entropy (double p)
{
  if (p == 0 || p == 1)
    return 0;
  return -(p * log (p) + (1-p) * log (1-p)) / log (2);
}

double ivec_unbalanced_factor(const int *hist, long n) {
  int vw;
  double tot = 0, uf = 0;

  for (vw = 0 ; vw < n ; vw++) {
    tot += hist[vw];
    uf += hist[vw] * (double) hist[vw];
  }

  uf = uf * n / (tot * tot);

  return uf;

}

/*---------------------------------------------------------------------------*/
/* Distances                                                                 */
/*---------------------------------------------------------------------------*/

long ivec_distance_hamming (const int * v1, const int * v2, long n)
{
  long i, dis = 0; 

  for (i = 0 ; i < n ; i++)
    dis += (v1[i] == v2[i] ? 0 : 1);

  return dis;
}


double fvec_distance_L2 (const float * v1, const float * v2, long n)
{
  return sqrt (fvec_distance_L2sqr (v1, v2, n));
}


double fvec_distance_L1 (const float * v1, const float * v2, long n) {
  long i;
  double dis = 0;

  for (i = 0 ; i < n ; i++) {
    dis += fabs(v1[i] - v2[i]);
  }

  return dis;  
}

double fvec_distance_L2sqr (const float * v1, const float * v2, long n)
{
  long i;
  double dis = 0, a;

  for (i = 0 ; i < n ; i++) {
    a = (double) v1[i] - v2[i];
    dis += a * a;
  }

  return dis;  
}


double fvec_inner_product (const float * v1, const float * v2, long n)
{
  long i;
  double res = 0;
  for (i = 0 ; i < n ; i++)
    res += v1[i] * v2[i];
  return res;
}


/*---------------------------------------------------------------------------*/
/* Sparse vector handling                                                    */
/*---------------------------------------------------------------------------*/

int fvec_to_spfvec (float * v, int n, int ** idx_out, float ** v_out)
{
  int i, ii = 0;
  int nz = fvec_nz (v, n);
  int * idx = ivec_new (nz);
  float * val = fvec_new (nz);

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      idx[ii] = i;
      val[ii] = v[i];
      ii++;
    }

  *idx_out = idx;
  *v_out = val;
  return nz;
}


int ivec_to_spivec (int * v, int n, int ** idx_out, int ** v_out)
{
  int i, ii = 0;
  int nz = ivec_nz (v, n);
  int * idx = ivec_new (nz);
  int * val = ivec_new (nz);

  for (i = 0 ; i < n ; i++) 
    if (v[i] != 0) {
      idx[ii] = i;
      val[ii] = v[i];
      ii++;
    }

  *idx_out = idx;
  *v_out = val;
  return nz;
}


float * spfvec_to_fvec (int * idx, float * v, int nz, int n)
{
  int i;
  float * ret = fvec_new_0 (n);
  for (i = 0 ; i < nz ; i++) 
    if(idx[i] >= 0)  /* ignore bad bins */
      ret[idx[i]] = v[i];

  return ret;
}


int * spivec_to_ivec (int * idx, int * v, int nz, int n)
{
  int i;
  int * ret = ivec_new_0 (n);
  for (i = 0 ; i < nz ; i++)
    ret[idx[i]] = v[i];

  return ret;
}


float spfvec_inner_product (int *idx1, float *val1, int nz1, 
			    int *idx2, float *val2, int nz2)
{
  double s = 0;
  long i1 = 0, i2 = 0;

  while (1) {

    if (i1 == nz1) 
      break;

    if (i2 == nz2) 
      break;

    if (idx1[i1] == idx2[i2]) {
      s += val1[i1] * val2[i2];
      i1++;
      i2++;
    }

    else {
      if (idx1[i1] < idx2[i2])
	i1++;

      else
	i2++;
    }
  }
  return s;
}


long ivec_index(const int * v, long n,int val) {
  long i;
  for(i=0;i<n;i++) if(v[i]==val) return i;
  return -1;
}


void fvec_cpy_subvectors (const float * v, int * idx, int d, int nout, float * vout)
{
  long i;
  for (i = 0 ; i < nout ; i++)
    fvec_cpy (vout + i * d, v + (long) idx[i] * d, d);
}


void ivec_to_fvec(const int *v, float *f, long n) {
  while(n--) f[n] = v[n];
}

/* copy a subset of byte vectors and cast them to float vectors in the same time */
void b2fvec_cpy_subvectors (const unsigned char * v, int * idx, int d, int nout, float * vout)
{
  long i, j;
  for (i = 0 ; i < nout ; i++)
    for (j = 0 ; j < d ; j++)
      vout[i* (long)d+j] = v[idx[i] * (long) d+j];
}

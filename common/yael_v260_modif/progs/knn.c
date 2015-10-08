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


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <yael/vector.h>
#include <yael/machinedeps.h>
#include <yael/nn.h>

#define FMT_FVEC  0
#define FMT_IVEC  1
#define FMT_BVEC  2
#define FMT_TEXT  3


void usage (const char * cmd)
{
  printf ("Usage: %s [options]\n", cmd);
  
  printf (
	  "  Input\n"
	  "    -b filename    database to be searched (fvec file format)\n"
	  "    -bb filename   database to be searched (bvec format)\n"
	  "    -bt filename   database to be searched (text format)\n"
	  "    -q filename    query vectors (file in fvec format)\n"
	  "    -qb filename   query vectors (file in bvec format)\n"
	  "    -qt filename   query vectors (text format)\n"
	  "    -nb #          use the first n database vectors (default: all)\n"
	  "    -nq #          use the first n query vectors (default: all)\n"
          "    -d #           dimension of the vector (used with -itxt only)\n"
	  "  Search parameters\n"                  
	  "    -k #           number of nearest neighbors to retrieve (default: 10)\n"
          "    -nt #          number of threads used. Default: nt=all\n"
	  "  Ouput\n"                          
	  "    -onn idxfile   file (fmt=ivec) of k-NN indexes, 1 vec/query (default: nn.out)\n"
	  "    -onnt idxfile  same as onn but in text format\n"
	  "    -odis disfile  square distances to k-NN, 1 vec/query (default: dis.out)\n"
	  "    -odist disfile same as dis but output in text format\n"
          "  Others:\n"
          "    -verbose       better to look at this than chatting on Internet\n"
          "    -silent        shut up!\n\n"
	  "Remark: in text format, all the parameters nb, nq and d must be defined\n"
	  );
  exit (0);
}

/* read a file containing floating vectors in various formats. Check the NaN */
float * my_fvec_read (const char * fname, int fmt, int verbose, 
		      int * nptr, int * dptr)
{
  int ret = -1;
  float * v = NULL;
  int n = *nptr;
  int d = *dptr;

  /* read the vectors from the input file, and sanitize them if needed */
  switch (fmt) {
  case FMT_FVEC:
    if (d == 0 || n == 0) {/* automatically read the dimension */
      ret = fvecs_fsize (fname, &d, &n);
      fprintf (stderr, "File %s contains (%d bytes) %d vectors of dimension %d\n", 
	       fname, ret, n, d);
      assert (ret);    
    }
    v = fvec_new (n * d);
    ret = fvecs_read (fname, d, n, v);
    break;

  case FMT_TEXT:
    assert (n !=0 && d != 0);
    v = fvec_new (n * d);
    ret = fvecs_read_txt (fname, d, n, v);
    break;

  case FMT_BVEC:
    if (d == 0 || n == 0) {/* automatically read the dimension */
      ret = bvecs_fsize (fname, &d, &n);
      fprintf (stderr, "File %s contains (%d bytes) %d vectors of dimension %d\n", 
	       fname, ret, n, d);
      assert (ret);    
    }
    v = fvec_new (n * d);
    ret = b2fvecs_read (fname, d, n, v);
    break;

  default: assert (1 || "Unknown input format\n");
  }

  assert (ret >= n);

  /* Remove the Nan values */
  int nNaN = fvec_purge_nans (v, n * d, 2);
  if (nNaN > 0 && verbose >= 1)
    fprintf (stderr, "found %d NaN values in what read from %s\n", 
	     nNaN, fname);

  *dptr = d;
  *nptr = n;
  return v;
}


int main (int argc, char ** argv)
{
  int i;
  int k = 10;
  int d = 0;
  int nb = 0;
  int nq = 0;
  int nt = count_cpu();
  int verbose = 1;
  int ret = 0;

  int fmt_b = FMT_FVEC;
  int fmt_q = FMT_FVEC;
  int fmt_nn = FMT_IVEC;
  int fmt_dis = FMT_FVEC;

  const char * fb_name = NULL;    /* database filename */
  const char * fq_name = NULL;    /* query filename */
  const char * fnn_name = "nn.out";   /* nn idx filename */
  const char * fdis_name = "dis.out";  /* nn dis filename */

  if (argc == 1)
    usage (argv[0]);

  for (i = 1 ; i < argc ; i++) {
    char *a = argv[i];

    if (!strcmp (a, "-h") || !strcmp (a, "--help"))
      usage (argv[0]);
    else if (!strcmp (a, "-silence")) {
      verbose = 0;
    }
    else if (!strcmp (a, "-verbose")) {
      verbose = 2;
    }
    else if (!strcmp (a, "-k") && i+1 < argc) {
      ret = sscanf (argv[++i], "%d", &k);
      assert (ret);
    }
    else if (!strcmp (a, "-d") && i+1 < argc) {
      ret = sscanf (argv[++i], "%d", &d);
      assert (ret);
    }
    else if (!strcmp (a, "-nt") && i+1 < argc) {
      ret = sscanf (argv[++i], "%d", &nt);
      assert (ret);
    }
    else if (!strcmp (a, "-nb") && i+1 < argc) {
      ret = sscanf (argv[++i], "%d", &nb);
      assert (ret);
    }
    else if (!strcmp (a, "-nq") && i+1 < argc) {
      ret = sscanf (argv[++i], "%d", &nq);
      assert (ret);
    }
    else if (!strcmp (a, "-b") && i+1 < argc) {
      fb_name = argv[++i];
      fmt_b = FMT_FVEC;
    }
    else if (!strcmp (a, "-bb") && i+1 < argc) {
      fb_name = argv[++i];
      fmt_b = FMT_BVEC;
    }
    else if (!strcmp (a, "-bt") && i+1 < argc) {
      fb_name = argv[++i];
      fmt_b = FMT_TEXT;
    }
    else if (!strcmp (a, "-q") && i+1 < argc) {
      fq_name = argv[++i];
      fmt_q = FMT_FVEC;
    }
    else if (!strcmp (a, "-qb") && i+1 < argc) {
      fq_name = argv[++i];
      fmt_q = FMT_BVEC;
    }
    else if (!strcmp (a, "-qt") && i+1 < argc) {
      fq_name = argv[++i];
      fmt_q = FMT_TEXT;
    }
    else if (!strcmp (a, "-onn") && i+1 < argc) {
      fnn_name = argv[++i];
      fmt_nn = FMT_IVEC;
    }
    else if (!strcmp (a, "-onnt") && i+1 < argc) {
      fnn_name = argv[++i];
      fmt_nn = FMT_TEXT;
    }
    else if (!strcmp (a, "-odis") && i+1 < argc) {
      fdis_name = argv[++i];
      fmt_dis = FMT_FVEC;
    }
    else if (!strcmp (a, "-odist") && i+1 < argc) {
      fdis_name = argv[++i];
      fmt_dis = FMT_TEXT;
    }
  }

  assert (fb_name && fq_name);

  fprintf (stderr, "k = %d\nd = %d\nnt = %d\n", k, d, nt);

  if (verbose) {
    fprintf (stderr, "fb = %s  (fmt = %s)\n", fb_name, 
	     (fmt_b == FMT_FVEC ? "fvec" : (fmt_b == FMT_BVEC ? "bvec" : "txt")));
    fprintf (stderr, "fq = %s  (fmt = %s)\n", fq_name, 
	     (fmt_q == FMT_FVEC ? "fvec" : (fmt_q == FMT_BVEC ? "bvec" : "txt")));
    fprintf (stderr, "fnn = %s  (fmt = %s)\n", fnn_name, 
	     (fmt_nn == FMT_IVEC ? "ivec" : "txt"));
    fprintf (stderr, "fdis = %s  (fmt = %s)\n", fdis_name, 
	     (fmt_dis == FMT_FVEC ? "fvec" : "txt"));
  }


  /* read the input vectors for database and queries */
  float * vb = my_fvec_read (fb_name, fmt_b, verbose, &nb, &d);
  float * vq = my_fvec_read (fq_name, fmt_q, verbose, &nq, &d);


  /* Search */
  int * idx = ivec_new (k * nq);
  float * dis = fvec_new (k * nq);

  knn_full_thread (2, nq, nb, d, k, vb, vq, NULL, idx, dis, nt, NULL, NULL);
  knn_reorder_shortlist (nq, nb, d, k, vb, vq, idx, dis);

  /* write the distance output file */
  if (fmt_dis == FMT_FVEC)
    ret = fvecs_write (fdis_name, k, nq, dis);
  else if (fmt_dis == FMT_TEXT)
    ret = fvecs_write_txt (fdis_name, k, nq, dis);
  else assert (0 || "Unknow output format\n");
  assert (ret == nq);
  
  /* write the distance output file */
  if (fmt_nn == FMT_IVEC)
    ret = ivecs_write (fnn_name, k, nq, idx);
  else if (fmt_nn == FMT_TEXT)
    ret = ivecs_write_txt (fnn_name, k, nq, idx);
  else assert (0 || "Unknow output format\n");
  assert (ret == nq);
  
  free (idx);
  free (dis);
  free (vb);
  free (vq);
  return 0;
}

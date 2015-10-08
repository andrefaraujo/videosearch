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


/* Display a set of vector stored in a vector float format */

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <assert.h>


void display_help (const char * progname)
{
  fprintf (stderr, "%s reads a fvecfile on standard input\n", progname);
  exit (0);
}


int main (int argc, char **argv) 
{
  int i, j, d, ret;
  int n = INT_MAX;          /* maximum number of vectors to be read */
  int maxd = 500000000;

  FILE * fi = stdin;
  FILE * fo = stdout;

  /* read the arguments */
  for (i = 1 ; i < argc ; i++) {
    char *a = argv[i];

    if (!strcmp (a, "-n") && i + 1 < argc) {
      ret = sscanf (argv[++i], "%d", &n);
      assert (ret == 1);
    } 
    else if (!strcmp (a, "-maxd") && i + 1 < argc) {
      ret = sscanf (argv[++i], "%d", &maxd);
      assert (ret == 1);
    } 
    else {
      fprintf (stderr, "could not parse argument %s\n", a);
      display_help (argv[0]);
    }
  }

  /* Read the values while there are some */
  int * v = malloc (sizeof (*v) * maxd);
  i = 0;
  while (!feof (fi) && i < n) {
    ret = fread (&d, sizeof (d), 1, fi);

    if (ret == 0)
      break;

    assert (d < maxd); 
    ret = fread (v, sizeof (*v), d, fi);
    assert (ret == d);
    
    fprintf (fo, "[");
    for (j = 0 ; j < d ; j++)
      fprintf (fo, "%d ", v[j]);
    fprintf (fo, "]\n");
    i++;
  }

  free (v);
  fprintf (stderr, "found %d vectors\n", i);
  return 0;
}



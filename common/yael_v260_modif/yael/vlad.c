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


#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "vlad.h"
#include "nn.h"
#include "vector.h"
#include "sorting.h"

void vlad_compute(int k, int d, const float *centroids, 
                  int n, const float *v,
                  float *desc) {
  
  int i,j;
  int *assign=ivec_new(n);
 
  nn(n,k,d,centroids,v,assign,NULL,NULL);

  fvec_0(desc,k*d);
      
  for(i=0;i<n;i++) {
    for(j=0;j<d;j++) 
      desc[assign[i]*d+j]+=v[i*d+j]-centroids[assign[i]*d+j];
  }      

  free(assign);
}


void vlad_compute_weighted(int k, int d, const float *centroids, 
                           int n, const float *v, const float *weights, 
                           float *desc) {
  
  int i,j;
  int *assign=ivec_new(n);
 
  nn(n,k,d,centroids,v,assign,NULL,NULL);

  fvec_0(desc,k*d);
      
  for(i=0;i<n;i++) {
    float w=weights[i];
    for(j=0;j<d;j++) 
      desc[assign[i]*d+j] += (v[i*d+j]-centroids[assign[i]*d+j])*w;
  }      

  free(assign);
}


void vlad_compute_subsets(int k, int d, const float *centroids, 
                          int n, const float *v,
                          int n_subset,
                          const int *subset_indexes, 
                          const int *subset_ends,
                          float *desc) {
  int j;
  int *assign=ivec_new(n);
 
  nn(n,k,d,centroids,v,assign,NULL,NULL);

  fvec_0(desc,k*d*n_subset);
      
  int ss,ss_begin=0;
  for(ss=0;ss<n_subset;ss++) {
    float *descss=desc+ss*k*d;
    int ss_end=subset_ends[ss],ii;
    for(ii=ss_begin;ii<ss_end;ii++) {
      int i=subset_indexes[ii];
      for(j=0;j<d;j++) 
        descss[assign[i]*d+j] += v[i*d+j]-centroids[assign[i]*d+j];
    }
    ss_begin=ss_end;
  }

  free(assign);
  
}


void bof_compute_subsets(int k, int d, const float *centroids, 
                         int n, const float *v,
                         int n_subset,
                         const int *subset_indexes, 
                         const int *subset_ends,
                         float *desc) 
{
  int *assign=ivec_new(n);
 
  nn(n,k,d,centroids,v,assign,NULL,NULL);

  fvec_0(desc,k*n_subset);
      
  int ss,ss_begin=0;
  for(ss=0;ss<n_subset;ss++) {
    float *descss=desc+ss*k;
    int ss_end=subset_ends[ss],ii;
    for(ii=ss_begin;ii<ss_end;ii++) {
      int i=subset_indexes[ii];
      descss[assign[i]] ++;
    }
    ss_begin=ss_end;
  }

  free(assign);
}


void bof_compute (int k, int d, const float *centroids, 
		  int n, const float *v, int *desc)
{
  int i;
  int *assign=ivec_new(n);
 
  nn(n,k,d,centroids,v,assign,NULL,NULL);

  ivec_0(desc,k);

  for(i=0;i<n;i++)
    desc[assign[i]]++;

  free(assign);
}


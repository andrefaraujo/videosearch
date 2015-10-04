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


#ifndef __hkm_h
#define __hkm_h


/*--------------------------------------------------------------
 * hierarchical clustering
 --------------------------------------------------------------*/

/*! the structure used for the quantization */
typedef struct hkm_s {
  int nlevel;            /* number of levels */
  int bf;                /* the branching factor */
  int k;                 /* the number of leaves (bf^nlevel) */
  int d;                 /* dimension of the input vectors */
  float ** centroids;    /* centroids for all levels */
} hkm_t;


/*! create/delete a hierarchical quantizer structure.
   nlevel is the number of levels in the tree and bf the branching factor */
hkm_t * hkm_learn (int n, int d, int nlevel, int bf, 
		   const float * v, int nb_iter_max, int nt, int verbose, 
		   int ** clust_assign_out);
		 
void hkm_delete (hkm_t * hkm);

/*! Quantization usign the hierarchical clustering */
void hkm_quantize (const hkm_t * hkm, int n, const float * v, int * idx);

/*! I/O function for hkm */
void hkm_write (const char * filename, const hkm_t * hkm);
hkm_t * hkm_read (const char * filename);

/*! retrieve the centroids from a particular level */
float * hkm_get_centroids (const hkm_t * hkm, int l, int no);

#endif

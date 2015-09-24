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

/*
 * sift_extractor.cc
 *
 * Library to extract SIFT keypoints and get descriptors, based on VLFEAT.
 *
 *  Created on: May 22, 2014
 *      Author: andrefaraujo
 */

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <string>

extern "C" {
#include "../../common/vlfeat-0.9.18/vl/mathop.h"
#include "../../common/vlfeat-0.9.18/vl/sift.h"
}

#include "sift_extractor.h"

using namespace cv;
using namespace std;

bool sift_keypoints_and_descriptors(string image_path, bool divide_512,	int verbose,
		bool display_image, vector<float*>& frames, vector<float*>& descr, unsigned int& number_desc) {

	struct timespec start, finish; //for measuring execution time
	double elapsed;

	if (verbose) {
		clock_gettime(CLOCK_MONOTONIC, &start);
	}

	// Loading image
	Mat image = imread(image_path.c_str(), CV_LOAD_IMAGE_GRAYSCALE);   // Read the file

	if (! image.data )                              // Check for invalid input
	{
		cout <<  "Could not open or find image " << image_path << endl ;
		return false;
	}

	int im_width = image.cols;
	int im_height = image.rows;

	// This is for debugging, checking the image was correctly loaded.
	if (display_image) {
		string window_name = "Image " + image_path;
		namedWindow( window_name, WINDOW_AUTOSIZE );// Create a window for display.
		imshow( window_name, image );               // Show our image inside it.
		waitKey(0);                                 // Wait for a keystroke in the window
	}

	// Transferring image to vlfeat structure
	unsigned int number_pixels = im_width*im_height;
	vl_sift_pix* data = new vl_sift_pix[number_pixels*sizeof(vl_sift_pix)];
	for (unsigned int ind = 0; ind < number_pixels; ind++) {
		data[ind] = static_cast<vl_sift_pix>(image.data[ind]);
	}

	// VLSIFT parameters
	int                O     = - 1 ;
	int                S     =   3 ;
	int                o_min =   -1 ;
	double             edge_thresh =  5;
	double             peak_thresh = 10 ;
	double             norm_thresh = -1 ;
	double             magnif      = -1 ;
	double             window_size = -1 ;

	bool            force_orientations = false ;

	VlSiftFilt* filt = vl_sift_new(im_width, im_height, O, S, o_min);

	int                nframes = 0, i,j,q ;

	if (peak_thresh >= 0) vl_sift_set_peak_thresh (filt, peak_thresh) ;
	if (edge_thresh >= 0) vl_sift_set_edge_thresh (filt, edge_thresh) ;
	if (norm_thresh >= 0) vl_sift_set_norm_thresh (filt, norm_thresh) ;
	if (magnif      >= 0) vl_sift_set_magnif      (filt, magnif) ;
	if (window_size >= 0) vl_sift_set_window_size (filt, window_size) ;

	if (verbose) {
	  printf("vl_sift: filter settings:\n") ;
	  printf("vl_sift:   image width           = %d\n",
				im_width) ;
	  printf("vl_sift:   image height          = %d\n",
				im_height) ;
	  printf("vl_sift:   octaves      (O)      = %d\n",
				vl_sift_get_noctaves      (filt)) ;
	  printf("vl_sift:   levels       (S)      = %d\n",
				vl_sift_get_nlevels       (filt)) ;
	  printf("vl_sift:   first octave (o_min)  = %d\n",
				vl_sift_get_octave_first  (filt)) ;
	  printf("vl_sift:   edge thresh           = %g\n",
				vl_sift_get_edge_thresh   (filt)) ;
	  printf("vl_sift:   peak thresh           = %g\n",
				vl_sift_get_peak_thresh   (filt)) ;
	  printf("vl_sift:   norm thresh           = %g\n",
				vl_sift_get_norm_thresh   (filt)) ;
	  printf("vl_sift:   window size           = %g\n",
				vl_sift_get_window_size   (filt)) ;

	  printf("vl_sift: will force orientations? %s\n",
				force_orientations ? "yes" : "no") ;
	}

	/* ...............................................................
	 *                                             Process each octave
	 * ............................................................ */

	i     = 0 ;
	bool first = true;
	while (true) {
		int                   err ;
		VlSiftKeypoint const *keys  = 0 ;
		int                   nkeys = 0 ;

		if (verbose) {
			printf ("vl_sift: processing octave %d\n",
					   vl_sift_get_octave_index (filt)) ;
		}

		/* Calculate the GSS for the next octave .................... */
		if (first) {
			err   = vl_sift_process_first_octave (filt, data) ;
			first = false;
		} else {
			err   = vl_sift_process_next_octave  (filt) ;
		}

		if (err) break ;

		if (verbose > 1) {
			printf("vl_sift: GSS octave %d computed\n",
					  vl_sift_get_octave_index (filt));
		}

		/* Run detector ............................................. */

		vl_sift_detect (filt) ;

		keys  = vl_sift_get_keypoints  (filt) ;
		nkeys = vl_sift_get_nkeypoints (filt) ;
		i     = 0 ;

		if (verbose > 1) {
		  printf ("vl_sift: detected %d (unoriented) keypoints\n", nkeys) ;
		}


		/* For each keypoint ........................................ */
		for (; i < nkeys ; ++i) {
			double                angles [4] ;
			int                   nangles ;
			VlSiftKeypoint const *k ;

			/* Obtain keypoint orientations ........................... */
			k = keys + i ;
			nangles = vl_sift_calc_keypoint_orientations(filt, angles, k) ;

			/* For each orientation ................................... */
			for (q = 0 ; q < nangles ; ++q) {
				vl_sift_pix rbuf [128] ;
				float* this_frame = new float[4*sizeof(float)];
				float* this_descr = new float[128*sizeof(float)];

				/* compute descriptor */
				vl_sift_calc_keypoint_descriptor (filt, rbuf, k, angles [q]) ;

				this_frame [0] = k -> x ;
				this_frame [1] = k -> y ;
				this_frame [2] = k -> sigma ;
				this_frame [3] = angles [q];

				frames.push_back(this_frame);

				for (j = 0 ; j < 128 ; ++j) {
					float x;
					if (divide_512) {
						x = rbuf [j] ;
					} else {
						x = 512.0F * rbuf [j] ;
					}
					this_descr [j] = x ;
				}
				descr.push_back(this_descr);

				++ nframes ;
			} /* next orientation */
		} /* next keypoint */
	} /* next octave */

	// Clean up
	/* release filter */
	if (filt) {
		vl_sift_delete(filt);
		filt = 0;
	}
	/* release image data */
	if (data) {
	  delete[] data;
	  data = 0 ;
	}

	// Write results
	if (verbose) {
		printf ("vl_sift: found %d keypoints\n", nframes) ;
		clock_gettime(CLOCK_MONOTONIC, &finish);
		elapsed = (finish.tv_sec - start.tv_sec);
		elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
		printf ("vl_sift: took %f seconds\n", elapsed) ;
	}

	number_desc = nframes;

	return true;
}

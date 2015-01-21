/*
 * sift_extractor.h
 *
 *  Created on: May 22, 2014
 *      Author: andrefaraujo
 */
#ifndef SIFT_EXTRACTOR_H
#define SIFT_EXTRACTOR_H

#include <iostream>
#include <vector>

using namespace std;

/*
 The function sift_keypoints_and_descriptors will get keypoints and descriptors
 from image data. It returns true if it works successfully.
 This is based on VLFEAT's implementation. It has been checked against VLFEAT
 calculations. Minor differences exist, due to different RGB to gray conversions
 between the opencv imread function and MATLAB's rgb2gray.
 Note: in this function, we fixed peak-thresh = 10 and edge-thresh = 5.
 INPUTS:
 - string image_path : path of the image which we will detect keypoints from.
 - bool divide_512 : descriptor elements are between 0 and 1 if this is set.
 - int verbose : verbose mode if verbose > 0
 - bool display_image : if true, grayscale image will be displayed. This is used
	 just to double-check image loading.
 - vector<float*>& frames : x, y, s and o of the keypoints, in this order, for each keypoint
	 (each element in the vector is a keypoint).
 - vector<float*>& descr : 128 dimensions of SIFT components, for each keypoint
 	 (each element in the vector is a decriptor).
 - unsigned int& number_desc : number of keypoints/descriptors that are detected in image.
 OUTPUT:
 It returns true if it works, and false otherwise.
 */
bool sift_keypoints_and_descriptors(string image_path, bool divide_512,	int verbose,
		bool display_image, vector<float*>& frames, vector<float*>& descr, unsigned int& number_desc);

#endif

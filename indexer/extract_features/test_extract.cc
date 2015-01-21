/*
 *
 * This program extracts features from one image.
 * This tests the sift_extractor library
 *
 *  Created on: May 22, 2014
 *      Author: andrefaraujo
 */

#include <iostream>

extern "C" {
#include "../../common/vlfeat-0.9.18/vl/sift.h"
}

#include "sift_extractor.h"
#include "../../common/feature_set/feature_set.h"

using namespace std;

const float DESC_DIVISOR = 512.0;

int main(int argc, const char * argv[]) {

	string image_path = "test_image.jpg";
	vector<float*> keypoints;
	vector<float*> descriptors;
	unsigned int number_desc;
	int verbose = 1;
	bool display_image = false;
	bool divide_512 = false;
	bool err =
		sift_keypoints_and_descriptors(image_path, divide_512, verbose,
					       display_image, keypoints, descriptors, number_desc);

	if (err) {
		cout << number_desc << " descriptors were detected successfully" << endl;

	}
	else {
		cout << "Descriptors were NOT detected successfully" << endl;
	}

	unsigned int desc_length = 128;
	unsigned int keypoint_length = 4;

	cout << "Adding to feature set!" << endl;

	FeatureSet* fSet = new FeatureSet(desc_length, keypoint_length);
	for (unsigned int count_desc = 0; count_desc < number_desc; count_desc++) {

		// Allocate new descriptor
		float* new_descriptor = new float[desc_length];
		for (unsigned int count_dim = 0; count_dim < desc_length; count_dim++) {
			// Division here is necessary to compensate for multiplier in vl_sift
			// (this division is already integrated into the reading function in 
			// common/file_io/file_io.cc)
			new_descriptor[count_dim] = descriptors.at(count_desc)[count_dim]
				/DESC_DIVISOR;
		}

		// Allocate new keypoint
		float* new_keypoint = new float[keypoint_length];
		for (unsigned int count_key = 0; count_key < keypoint_length; count_key++) {
			new_keypoint[count_key] = keypoints.at(count_desc)[count_key];
		}

		// Adding to feature set
		fSet->addFeature(new_descriptor, new_keypoint);

	}
	cout << "Added to feature set!" << endl;
	assert(number_desc == fSet->m_nNumFeatures);

	fSet->print();

	// Clean up
	for (unsigned int count_desc = 0; count_desc < number_desc; count_desc++) {
		if (keypoints.at(count_desc)) {
			delete[] keypoints.at(count_desc);
			keypoints.at(count_desc) = 0;
		}
		if (descriptors.at(count_desc)) {
			delete[] descriptors.at(count_desc);
			descriptors.at(count_desc) = 0;
		}
	}
	if (fSet) {
		delete fSet;
	}

	return EXIT_SUCCESS;
}

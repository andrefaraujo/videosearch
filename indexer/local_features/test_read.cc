/*
 *
 * This program illustrates a way to read the SIFT descriptors
 * that are saved by the extract program.
 *
 *  Created on: Jan 20, 2015
 *      Author: andrefaraujo
 */

#include <iostream>

#include "../../common/file_io/file_io.h"
#include "../../common/feature_set/feature_set.h"

using namespace std;

int main(int argc, const char * argv[]) {

	string descriptor_path = "test_image.siftb";

	uint sift_desc_dim = 128;
	uint sift_keypoint_dim = 4;

	FeatureSet* fSet = readSIFTFile(descriptor_path,
					sift_keypoint_dim,
					sift_desc_dim);

	fSet->print();

	// Clean up
	if (fSet) {
		delete fSet;
	}

	return EXIT_SUCCESS;
}

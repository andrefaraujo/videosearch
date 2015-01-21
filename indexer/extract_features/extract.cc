/*
 * main.cc
 *
 * This program extracts features from a list of images.
 *
 *  Created on: June 2014
 *      Author: andrefaraujo
 */

#include <iostream>
#include <fstream>
#include <thread>
#include <sys/stat.h>

#include "sift_extractor.h"
#include "../../common/feature_set/feature_set.h"

extern "C" {
#include "../../common/vlfeat-0.9.18/vl/sift.h"
}

using namespace std;

// Descriptor modes
const int SIFT = 0;

// SIFT constants
const unsigned int SIFT_DESCRIPTOR_LENGTH = 128;
const unsigned int SIFT_KEYPOINT_LENGTH = 4;

// Progress update constant
const unsigned int PROGRESS_INTERVAL = 100;

void save_sift_descriptor(vector<float*> keypoints, vector<float*> descriptors, unsigned int number_desc,
		string desc_file_name) {
	// Initialize file
	ofstream out_file;
	out_file.open(desc_file_name.c_str(), ios::out | ios::binary);
	if (!out_file.is_open()) {
		cout << "File " << desc_file_name << " could not be open." << endl;
		exit(EXIT_FAILURE);
	}


	// Save number of descs.
	int number_desc_int = static_cast<int>(number_desc); // reading format in REVV file is int
	out_file.write((char*)&number_desc_int, sizeof(int));

	// Save each desc. and keypoint
	for (unsigned int count_d = 0; count_d < number_desc; count_d++) {
		out_file.write((char*)keypoints.at(count_d),
				SIFT_KEYPOINT_LENGTH*sizeof(float)); //writes frame
		out_file.write((char*)descriptors.at(count_d),
				SIFT_DESCRIPTOR_LENGTH*sizeof(float)); //writes desc.
	}

	out_file.close();
}

void process_images(unsigned int start_index, unsigned int number_images, vector<string> all_images,
		int descriptor_mode, unsigned int thread_ind, bool skip_existing_files) {
	unsigned int count_proc_thread = 0;
	struct timespec start, finish; //for measuring execution time
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC, &start);

	for (unsigned int count_image = start_index; count_image < start_index + number_images;
			count_image++) {
		// Progress update
		if ((count_proc_thread++ % PROGRESS_INTERVAL) == 0) {
			cout << "Thread " << thread_ind << " is at image " << count_image
					<< ". First was " << start_index << " and last will be "
					<< start_index + number_images - 1 << endl;
			clock_gettime(CLOCK_MONOTONIC, &finish);
			elapsed = (finish.tv_sec - start.tv_sec);
			elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
			cout << "------> Last " << PROGRESS_INTERVAL << " took " << elapsed << " secs." << endl;
			clock_gettime(CLOCK_MONOTONIC, &start);
		}

		//Extract descriptors and keypoints and save them
		if (descriptor_mode == SIFT) {
			int lastindex = all_images.at(count_image).find_last_of(".");
			string desc_file_name = all_images.at(count_image).substr(0, lastindex) + ".siftb";

            // Check if file already exists
            struct stat stat_buffer;
            if (skip_existing_files && stat (desc_file_name.c_str(), &stat_buffer) == 0)
                continue;


			// Extracting
			vector<float*> keypoints;
			vector<float*> descriptors;
			unsigned int number_desc;
			// We will not divide by 512, since REVV code assumes it was multiplied by 512
			if (!sift_keypoints_and_descriptors(all_images.at(count_image), false, 0,
					false, keypoints, descriptors, number_desc)) {
				cout << "Descriptor extraction of " << all_images.at(count_image) << " did not work..." << endl;
				exit(EXIT_FAILURE);
			}

			// Saving
			save_sift_descriptor(keypoints, descriptors, number_desc, desc_file_name);

			// Clean up
			for (unsigned int count_d = 0; count_d < number_desc; count_d++) {
				if (keypoints.at(count_d)) {
					delete[] keypoints.at(count_d);
					keypoints.at(count_d) = 0;
				}
				if (descriptors.at(count_d)) {
					delete[] descriptors.at(count_d);
					descriptors.at(count_d) = 0;
				}
			}
		} else {
			cout << "descriptor_mode = " << descriptor_mode << " is not recognized" << endl;
			exit(EXIT_FAILURE);
		}

	}
}

void usage() {
	cout << "Extracts descriptors. Usage: "<< endl;
	cout << "./extract [options] -i image_list.txt "<< endl;
	cout << "image_list.txt is a file with one image path per line." << endl;
	cout << "We know it works with .jpg files; other extensions have currently not been tested." << endl;
	cout << "Options: " << endl;
	cout << "-h: print this help" << endl;
	cout << "-t X: number of threads to be used = X (default: 1)" << endl;
	cout << "-d X: descriptor mode = X (SIFT = 0, default; other modes currently not supported)" << endl;
	cout << "-s: skip an image if a feature file already exists for it." << endl;
}

int main(int argc, const char * argv[]) {

	string image_list_path;
	unsigned int number_threads = 1;
	int descriptor_mode = SIFT;
	bool skip = false;

	// Parse command-line
	for (unsigned int i = 1; i < argc; i++) {
	    if (strcmp(argv[i],"-h")==0) {
	    	usage();
	    	exit(EXIT_FAILURE);
	    } else if (strcmp(argv[i],"-i")==0) {
	        if (i + 2 > argc) {
	        	cout << "Command-line error with option -i" << endl;
	        	exit(EXIT_FAILURE);
	        } else {
	        	image_list_path = string(argv[i+1]);
	        }
	        i++;
	    } else if (strcmp(argv[i],"-t")==0) {
	        if (i + 2 > argc) {
	        	cout << "Command-line error with option -t" << endl;
	        	exit(EXIT_FAILURE);
	        } else {
	        	number_threads = static_cast<unsigned int>(atoi(argv[i + 1]));
	        }
	        i++;
	    } else if (strcmp(argv[i],"-d")==0) {
	        if (i + 2 > argc) {
	        	cout << "Command-line error with option -d" << endl;
	        	exit(EXIT_FAILURE);
	        } else {
	        	descriptor_mode = atoi(argv[i + 1]);
	        }
	        i++;
	    } else if (strcmp(argv[i],"-s")==0) {
		    skip = true;
	    }
	}
	if (image_list_path.empty()) {
		cout << "Wrong arguments!! Look at options below:" << endl;
		cout << endl;
		cout << endl;
		cout << endl;
		usage();
		exit(EXIT_FAILURE);
	}
	cout << "Input filename: " << image_list_path << endl;
	cout << "Number of threads: " << number_threads << endl;
	cout << "Descriptor mode: " << descriptor_mode << endl;
	cout << "Skip existing feature files: " << skip << endl;

	// Process list of files
	ifstream in_file(image_list_path.c_str());
	string line;
	vector<string> all_images;
	while (in_file >> line) {
		all_images.push_back(line);
	}
	in_file.close();

	unsigned int number_images = all_images.size();

	// Define number of images per thread
	unsigned int computations_per_thread = static_cast<unsigned int>(floor(number_images/number_threads));
	int rest_computations = number_images % number_threads;
	unsigned int number_threads_to_use = number_threads;

	if (computations_per_thread == 0) { //less work than thread..
		computations_per_thread = 1;
		rest_computations = 0;
		number_threads_to_use = number_images;
	}

	// Start threads
	vector<thread> extract_threads;
	unsigned int number_images_so_far = 0;
	for (unsigned int count_thread = 0; count_thread < number_threads_to_use; count_thread++) {
		unsigned int start_index = number_images_so_far;
		unsigned int number_images_per_thread;

		// Casting is necessary here since we're comparing numbers that might be negative:
		if (static_cast<int>(count_thread) > rest_computations - 1) {
			number_images_per_thread = computations_per_thread;
		}
		else {
			number_images_per_thread = computations_per_thread + 1;
		}

		extract_threads.push_back(thread(process_images, start_index, number_images_per_thread,
				all_images, descriptor_mode, count_thread, skip));

		number_images_so_far += number_images_per_thread;
	}

	for (unsigned int count_thread = 0; count_thread < number_threads_to_use; count_thread++) {
		extract_threads.at(count_thread).join();
	}

	return EXIT_SUCCESS;
}


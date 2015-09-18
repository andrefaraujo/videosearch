/**********************************************************

This program will detect shots from a sequence of frames.
It will return a text file with one number per line. This number
indicates the first frame in a shot.

**********************************************************/

#include <fstream>
#include <iostream>
#include <opencv2/opencv.hpp>

using namespace std;
using namespace cv;

const int MAX_FRAME_WIDTH = 480;

void usage() {
	cout << "Detect shots from a list of video frames" << endl;
	cout << "Usage:" << endl;
	cout << "./shot_detector [options] --list[-l] database_list_file --output[-o] output_file_path" << endl;
	cout << "Command-line example (be sure to extract keyframes beforehand for this example to work):" << endl;
    cout << "./shot_detector -l ../test_db/test_video_0_keyframes.txt -t 0.8 -o ../test_db/test_video_0_keyframes.shot_t0.8 -v 0"<< endl;
	cout << "The program will output a text file, containing one number per line" << endl;
	cout << "The number corresponds to the first frame number of the shot (starting at 0)" << endl;
	cout << "Options:" << endl;
	cout << "--verbose_level[-v] ARG: (default: 1)" << endl;
	cout << "--threshold[-t] ARG: a number between 0 and 1; the higher, the less strict (more shots are found), and vice-versa (default: 0.7)" << endl;
}

bool shot_bound(Mat& curr_frame, Mat& prev_frame, int verbose_level, double shot_detector_thresh);

int main(int argc, char* * argv) {
	// Mandatory arguments
	string db_list_path = "";
	string output_file_path = "";
	
	// Default values for options
	int verbose_level = 1;
	double shot_detector_thresh = 0.7;

	if (argc < 3) {
		cout << "Wrong usage!!!" << endl;
		cout << "***********************************" << endl;
		cout << "***********************************" << endl;
		cout << "***********************************" << endl;
		cout << "See usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	} else {
		for (int count_arg = 1; count_arg < argc; count_arg++) {
			if ((!strcmp(argv[count_arg], "--list")) || (!strcmp(argv[count_arg], "-l"))) {
				db_list_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--output")) || (!strcmp(argv[count_arg], "-o"))) {
				output_file_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--verbose_level")) || (!strcmp(argv[count_arg], "-v"))) {
				verbose_level = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--threshold")) || (!strcmp(argv[count_arg], "-t"))) {
				shot_detector_thresh = atof(argv[count_arg + 1]);
				count_arg++;
			}			
		}
	}

	// Check that all needed arguments were provided
	if (db_list_path == "") {
		cout << "database frames list argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}
	if (output_file_path == "") {
		cout << "output argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}

	if (verbose_level) {
		cout << "Starting shot detection using:" << endl;
		cout << "------>db_list_path = " << db_list_path  << endl;
		cout << "------>output_file_path = " << output_file_path  << endl;
		cout << "------>verbose_level = " << verbose_level  << endl;
		cout << "------>shot_detector_thresh = " << shot_detector_thresh  << endl;
	}

	ifstream list_file(db_list_path.c_str());
	ofstream out_file;
	out_file.open(output_file_path.c_str());
	string list_line;
	unsigned int frame_number = 0;
	Mat prev_frame;
	while (list_file >> list_line) {
		// Loading image
		if (verbose_level >= 2) cout << "Doing frame " << frame_number << endl;
		if (verbose_level >=3) cout << "Loading image " << list_line << endl;
		Mat curr_frame = imread(list_line.c_str(), CV_LOAD_IMAGE_COLOR);   // Read the file
		if (! curr_frame.data )                              // Check for invalid input
		{
			cout <<  "Could not open or find image " << list_line << endl ;
			exit(EXIT_FAILURE);
		}
		
		if (frame_number) {
			if (shot_bound(curr_frame, prev_frame, verbose_level, shot_detector_thresh)) {
				out_file << frame_number << endl;
				if (verbose_level) cout << "Found a shot starting at frame " << frame_number << endl;
			}
		} else {
			// First frame, so it will be the beginning of a shot
			out_file << frame_number << endl;
			if (verbose_level) cout << "Found a shot starting at frame " << frame_number << endl;
		}
		curr_frame.copyTo(prev_frame);
		frame_number++;
	}

	out_file.close();

	return EXIT_SUCCESS;
}

bool shot_bound(Mat& curr_frame, Mat& prev_frame, int verbose_level, double shot_detector_thresh) {
	// Resize if necessary	
	if (curr_frame.cols > MAX_FRAME_WIDTH) {
		resize(curr_frame, curr_frame, Size(static_cast<float>(MAX_FRAME_WIDTH), 
						    round((static_cast<float>(MAX_FRAME_WIDTH)/static_cast<float>(curr_frame.cols))*static_cast<float>(curr_frame.rows))));
	}
	if (prev_frame.cols > MAX_FRAME_WIDTH) {
		resize(prev_frame, prev_frame, Size(static_cast<float>(MAX_FRAME_WIDTH), 
						    round(static_cast<float>(MAX_FRAME_WIDTH)/static_cast<float>(prev_frame.cols))*static_cast<float>(prev_frame.rows)));
	}

	// Convert to HSV
	Mat hsv1, hsv2;
	cvtColor(curr_frame, hsv1, CV_BGR2HSV);
	cvtColor(prev_frame, hsv2, CV_BGR2HSV);

	// Making sure resizing and color conversion were done correctly
	assert(hsv1.data);
	assert(hsv2.data);

	if (verbose_level >= 3) {
		// Writing some values just to test
		unsigned char *hsv1_data_ptr = (unsigned char*)(hsv1.data);
		cout << "curr image" << endl;
		cout << "H at 0,0: " << hsv1_data_ptr[0] << endl;
		cout << "S at 0,0: " << hsv1_data_ptr[1] << endl;
		cout << "V at 0,0: " << hsv1_data_ptr[2] << endl;
		cout << "H at row 100, col 100: " << hsv1_data_ptr[100*hsv1.cols*hsv1.channels() + 100*hsv1.channels() + 0] << endl;
		cout << "S at row 100, col 100: " << hsv1_data_ptr[100*hsv1.cols*hsv1.channels() + 100*hsv1.channels() + 1] << endl;
		cout << "V at row 100, col 100: " << hsv1_data_ptr[100*hsv1.cols*hsv1.channels() + 100*hsv1.channels() + 2] << endl;

		cout << "prev image" << endl;
		cout << "H at 0,0: " << hsv2.data[0] << endl;
		cout << "S at 0,0: " << hsv2.data[1] << endl;
		cout << "V at 0,0: " << hsv2.data[2] << endl;
		cout << "H at row 100, col 100: " << hsv2.data[100*hsv2.cols*hsv2.channels() + 100*hsv2.channels() + 0] << endl;
		cout << "S at row 100, col 100: " << hsv2.data[100*hsv2.cols*hsv2.channels() + 100*hsv2.channels() + 1] << endl;
		cout << "V at row 100, col 100: " << hsv2.data[100*hsv2.cols*hsv2.channels() + 100*hsv2.channels() + 2] << endl;

	}

	// Compute HSV histograms
	int hbins = 16, sbins = 16, vbins = 16;
	int histSize[] = {hbins, sbins, vbins};
	// hue varies from 0 to 179, see cvtColor
	float hranges[] = { 0, 180 };
	// saturation and value vary from 0 to 255
	float sranges[] = { 0, 256 };
	float vranges[] = { 0, 256 };
	const float* ranges[] = { hranges, sranges, vranges };
	MatND hist1, hist2;
	int channels[] = {0, 1, 2};
	bool bUniform = true, bAccumulate = false;
	if (verbose_level >= 3) cout << "Calculating histogram for current frame..." << endl;
	calcHist( &hsv1, 1, channels, Mat(), hist1, 3, histSize, ranges, bUniform, bAccumulate);
	if (verbose_level >= 3) cout << "done!" << endl;

	if (verbose_level >= 3) cout << "Calculating histogram for previous frame..." << endl;
	calcHist( &hsv2, 1, channels, Mat(), hist2, 3, histSize, ranges, bUniform, bAccumulate);
	if (verbose_level >= 3) cout << "done!" << endl;

	if (verbose_level >= 3) cout << "Normalizing histogram for current frame..." << endl;
	normalize(hist1, hist1, 1, 0, NORM_L1, -1, Mat());
	if (verbose_level >= 3) cout << "done!" << endl;

	if (verbose_level >= 3) cout << "Normalizing histogram for previous frame..." << endl;
	normalize(hist2, hist2, 1, 0, NORM_L1, -1, Mat());
	if (verbose_level >= 3) cout << "done!" << endl;

	// Compare histograms
	double dCorr = compareHist(hist1, hist2, CV_COMP_INTERSECT);
	bool bShotBound = dCorr < shot_detector_thresh;

	if (verbose_level >=2 ) cout << "Hist comp. score = " << dCorr << endl;
	if (verbose_level >=2 ) cout << "Did we find a new shot? = " << bShotBound << endl;

	return bShotBound;

}

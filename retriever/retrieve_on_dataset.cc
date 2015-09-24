/**********************************************************
This program will generate retrieval results for a specific
dataset, with a specified set of parameters
**********************************************************/

#include <cassert>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "retriever.h"

using namespace std;

// Useful constants for file path construction
const string MP4_EXTENSION = ".mp4";
const string CENTROIDS_PRE_EXTENSION_FV = "_keyframes.sift_fv_idx_k";
const string CENTROIDS_PRE_EXTENSION_SCFV = "_keyframes.sift_scfv_idx_k";
const string SHOT_PRE_EXTENSION = "_shot_t";
const string FPSHOT_PRE_EXTENSION = "_n";
const string SHOT_MODE_PRE_EXTENSION = "_m";

int str2Int(string string);
string getSceneSigPath(string scene_path, uint number_centroids, 
		       float shot_threshold, int fpshot, uint shot_mode, bool use_fv);
vector<string> getVectorOfStringsFromFileLines(string file_name);
void parseAggList(string agg_list_path, int verbose_level, 
		  vector < vector < uint > >& agg_list_vec);
void parseShotInfo(string shot_helper_path, int verbose_level, uint number_centroids,
		   int fpshot, float shot_threshold, int shot_mode, bool use_fv,
		   vector < pair < string, pair < uint, uint > > >& shot_info);

void usage() {
	cout << "Perform retrieval using a specific dataset" << endl;
	cout << "Usage:" << endl;
	cout << "./retrieve_on_dataset [options] --index[-i] gdindex_file --db_list[-d] database_list_file --query_list[-q] query_list_file --output[-o] output_base_path" << endl;
	cout << "For example:" << endl;
    // TODO: example
	cout << "./retrieve_on_dataset -i YOUR_INDEX_FILE -d ../indexer/ -q ../indexer/extract_features/cnn12h_lists/all_cnn12h_queries.txt -g /home/andrefaraujo/datasets/VideoSearchImageQueries/Queries_12hCNN/groundTruth.txt -o results_cnn12h/test_system" << endl;
	cout << "The system will output 2 different files (one with log, one with results), using as the base of the name the output argument." << endl;
	cout << "Options:" << endl;
	cout << "--feature_mode[-f] ARG: feature (local descriptor) mode (default: 0 = SIFT)" << endl;
	cout << "--centroids[-c] ARG: number of centroids/Gaussians to use in global signatures (default: 512)" << endl;
	cout << "--keyframe_numbers[-e] ARG: path to file containing frame numbers to use. This is particularly useful when using some shot modes. (default: not using it)" << endl;
	cout << "--gdindex_parameters_path ARG: path where GDIndex pretrained parameters are saved (default: ../indexer/global_descriptors/trained_parameters)" << endl;
	cout << "--number_output_results[-s] ARG: size of list of results to output (default: 100)" << endl;
	cout << "--verbose_level[-v] ARG: (default: 1)" << endl;
	cout << "--shot_mode[-m] ARG: -1 = not using this; 0 = using indep. keyframes in a shot; 1 = aggregating local features from several frames into a global signature; 2 = aggregating global signatures from several frames into a global signature; 3 = aggregating local features from several frames, using tracking, into a global signature (default: -1)" << endl;
	cout << "--shot_list[-l] ARG: path of list of shots, used only if shot_mode == 1 (default: not using it)" << endl;
    cout << "--min_number_words_visited ARG: ARG = minimum number of centroids/Gaussians necessary to consider a given database item. Default: 20." << endl;
    cout << "--word_selection_mode ARG: Options: (ARG = 0 : use L1 norm mode), (ARG = 1 : use total soft assgn mode). Default: 0." << endl;
	cout << "--number_scenes_rerank ARG: when using two-stage retrieval, ARG is the number ofscenes to re-rank. Default: 0" << endl;
	cout << "--number_centroids_rerank ARG: when using two-stage retrieval, ARG is the number of centroids/Gaussians to use in the global signatures of the re-ranking stage. Default: 0" << endl;
	cout << "--group_lists_path ARG: ARG = path to file containing group configurations, when using two-stage retrieval. Default: empty" << endl;
	cout << "--shot_info_path ARG: ARG = path to file containing shot information, when using two-stage retrieval. Default: empty" << endl;
	cout << "--fpshot ARG: ARG = Number of frames per shot used in the shots when using two-stage retrieval. Default: 10" << endl;
	cout << "--shot_threshold ARG: ARG = Shot threshold used in the shots when using two-stage retrieval. Default: 0.8" << endl;
	cout << "--word_selection_thresh ARG: ARG = Threshold for asymmetric distance computation. Default: 7.0" << endl;
	cout << "--word_selection_thresh_second_stage ARG: ARG = Threshold for asymmetric distance computation, used in re-ranking stage if using two-stage retrieval. Default: 8.0" << endl;
	cout << "--gdindex_path_other ARG: ARG = path to index file for index used in second stage, used when using two-stage retrieval." << endl;
    cout << "--avoid_redundant_scene_results: flag that, if set, will make the program write more than number_output_results to the results file, such that it contains exactly number_output_results scenes (clips)"<< endl;
}

int main(int argc, char* * argv) {
	// Mandatory arguments
	string gdindex_path = "";
	string db_list_path = "";
	string query_list_path = "";
	string output_base_path = "";
	
	// Default values for options
    int feat_mode = 0;
	uint number_centroids = 512;
	string keyframe_numbers_path = "";
    string gdindex_parameters_path = "../indexer/global_descriptors/trained_parameters";
	uint number_output_results = 100;
	int verbose_level = 1;
	int shot_mode = -1;
	string shot_list_path = "";
    uint min_number_words_visited = 20;
    int word_selection_mode = 0;
	uint number_scenes_rerank = 0;
	uint number_centroids_rerank = 0;
	string group_lists_path = "";
	string shot_info_path = "";
	int fpshot = 10;
	float shot_threshold = 0.8;
	float word_selection_thresh = 7.0;
	float word_selection_thresh_second_stage = 8.0;
	string gdindex_path_other = "";
    bool avoid_redundant_scene_results = true;

	if (argc < 9) {
		cout << "Wrong usage!!!" << endl;
		cout << "***********************************" << endl;
		cout << "***********************************" << endl;
		cout << "***********************************" << endl;
		cout << "See usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	} else {
		for (int count_arg = 1; count_arg < argc; count_arg++) {
			if ((!strcmp(argv[count_arg], "--index")) || (!strcmp(argv[count_arg], "-i"))) {
				gdindex_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--db_list")) || (!strcmp(argv[count_arg], "-d"))) {
				db_list_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--query_list")) || (!strcmp(argv[count_arg], "-q"))) {
				query_list_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--output")) || (!strcmp(argv[count_arg], "-o"))) {
				output_base_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--feature_mode")) || (!strcmp(argv[count_arg], "-f"))) {
				feat_mode = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--number_output_results")) || (!strcmp(argv[count_arg], "-s"))) {
				number_output_results = static_cast<uint>(atoi(argv[count_arg + 1]));
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--centroids")) || (!strcmp(argv[count_arg], "-c"))) {
				number_centroids = static_cast<uint>(atoi(argv[count_arg + 1]));
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--verbose_level")) || (!strcmp(argv[count_arg], "-v"))) {
				verbose_level = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--gdindex_parameters_path")) {
				gdindex_parameters_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--keyframe_numbers")) || (!strcmp(argv[count_arg], "-e"))) {
				keyframe_numbers_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--shot_list")) || (!strcmp(argv[count_arg], "-l"))) {
				shot_list_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if ((!strcmp(argv[count_arg], "--shot_mode")) || (!strcmp(argv[count_arg], "-m"))) {
				shot_mode = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--word_selection_mode")) {
				word_selection_mode = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--min_number_words_visited")) {
				min_number_words_visited = static_cast<uint>(atoi(argv[count_arg + 1]));
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--number_scenes_rerank")) {
				number_scenes_rerank = static_cast<uint>(atoi(argv[count_arg + 1]));
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--number_centroids_rerank")) {
				number_centroids_rerank = static_cast<uint>(atoi(argv[count_arg + 1]));
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--fpshot")) {
				fpshot = atoi(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--shot_threshold")) {
				shot_threshold = atof(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--word_selection_thresh")) {
				word_selection_thresh = atof(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--word_selection_thresh_second_stage")) {
				word_selection_thresh_second_stage = atof(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--group_lists_path")) {
				group_lists_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--shot_info_path")) {
				shot_info_path = string(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--gdindex_path_other")) {
				gdindex_path_other = string(argv[count_arg + 1]);
				count_arg++;
			} else if (!strcmp(argv[count_arg], "--avoid_redundant_scene_results")) {
				avoid_redundant_scene_results = true;
			} else {
				cout << "Unrecognized argument " << argv[count_arg] 
				     << " , quitting..." << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	// Check that all needed arguments were provided
	if (gdindex_path == "") {
		cout << "index argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}
	if (db_list_path == "") {
		cout << "database list argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}
	if (query_list_path == "") {
		cout << "query list argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}
	if (output_base_path == "") {
		cout << "output argument was not provided... see usage below:" << endl;
		usage();
		exit(EXIT_FAILURE);
	}

	if (verbose_level) {
		cout << "Starting evaluation using:" << endl;
		cout << "------>gdindex_path = " << gdindex_path  << endl;
		cout << "------>db_list_path = " << db_list_path  << endl;
		cout << "------>query_list_path = " << query_list_path  << endl;
		cout << "------>output_base_path = " << output_base_path  << endl;
		cout << "------>feat_mode = " << feat_mode  << endl;
		cout << "------>number_centroids = " << number_centroids << endl;
		cout << "------>keyframes_numbers_path = " << keyframe_numbers_path << endl;
		cout << "------>gdindex_parameters_path = " << gdindex_parameters_path << endl;
		cout << "------>number_output_results = " << number_output_results  << endl;
		cout << "------>verbose_level = " << verbose_level  << endl;
		cout << "------>shot_mode = " << shot_mode << endl;
		cout << "------>shot_list_path = " << shot_list_path << endl;
		cout << "------>min_number_words_visited = " << min_number_words_visited << endl;
		cout << "------>word_selection_mode = " << word_selection_mode << endl;
		cout << "------>number_scenes_rerank = " << number_scenes_rerank << endl;
		cout << "------>number_centroids_rerank = " << number_centroids_rerank << endl;
		cout << "------>group_lists_path = " << group_lists_path << endl;
		cout << "------>shot_info_path = " << shot_info_path << endl;
		cout << "------>fpshot = " << fpshot << endl;
		cout << "------>shot_threshold = " << shot_threshold << endl;
		cout << "------>word_selection_thresh = " << word_selection_thresh << endl;
		cout << "------>word_selection_thresh_second_stage = " << word_selection_thresh_second_stage << endl;
		cout << "------>gdindex_path_other = " << gdindex_path_other << endl;
        cout << "------>avoid_redundant_scene_results = " << avoid_redundant_scene_results << endl;
	}

	// Instantiate retriever
	Retriever r;

	// Setting parameters...
	r.set_verbose_level(verbose_level);
	r.set_gdindex_path(gdindex_parameters_path);
	r.set_local_descriptor_mode(feat_mode);
	r.set_number_output_results(number_output_results);
	r.set_number_gaussians_global_descriptor(number_centroids);
    r.set_word_selection_mode(word_selection_mode);
    r.set_word_selection_thresh(word_selection_thresh);
    r.set_min_num_words_visited(min_number_words_visited);

	// Starting retrieval
    vector < vector < uint > > group_lists;
    vector < pair < string, pair < uint, uint > > > shot_info;
    if (group_lists_path != "") {
        // Parse group lists
        parseAggList(group_lists_path, 0, group_lists);			
    }
    if (shot_info_path != "") {
        // Parse shot info
        parseShotInfo(shot_info_path, 0, number_centroids_rerank,
				      fpshot, shot_threshold, shot_mode, false, shot_info);
    }
    r.retrieve_on_specific_dataset(gdindex_path, 
                                   db_list_path, 
                                   query_list_path,
                                   output_base_path,
                                   keyframe_numbers_path, 
                                   shot_list_path, 
                                   shot_mode, 
                                   number_scenes_rerank, 
                                   number_centroids_rerank, 
                                   group_lists, 
                                   shot_info, 
                                   word_selection_thresh_second_stage,
                                   gdindex_path_other, 
                                   avoid_redundant_scene_results);
		
	return EXIT_SUCCESS;
}

void parseShotInfo(string shot_helper_path, int verbose_level, uint number_centroids,
		   int fpshot, float shot_threshold, int shot_mode, bool use_fv,
		   vector < pair < string, pair < uint, uint > > >& shot_info) {
	// Get all lines
	vector <string> file_lines = getVectorOfStringsFromFileLines(shot_helper_path);
	uint number_lines = file_lines.size();
	if (verbose_level >=3) cout << "Read " << number_lines << " from " 
				    << shot_helper_path << endl;

	// Loop over all lines and parse things into shot_info
	for (uint count_line = 0; count_line < number_lines; count_line++) {
		if (verbose_level >= 3) cout << "count_line = " << count_line << endl;
		// Line under consideration
		string this_line = file_lines.at(count_line);
		if (verbose_level >= 3) cout << "Considering line " << this_line << endl;

		// Split string
		vector<string> line_elements;
		stringstream ss(this_line);
		string token;
		while(getline(ss, token, ' ')) {
			line_elements.push_back(token);
		}
		
		// Get different elements in that line
		uint shot_number = static_cast<uint>(str2Int(line_elements.at(0)));
		assert(shot_number == count_line);
		uint start_shot = static_cast<uint>(str2Int(line_elements.at(2)));
		uint end_shot = static_cast<uint>(str2Int(line_elements.at(3)));
		if (verbose_level >= 3) cout << "Read start_shot =  " << start_shot 
					     << " and end_shot = " << end_shot << endl;
		assert(start_shot <= end_shot);
		pair < uint, uint > pair_shot_numbers;
		pair_shot_numbers.first = start_shot;
		pair_shot_numbers.second = end_shot;
		string scene_path = getSceneSigPath(line_elements.at(1), number_centroids, 
						    shot_threshold, fpshot, shot_mode, use_fv);
		if (verbose_level >=3) cout << "Parsed scene_index_path = " << scene_path << endl;
		pair < string, pair < uint, uint> > pair_shot_info;
		pair_shot_info.first = scene_path;
		pair_shot_info.second = pair_shot_numbers;
		shot_info.push_back(pair_shot_info);
	}
}

void parseAggList(string agg_list_path, int verbose_level, 
		  vector < vector < uint > >& agg_list_vec) {
	// Get all lines
	vector <string> file_lines = getVectorOfStringsFromFileLines(agg_list_path);
	uint number_lines = file_lines.size();
	if (verbose_level >=3) cout << "Read " << number_lines << " from " 
				    << agg_list_path << endl;
	
	// Loop over all lines and parse things into agg_list_vec
	for (uint count_line = 0; count_line < number_lines; count_line++) {
		if (verbose_level >= 3) cout << "count_line = " << count_line << endl;
		// Line under consideration
		string this_line = file_lines.at(count_line);
		if (verbose_level >= 3) cout << "Considering line " << this_line << endl;

		// Split string
		vector<string> line_elements;
		stringstream ss(this_line);
		string token;
		while(getline(ss, token, ' ')) {
			line_elements.push_back(token);
		}
		
		// Get different elements in that line
		uint number_elements = line_elements.size();
		vector < uint > this_vec;
		if (verbose_level >= 3) cout << "Starting to add elements..." << endl;
		for (uint count_el = 0; count_el < number_elements; count_el++) {			       
			uint item_to_add = static_cast<uint>(str2Int(line_elements.at(count_el)));
			this_vec.push_back(item_to_add);
			if (verbose_level >= 3) cout << "Added  " << item_to_add << endl;
		}
		agg_list_vec.push_back(this_vec);
	}
	
	
}
vector<string> getVectorOfStringsFromFileLines(string file_name) {

	ifstream in_file(file_name.c_str());
	string line;
	vector<string> output;
	while (!in_file.eof()) {
		if (getline(in_file, line)) output.push_back(line);
	}
	return output;
}
int str2Int(string string){
    stringstream parser(string);
    int i;
    parser>>i;
    return i;
}

string getSceneSigPath(string scene_path, uint number_centroids, 
		       float shot_threshold, int fpshot, uint shot_mode, bool use_fv) {
	// Replace extension by shot extension
	uint pos = scene_path.find(MP4_EXTENSION);
	assert(pos != string::npos); // scene_path should contain MP4_EXTENSION
	stringstream ss(scene_path.substr(0, pos));
	string base_path;
	ss >> base_path;
	ostringstream ss_thresh;
	ss_thresh << setprecision(1) << shot_threshold;
	string shot_str = ss_thresh.str();
	ostringstream ss_centroids;
	ss_centroids << number_centroids;
	string centroids_str = ss_centroids.str();
	ostringstream ss_fpshot;
	ss_fpshot << fpshot;
	string fpshot_str = ss_fpshot.str();
	ostringstream ss_shot_mode;
	ss_shot_mode << shot_mode;
	string shot_mode_str = ss_shot_mode.str();
	string shot_path;
	if (use_fv) {
		shot_path = base_path + CENTROIDS_PRE_EXTENSION_FV + centroids_str 
			+ SHOT_PRE_EXTENSION + shot_str
			+ FPSHOT_PRE_EXTENSION + fpshot_str
			+ SHOT_MODE_PRE_EXTENSION + shot_mode_str;
	} else {
		shot_path = base_path + CENTROIDS_PRE_EXTENSION_SCFV + centroids_str 
			+ SHOT_PRE_EXTENSION + shot_str
			+ FPSHOT_PRE_EXTENSION + fpshot_str
			+ SHOT_MODE_PRE_EXTENSION + shot_mode_str;
	}
	return shot_path;		
}

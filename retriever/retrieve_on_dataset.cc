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

int str2int(const string string);
void get_vector_of_strings_from_file_lines(const string file_name,
                                           vector<string>& out);
void parse_agg_list(const string agg_list_path, const int verbose_level, 
                    vector < vector < uint > >& agg_list_vec);

void usage() {
    cout << "Perform retrieval using a specific dataset" << endl;
    cout << "Usage:" << endl;
    cout << "./retrieve_on_dataset [options] --index[-i] gdindex_file --db_list[-d] database_list_file --query_list[-q] query_list_file --output[-o] output_base_path" << endl;
    cout << "For example:" << endl;
    cout << "./retrieve_on_dataset -i YOUR_INDEX_FILE -d ../indexer/ -q ../indexer/extract_features/cnn12h_lists/all_cnn12h_queries.txt -g /home/andrefaraujo/datasets/VideoSearchImageQueries/Queries_12hCNN/groundTruth.txt -o results_cnn12h/test_system" << endl;
    cout << "The system will output 2 different files (one with log, one with results), using as the base of the name the output argument." << endl;
    cout << "Options:" << endl;
    cout << "--query_index ARG: index of global descriptors for queries (default: not using this, and global descriptors are computed on the fly from local descriptors)" << endl;
    cout << "--feature_mode[-f] ARG: feature (local descriptor) mode (default: 0 = SIFT); another option is 1 = SIFTGEO" << endl;
    cout << "--gd_intra_normalization: Boolean that sets usage of intra-normalization mode for FVs (default: false)" << endl;
    cout << "--gd_unbinarized: Boolean that sets usage of standard FVs, without binarization (default: false)" << endl;
    cout << "--centroids[-c] ARG: number of centroids/Gaussians to use in global signatures (default: 512)" << endl;
    cout << "--keyframe_numbers[-e] ARG: path to file containing frame numbers to use. This is particularly useful when using some shot modes. (default: not using it)" << endl;
    cout << "--gdindex_parameters_path ARG: path where GDIndex pretrained parameters are saved (default: ../indexer/global_descriptors/trained_parameters)" << endl;
    cout << "--number_output_results[-s] ARG: size of list of results to output (default: 100)" << endl;
    cout << "--verbose_level[-v] ARG: (default: 1)" << endl;
    cout << "--shot_mode[-m] ARG: -1 = not using this; 0 = using indep. keyframes in a shot; 1 = aggregating local features from several frames into a global signature; 2 = aggregating global signatures from several frames into a global signature; 3 = aggregating local features from several frames, using tracking, into a global signature (default: -1)" << endl;
    cout << "--shot_list[-l] ARG: path of list of shots, used only if shot_mode == 1, 2 or 3 (default: not using it)" << endl;
    cout << "--min_number_words_visited ARG: ARG = minimum number of centroids/Gaussians necessary to consider a given database item. Default: 0." << endl;
    cout << "--word_selection_mode ARG: Options: (ARG = 0 : use L1 norm mode), (ARG = 1 : use total soft assgn mode). Default: 0." << endl;
    cout << "--number_scenes_rerank ARG: when using two-stage retrieval, ARG is the number ofscenes to re-rank. Default: 0" << endl;
    cout << "--number_centroids_rerank ARG: when using two-stage retrieval, ARG is the number of centroids/Gaussians to use in the global signatures of the re-ranking stage. Default: 0" << endl;
    cout << "--group_lists_rerank_path ARG: ARG = path to file containing group configurations, when using two-stage retrieval. Default: empty" << endl;
    cout << "--asym_scoring_mode ARG: ARG = asymmetric scoring mode. Options: (0 = no asymmetric scoring), (1 = QAGS), (2 = DAGS), (3 = SGS). Default: 1 (QAGS)" << endl;
    cout << "--asym_scoring_mode_rerank ARG: ARG = asymmetric scoring mode for re-ranking, when using two-stage retrieval. Options: (0 = no asymmetric scoring), (1 = QAGS), (2 = DAGS), (3 = SGS). Default: 1 (QAGS)" << endl;
    cout << "--score_den_power_norm ARG: ARG = power normalization for score denominator. Default: 0.5" << endl;
    cout << "--score_den_power_norm_rerank ARG: ARG = power normalization for score denominator in re-ranking stage, when using two-stage retrieval. Default: 0.5" << endl;
    cout << "--word_selection_thresh ARG: ARG = Threshold for asymmetric distance computation. Default: 7.0" << endl;
    cout << "--word_selection_thresh_rerank ARG: ARG = Threshold for asymmetric distance computation, used in re-ranking stage if using two-stage retrieval. Default: 8.0" << endl;
    cout << "--gdindex_path_rerank ARG: ARG = path to index file for index used in second stage, used when using two-stage retrieval." << endl;
    cout << "--avoid_redundant_scene_results: flag that, if set, will make the program write more than number_output_results to the results file, such that it contains exactly number_output_results scenes (clips)"<< endl;
}

int main(int argc, char* * argv) {
    // Mandatory arguments
    string gdindex_path = "";
    string db_list_path = "";
    string query_list_path = "";
    string output_base_path = "";
    
    // Default values for options
    string query_index_path = "";
    int feat_mode = 0;
    uint number_centroids = 512;
    string keyframe_numbers_path = "";
    string gdindex_parameters_path = "../indexer/global_descriptors/trained_parameters";
    uint number_output_results = 100;
    int verbose_level = 1;
    int shot_mode = -1;
    string shot_list_path = "";
    uint min_number_words_visited = 0;
    int word_selection_mode = 0;
    uint number_scenes_rerank = 0;
    uint number_centroids_rerank = 0;
    string group_lists_rerank_path = "";
    int asym_scoring_mode = 1;
    int asym_scoring_mode_rerank = 1;
    float word_selection_thresh = 7.0;
    float word_selection_thresh_rerank = 8.0;
    string gdindex_path_rerank = "";
    bool avoid_redundant_scene_results = false;
    bool gd_intra_normalization = false;
    bool gd_unbinarized = false;
    float score_den_power_norm = 0.5;
    float score_den_power_norm_rerank = 0.5;
    
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
            } else if ((!strcmp(argv[count_arg], "--query_index"))) {
                query_index_path = string(argv[count_arg + 1]);
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
            } else if ((!strcmp(argv[count_arg], "--asym_scoring_mode"))) {
                asym_scoring_mode = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if ((!strcmp(argv[count_arg], "--asym_scoring_mode_rerank"))) {
                asym_scoring_mode_rerank = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--score_den_power_norm")) {
                score_den_power_norm = atof(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--score_den_power_norm_rerank")) {
                score_den_power_norm_rerank = atof(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--word_selection_thresh")) {
                word_selection_thresh = atof(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--word_selection_thresh_rerank")) {
                word_selection_thresh_rerank = atof(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--group_lists_rerank_path")) {
                group_lists_rerank_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--gdindex_path_rerank")) {
                gdindex_path_rerank = string(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--avoid_redundant_scene_results")) {
                avoid_redundant_scene_results = true;
            } else if (!strcmp(argv[count_arg], "--gd_intra_normalization")) {
                gd_intra_normalization = true;
            } else if (!strcmp(argv[count_arg], "--gd_unbinarized")) {
                gd_unbinarized = true;
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
        cout << "------>query_index_path = " << query_index_path  << endl;
        cout << "------>feat_mode = " << feat_mode  << endl;
        cout << "------>gd_intra_normalization = " << gd_intra_normalization << endl;
        cout << "------>gd_unbinarized = " << gd_unbinarized << endl;
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
        cout << "------>group_lists_rerank_path = " << group_lists_rerank_path << endl;
        cout << "------>asym_scoring_mode = " << asym_scoring_mode  << endl;
        cout << "------>asym_scoring_mode_rerank = " << asym_scoring_mode_rerank  << endl;
        cout << "------>word_selection_thresh = " << word_selection_thresh << endl;
        cout << "------>word_selection_thresh_rerank = " << word_selection_thresh_rerank << endl;
        cout << "------>score_den_power_norm = " << score_den_power_norm << endl;
        cout << "------>score_den_power_norm_rerank = " << score_den_power_norm_rerank << endl;
        cout << "------>gdindex_path_rerank = " << gdindex_path_rerank << endl;
        cout << "------>avoid_redundant_scene_results = " << avoid_redundant_scene_results << endl;
    }

    // Instantiate retriever
    Retriever r;

    // Setting parameters...
    r.set_verbose_level(verbose_level);
    r.set_number_output_results(number_output_results);

    // Starting retrieval
    vector < vector < uint > > group_lists_rerank;
    if (group_lists_rerank_path != "") {
        // Parse rerank lists
        parse_agg_list(group_lists_rerank_path, 0, group_lists_rerank);            
    }
    r.retrieve_on_specific_dataset(gdindex_path,
                                   gdindex_parameters_path,
                                   feat_mode,
                                   db_list_path, 
                                   query_index_path,
                                   query_list_path,
                                   output_base_path,
                                   keyframe_numbers_path, 
                                   shot_list_path, 
                                   shot_mode, 
                                   number_scenes_rerank, 
                                   number_centroids,
                                   number_centroids_rerank, 
                                   group_lists_rerank, 
                                   word_selection_mode,
                                   word_selection_thresh,
                                   word_selection_thresh_rerank,
                                   min_number_words_visited,
                                   gdindex_path_rerank, 
                                   avoid_redundant_scene_results,
                                   gd_intra_normalization,
                                   gd_unbinarized,
                                   asym_scoring_mode,
                                   asym_scoring_mode_rerank,
                                   score_den_power_norm,
                                   score_den_power_norm_rerank);
        
    return EXIT_SUCCESS;
}

void parse_agg_list(const string agg_list_path, const int verbose_level, 
                    vector < vector < uint > >& agg_list_vec) {
    // Get all lines
    vector <string> file_lines;
    get_vector_of_strings_from_file_lines(agg_list_path, file_lines);
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
            uint item_to_add = static_cast<uint>(str2int(line_elements.at(count_el)));
            this_vec.push_back(item_to_add);
            if (verbose_level >= 3) cout << "Added  " << item_to_add << endl;
        }
        agg_list_vec.push_back(this_vec);
    }
}
void get_vector_of_strings_from_file_lines(const string file_name,
                                           vector<string>& out) {
    ifstream in_file(file_name.c_str());
    string line;
    out.clear();
    while (!in_file.eof()) {
        if (getline(in_file, line)) out.push_back(line);
    }
}
int str2int(const string string){
    stringstream parser(string);
    int i;
    parser>>i;
    return i;
}


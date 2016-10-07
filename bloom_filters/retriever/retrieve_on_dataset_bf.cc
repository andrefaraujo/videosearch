/**********************************************************
This program will generate retrieval results for a specific
dataset, with a specified set of parameters
**********************************************************/

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../bfindex.h"
#include "../point_indexed/point_index_io.h"

using namespace std;

typedef unsigned int uint;

const uint RESIDUAL_LENGTH = 32;
const string STR_INDEX_1 = "_keyframes.";
const string STR_INDEX_2 = "_bfv_point_idx_k";
const int SIFT_MODE = 0;
const string SIFT_NAME = "sift";
const int SIFTGEO_MODE = 1;
const string SIFTGEO_NAME = "siftgeo";

void get_vector_of_strings_from_file_lines(const string file_name,
                                           vector<string>& out);
void get_index_paths_from_clip_paths(const vector<string>& in,
                                     const string feat_name,
                                     const uint number_gaussians,
                                     vector<string>& out);
void get_index_path_from_query_path(const string query_list_path,
                                    const string feat_name,
                                    const uint number_gaussians,
                                    string& query_index_path);

void usage() {
    cout << "Perform retrieval using a specific dataset, using BF for indexing each video clip" << endl;
    cout << "Usage:" << endl;
    cout << "./retrieve_on_dataset_bf [options] --list_path clip_list_path --query_list query_list_path --output out_dir --n_bits number_bits --n_hashers number_hashers" << endl;
    cout << "For example:" << endl;
    cout << "./retrieve_on_dataset_bf --list_path /path/to/clip/list --query_list /path/to/query/list --output results_dir/ --n_bits 12 --n_hashers 512" << endl;
    cout << "The system will output 2 different files (one with log, one with results), using as the dir from out_dir." << endl;
    cout << "Options:" << endl;
    cout << "--feature_mode ARG: feature (local descriptor) mode (default: 0 = SIFT); another option is 1 = SIFTGEO" << endl;
    cout << "--number_gaussians ARG: number of Gaussians to use in point-indexed signatures (default: 512)" << endl;
    cout << "--results_per_query ARG: size of list of results to output (default: 100)" << endl;
    cout << "--alpha ARG: power normalization used in TF-IDF scoring (default: 0.75)" << endl;
    cout << "--verbose_level ARG: (default: 1)" << endl;

}

int main(int argc, char* * argv) {
    // Seeding random number generator
    srand(static_cast<uint>(time(0)));

    // Mandatory arguments
    string clip_list_path = "";
    string query_list_path = "";
    string out_dir = "";
    size_t number_bits = 0;
    size_t number_hashers = 0;
    
    // Default values for options
    int feat_mode = 0;
    uint number_gaussians = 512;
    uint results_per_query = 100;
    float alpha = 0.75;
    int verbose_level = 1;
    
    if (argc < 11) {
        cout << "Wrong usage!!!" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "***********************************" << endl;
        cout << "See usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    } else {
        for (int count_arg = 1; count_arg < argc; count_arg++) {
            if (!strcmp(argv[count_arg], "--list_path")) {
                clip_list_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--query_list")) {
                query_list_path = string(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--output")) {
                out_dir = string(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--n_bits")) {
                number_bits = static_cast<size_t>(atoi(argv[count_arg + 1]));
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--n_hashers")) {
                number_hashers = static_cast<size_t>(atoi(argv[count_arg + 1]));
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--feature_mode")) {
                feat_mode = atoi(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--results_per_query")) {
                results_per_query = static_cast<uint>(atoi(argv[count_arg + 1]));
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--number_gaussians")) {
                number_gaussians = static_cast<uint>(atoi(argv[count_arg + 1]));
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--alpha")) {
                alpha = atof(argv[count_arg + 1]);
                count_arg++;
            } else if (!strcmp(argv[count_arg], "--verbose_level")) {
                verbose_level = atoi(argv[count_arg + 1]);
                count_arg++;
            } else {
                cout << "Unrecognized argument " << argv[count_arg] 
                     << " , quitting..." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    // Check that all needed arguments were provided
    if (clip_list_path == "") {
        cout << "clip_list_path argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (query_list_path == "") {
        cout << "query_list_path argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (out_dir == "") {
        cout << "out_dir argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (number_bits == 0) {
        cout << "number_bits argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }
    if (number_hashers == 0) {
        cout << "number_hashers argument was not provided... see usage below:" << endl;
        usage();
        exit(EXIT_FAILURE);
    }

    if (verbose_level) {
        cout << "Starting evaluation using:" << endl;
        cout << "------>clip_list_path = " << clip_list_path  << endl;
        cout << "------>query_list_path = " << query_list_path  << endl;
        cout << "------>out_dir = " << out_dir  << endl;
        cout << "------>number_bits = " << number_bits  << endl;
        cout << "------>number_hashers = " << number_hashers  << endl;
        cout << "------>feat_mode = " << feat_mode  << endl;
        cout << "------>number_gaussians = " << number_gaussians << endl;
        cout << "------>results_per_query = " << results_per_query  << endl;
        cout << "------>alpha = " << alpha  << endl;
        cout << "------>verbose_level = " << verbose_level  << endl;
    }

    // Open log file
    ofstream log_file;
    string out_log = out_dir + "/out_log.txt";
    log_file.open(out_log.c_str());

    // Get feature extension
    string feat_name = "";
    if (feat_mode == SIFT_MODE) {
        feat_name = SIFT_NAME;
    } else if (feat_mode == SIFTGEO_MODE) {
        feat_name = SIFTGEO_NAME;
    } else {
        cout << "Error! feat_mode = " << feat_mode
             << " is not supported" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Instantiate & initialize BFIndex
    log_file << "main: Initializing BFIndex..." << endl;
    BFIndex b;    
    b.initialize(number_bits, number_hashers, RESIDUAL_LENGTH,
                 alpha, verbose_level);
    log_file << "main: done!" << endl;

    // Insert items into BFIndex
    // -- first, get appropriate index names
    vector<string> clip_paths, index_paths;
    get_vector_of_strings_from_file_lines(clip_list_path,
                                          clip_paths);    
    get_index_paths_from_clip_paths(clip_paths, feat_name,
                                    number_gaussians, index_paths);
    // -- second, do insertion
    log_file << "main: Inserting items into BFIndex..." << endl;
    b.insert_from_indexes(index_paths);
    log_file << "main: done!" << endl;

    // Open results files
    string out_results = out_dir + "/out_results.txt";
    ofstream results_file;
    results_file.open(out_results.c_str());
    
    // Process queries
    // -- first, load query index
    string query_index_path;
    get_index_path_from_query_path(query_list_path, feat_name,
                                   number_gaussians, query_index_path);
    vector < vector < uint > > vec_feat_assgns;
    vector < vector < float > > dummy;
    vector < vector < uint > > vec_feat_residuals_binarized;
    log_file << "main: Reading query point indexes..." << endl;
    read_bfv_point_index(query_index_path, vec_feat_assgns,
                         dummy, vec_feat_residuals_binarized);
    size_t number_queries = vec_feat_residuals_binarized.size();
    assert(number_queries == vec_feat_assgns.size());
    log_file << "main: done!" << endl;
    for (uint count_query = 0; count_query < number_queries; count_query++) {
        log_file << "main: Starting query " << count_query << "... (0-indexed)" << endl;
        // Search database
        vector< pair<float,uint> > results;
        log_file << "main: Querying..." << endl;
        b.perform_query(vec_feat_residuals_binarized.at(count_query),
                        vec_feat_assgns.at(count_query), results);
        log_file << "main: done!" << endl;

        // Write results to output files
        results_file << "Query " << count_query << endl;
        uint number_to_write = 
            min(results_per_query, static_cast<uint>(results.size()));
        for (size_t r = 0; r < number_to_write; r++) {
            // -- log file
            log_file << "#" << r+1 << ": "
                     << clip_paths.at(results.at(r).second)
                     << ", score = " << results.at(r).first
                     << endl;

            // -- results file
            results_file << clip_paths.at(results.at(r).second) << endl;
        }
        log_file << "main: Finished query " << count_query << "! (0-indexed)" << endl;
    }

    // Close files
    log_file.close();
    results_file.close();
    
    return EXIT_SUCCESS;
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
void get_index_paths_from_clip_paths(const vector<string>& in,
                                     const string feat_name,
                                     const uint number_gaussians,
                                     vector<string>& out) {
    size_t number_items = in.size();
    out.resize(number_items);
    for (size_t i = 0; i < number_items; i++) {
        int last_dot = in.at(i).find_last_of(".");
        string out_name = in.at(i).substr(0, last_dot)
            + STR_INDEX_1 + feat_name + STR_INDEX_2
            + to_string(number_gaussians);
        out.at(i) = out_name;
    }
}
void get_index_path_from_query_path(const string query_list_path,
                                    const string feat_name,
                                    const uint number_gaussians,
                                    string& query_index_path) {
    int last_dot = query_list_path.find_last_of(".");
    query_index_path = query_list_path.substr(0, last_dot+1)
        + feat_name + STR_INDEX_2
        + to_string(number_gaussians);
}


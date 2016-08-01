#include <algorithm>
#include <cassert>
#include <fstream>
#include <sstream>
#include <time.h>
#include <stdlib.h>

#include "retriever.h"
#include "../indexer/global_descriptors/gdindex.h"

using namespace std;

/********************************
HELPER FUNCTIONS
********************************/

static bool file_exists(string filename)
{
    ifstream ifile(filename.c_str());
    return bool(ifile);
}

/********************************
PUBLIC FUNCTIONS
********************************/

Retriever::Retriever() {
    // GDIndex-related variables
    gdindex_ptr_ = NULL;
    gdindex_ptr_rerank_ = NULL;
    query_index_ptr_ = NULL;
    local_descriptor_mode_ = GDIndex::SIFT_LOCAL_DESCRIPTOR;
    number_gaussians_global_descriptor_ = GD_NUMBER_GAUSSIANS_DEFAULT;
    min_number_words_visited_ = MIN_NUMBER_WORDS_SELECTED_DEFAULT;
    word_selection_mode_ = WORD_SELECTION_MODE_DEFAULT;
    word_selection_thresh_ = WORD_SELECTION_THRESH_DEFAULT;
    gdindex_trained_parameters_path_ = "../indexer/global_descriptors/trained_parameters/";

    number_output_results_ = DEFAULT_NUMBER_OUTPUT_RESULTS;
    verbose_level_ = DEFAULT_VERBOSE_LEVEL;

    // Clear member vectors
    db_list_.clear();
    keyframe_ids_for_eval_.clear();
    shot_first_frames_.clear();
}

Retriever::~Retriever()
{
    // Clear data structures
    db_list_.clear();
    keyframe_ids_for_eval_.clear();
    shot_first_frames_.clear();

    if (gdindex_ptr_ != NULL) {
        delete gdindex_ptr_;
        gdindex_ptr_ = NULL;
    }
    if (gdindex_ptr_rerank_ != NULL) {
        delete gdindex_ptr_rerank_;
        gdindex_ptr_rerank_ = NULL;
    }
}

void Retriever::retrieve_on_specific_dataset(const string gdindex_path, 
                                             const string db_list_path, 
                                             const string query_index_path,
                                             const string query_list_path, 
                                             const string output_base_path,
                                             const string keyframe_numbers_path,
                                             const string shot_list_path,
                                             const int shot_mode,
                                             const uint number_scenes_to_rerank,
                                             const uint number_gaussians_rerank,
                                             const vector < vector < uint > >& group_lists_rerank,
                                             const float word_selection_thresh_rerank,
                                             const string gdindex_path_rerank,
                                             const bool avoid_redundant_scene_results) {
    // Open files that will be written: log file and results file
    string log_file_name = output_base_path + "_log.txt";
    log_file_.open(log_file_name.c_str());
    log_file_ << "Starting retrieving using list of DB images from " << db_list_path << endl;

    ofstream results_file;
    string results_file_name = output_base_path + "_results.txt";
    results_file.open(results_file_name.c_str());     

    // Get queries
    vector<string> query_list;
    get_vector_of_strings_from_file_lines(query_list_path, query_list);
    uint number_queries = query_list.size();

    // Get list of database frames
    get_vector_of_strings_from_file_lines(db_list_path, db_list_);
    uint number_keyframes_in_database = db_list_.size();
    log_file_ << "This database contains " << number_keyframes_in_database << " frames" << endl;
    log_file_ << "The query set contains " << number_queries << " query images" << endl;

    // Set keyframe_ids_for_eval_, helpful for the usage of different types of indexing
    if ((shot_mode == GDIndex::SHOT_MODE_INDEP_KEYF) && (keyframe_numbers_path != "")) {
        // Shots with aggregation per frames in shots; in this case, we just need to 
        // use the correct indices for the frames; other than that, the processing works 
        // as usual
        get_vector_of_uints_from_file_lines(keyframe_numbers_path, 
                                            keyframe_ids_for_eval_);
    } else if ((shot_mode == GDIndex::SHOT_MODE_SHOT_AGG || shot_mode == GDIndex::SHOT_MODE_GLOBAL_AGG
                || shot_mode == GDIndex::SHOT_MODE_TRACK_AGG) && (shot_list_path != "")) {
        // Shots with one signature per shot, which aggregates features from many
        // frames; in this case, we'll use shot numbers in "keyframe_ids_for_eval_"
        // and we'll process the results from querying separately
        get_vector_of_uints_from_file_lines(shot_list_path,
                                            shot_first_frames_);
        for (uint count = 0; count < shot_first_frames_.size(); count++) {
            // Include indices corresponding to shot numbers
            keyframe_ids_for_eval_.push_back(count);
        }
    } else if (shot_mode == -1) {
        for (uint count_keyf = 0; count_keyf < number_keyframes_in_database; count_keyf++) {
            keyframe_ids_for_eval_.push_back(count_keyf);
        }
    } else {
        cout << "Problem! The combination of shot parameters you used is not allowed. Quitting..." << endl;
        exit(EXIT_FAILURE);
    }

    // Instantiate GDIndex, set parameters
    if (verbose_level_ >= 2) cout << "Creating GDIndex database..." << endl;
    gdindex_ptr_ = new GDIndex();
    uint ld_length, ld_frame_length;
    string ld_extension, ld_name;
    if (local_descriptor_mode_ == GDIndex::SIFT_LOCAL_DESCRIPTOR) {
        ld_length = GDIndex::SIFT_LENGTH;
        ld_frame_length = GDIndex::SIFT_FRAME_LENGTH;
        ld_extension = GDIndex::SIFT_EXTENSION;
        ld_name = GDIndex::SIFT_NAME;
    } else if (local_descriptor_mode_ == GDIndex::SIFTGEO_LOCAL_DESCRIPTOR) {
        ld_length = GDIndex::SIFTGEO_LENGTH;
        ld_frame_length = GDIndex::SIFTGEO_FRAME_LENGTH;
        ld_extension = GDIndex::SIFTGEO_EXTENSION;
        ld_name = GDIndex::SIFTGEO_NAME;
    } else {
        cout << "Problem! local_descriptors_mode_ = " 
             << local_descriptor_mode_
             << " is not recognized"
             << endl;
        exit(EXIT_FAILURE);
    }
    gdindex_ptr_->set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                       GDIndex::LD_PCA_DIM, GDIndex::LD_PRE_PCA_POWER, number_gaussians_global_descriptor_,
                                       GDIndex::GD_POWER, gdindex_trained_parameters_path_,
                                       verbose_level_);
    gdindex_ptr_->set_query_parameters(min_number_words_visited_, word_selection_mode_,
                                       word_selection_thresh_, gdindex_trained_parameters_path_,
                                       verbose_level_);

    // We instantiate another GDIndex object, in case we're using a two-step
    // scoring (eg, scene + shot)
    if (number_scenes_to_rerank != 0) {
        gdindex_ptr_rerank_ = new GDIndex();
        gdindex_ptr_rerank_->set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                                  GDIndex::LD_PCA_DIM, GDIndex::LD_PRE_PCA_POWER, number_gaussians_rerank,
                                                  GDIndex::GD_POWER, gdindex_trained_parameters_path_,
                                                  verbose_level_);
        gdindex_ptr_rerank_->set_query_parameters(min_number_words_visited_, word_selection_mode_,
                                                  word_selection_thresh_rerank, gdindex_trained_parameters_path_,
                                                  verbose_level_);
    }
    if (verbose_level_ >= 2) cout << "done!" << endl;
    if (verbose_level_ >= 2) cout << "It will use feature extension = " 
                                  << ld_extension << endl;
    log_file_ << "GDIndex database was instantiated. It will use feature extension = " 
              << ld_extension << endl;

    // Instantiate query index, if necessary
    if (query_index_path != "") {
        if (verbose_level_ >= 2) cout << "Creating GDIndex-query..." << endl;
        query_index_ptr_ = new GDIndex();
        query_index_ptr_->set_index_parameters(ld_length, ld_frame_length, ld_extension, ld_name,
                                               GDIndex::LD_PCA_DIM, GDIndex::LD_PRE_PCA_POWER, number_gaussians_global_descriptor_,
                                               GDIndex::GD_POWER, gdindex_trained_parameters_path_,
                                               verbose_level_);
        query_index_ptr_->set_query_parameters(min_number_words_visited_, word_selection_mode_,
                                               word_selection_thresh_, gdindex_trained_parameters_path_,
                                               verbose_level_);        
        if (verbose_level_ >= 2) cout << "done!" << endl;
        log_file_ << "GDIndex-query was instantiated" << endl;
    }

    // Reading relevant indexes
    gdindex_ptr_->read(gdindex_path);
    if (verbose_level_ >= 2) cout << "GDIndex index was loaded from " << gdindex_path << endl;
    log_file_ << "GDIndex index was loaded from " << gdindex_path << endl;

    if (number_scenes_to_rerank != 0 && gdindex_path_rerank != "") {
        gdindex_ptr_rerank_->read(gdindex_path_rerank);
        if (verbose_level_ >= 2) cout << "GDIndex-rerank index was loaded from " << gdindex_path_rerank << endl;
        log_file_ << "GDIndex-rerank index was loaded from " << gdindex_path_rerank << endl;
    }

    if (query_index_path != "") {
        query_index_ptr_->read(query_index_path);
        if (verbose_level_ >= 2) cout << "GDIndex-query index was loaded from " << query_index_path << endl;
        log_file_ << "GDIndex-query index was loaded from " << query_index_path << endl;
        if (query_index_ptr_->get_number_global_descriptors() != number_queries) {
            cout << "Error: number of query signatures is different from number of queries" << endl;
            cout << "Quitting..." << endl;

            log_file_ << "Error: number of query signatures is different from number of queries" << endl;
            log_file_ << "Quitting..." << endl;            
            exit(EXIT_FAILURE);
        }
    }

    if (shot_mode == -1) assert(number_keyframes_in_database 
                                == gdindex_ptr_->get_number_global_descriptors());

    for (unsigned int count_query = 0; count_query < number_queries; count_query++) {
        string this_query_name = query_list.at(count_query);

        log_file_ << "************************************" << endl;
        log_file_ << "************************************" << endl;
        log_file_ << "Starting query " << count_query
                  << " (first one was 0, not 1)" << endl;
        log_file_ << "Image file is " << this_query_name << endl;
        if (verbose_level_ >= 3) cout << "Starting query " << count_query << endl;

        // Declaring duration variables
        float query_duration = 0;

        // Declaring results variables for clip-based scoring
        vector< pair<float,uint> > results_query;

        string feature_file;
        if (query_index_path == "") {
            // In this case, query global descriptor will be computed in "perform_query"

            // Note: feature_file will have the full path
            feature_file = get_local_descriptor_filename(this_query_name,
                                                         ld_extension);
            if (!file_exists(feature_file)) {
                cout << "Feature file " << feature_file << " does not exist." << endl;
                cout << "Quitting..." << endl;

                log_file_ << "Feature file " << feature_file << " does not exist." << endl;
                log_file_ << "Quitting..." << endl;            
                exit(EXIT_FAILURE);
            }
            if (verbose_level_ >= 2) cout << "Using features from " << feature_file << endl;
            log_file_ << "Using features from " << feature_file << endl;
        }

        // Time search
        clock_t begin_clock = clock();

        if (verbose_level_ >= 2) cout << "Starting querying GDIndex..." << endl;
        log_file_ << "Starting querying GDIndex..." << endl;
        // Note: results_query will be sorted by scores
        gdindex_ptr_->perform_query(feature_file,
                                    query_index_ptr_,
                                    count_query,
                                    keyframe_ids_for_eval_,
                                    results_query, 
                                    number_scenes_to_rerank, 
                                    gdindex_ptr_rerank_,
                                    group_lists_rerank,
                                    verbose_level_);
        clock_t end_clock = clock();
        query_duration = double(end_clock - begin_clock)/CLOCKS_PER_SEC;
        if (verbose_level_ >= 2) cout << "done! Took "
                                      << query_duration << " secs." << endl;
        log_file_ << "done! Took "
                  << query_duration << " secs." << endl;


        if (!avoid_redundant_scene_results) {
            results_file << "Query " << count_query << endl;
            write_results(results_query, 
                          shot_mode,
                          results_file);
        } else {
            results_file << "Query " << count_query << endl;
            write_results_no_redundancy(results_query, 
                                        shot_mode,
                                        results_file);
        }

        if (verbose_level_) cout << "-----> Query " << count_query <<
                                ", query duration = " << query_duration <<
                                " secs" << endl;
        log_file_ << "-----> Query " << count_query <<
            ", query duration = " << query_duration <<
            " secs" << endl;
    }

    // Close files
    log_file_.close();
    results_file.close();
}

void Retriever::set_verbose_level(uint level) {
    verbose_level_ = level;
    if (verbose_level_ >= 3) cout << "Just set verbose_level_ to " 
                                  << verbose_level_ << endl;
}

void Retriever::set_local_descriptor_mode(int mode) {
    local_descriptor_mode_ = mode;
    if (verbose_level_ >= 3) cout << "Just set local_descriptor_mode_ to " << 
                                 local_descriptor_mode_ << endl;
}

void Retriever::set_number_output_results(uint n) {
    number_output_results_ = n;
    if (verbose_level_ >= 3) cout << "Just set number_output_results_ to " 
                                  << number_output_results_ << endl;
}

void Retriever::set_number_gaussians_global_descriptor(uint n) {
    number_gaussians_global_descriptor_ = n;
    if (verbose_level_ >= 3) cout << "Just set number_gaussians_global_descriptor_ to " 
                                  << number_gaussians_global_descriptor_ << endl;
}

void Retriever::set_gdindex_path(string path) {
    gdindex_trained_parameters_path_ = path;
    if (verbose_level_ >= 3) cout << "Just set gdindex_trained_parameters_path_ to " 
                                  << gdindex_trained_parameters_path_ << endl;
}

void Retriever::set_word_selection_mode(int mode) {
    word_selection_mode_ = mode;
    if (verbose_level_ >= 3) cout << "Just set word_selection_mode_ to " << 
                                 word_selection_mode_ << endl;
}

void Retriever::set_word_selection_thresh(float t) {
    word_selection_thresh_ = t;
    if (verbose_level_ >= 3) cout << "Just set word_selection_thresh_ to " << 
                                 word_selection_thresh_ << endl;
}

void Retriever::set_min_num_words_visited(uint n) {
    min_number_words_visited_ = n;
    if (verbose_level_ >= 3) cout << "Just set min_number_words_visited_ to " << 
                                 min_number_words_visited_ << endl;
}

/********************************
PRIVATE FUNCTIONS
********************************/


void Retriever::get_vector_of_strings_from_file_lines(const string file_name,
                                                      vector<string>& out) {
    ifstream in_file(file_name.c_str());
    string line;
    out.clear();
    while (in_file >> line) {
        out.push_back(line);
    }
}

void Retriever::get_vector_of_uints_from_file_lines(const string file_name,
                                                    vector<uint>& out) {
    ifstream in_file(file_name.c_str());
    string line;
    out.clear();
    while (in_file >> line) {
        istringstream s(line);
        uint val;
        s >> val;
        out.push_back(val);
    }
}

string Retriever::get_local_descriptor_filename(const string image_path,
                                                const string ext) {
    uint pos = image_path.find_last_of(".");
    return (image_path.substr(0, pos) + ext);
}

void Retriever::write_results(const vector< pair < float, uint > >& results,
                              const int shot_mode,
                              ofstream& file) {
    uint number_to_write = 
        min(number_output_results_, static_cast<uint>(results.size()));

    for (uint count_result = 0; count_result < number_to_write; count_result++) {
        uint result_number;
        if (shot_mode == GDIndex::SHOT_MODE_SHOT_AGG 
            || shot_mode == GDIndex::SHOT_MODE_GLOBAL_AGG
            || shot_mode == GDIndex::SHOT_MODE_TRACK_AGG) {
            result_number = 
                shot_first_frames_.at(results.at(count_result).second);
        } else {
            result_number = results.at(count_result).second;
        }
        float score = results.at(count_result).first;
        string path_to_write = db_list_.at(result_number);
        log_file_ << "#" << count_result + 1 << ", number " << result_number << ", file name is " << path_to_write << ", score = " << score << endl;
        file << path_to_write << endl;
    }
}

void Retriever::write_results_no_redundancy(const vector< pair < float, uint > >& results,
                                            const int shot_mode,
                                            ofstream& file) {
    vector<string> already_scored_clips;
    uint count_result = 0, count_used_result = 0;
    while (count_used_result < number_output_results_ 
           && count_result < results.size()) {
        uint result_number;
        if (shot_mode == GDIndex::SHOT_MODE_SHOT_AGG 
            || shot_mode == GDIndex::SHOT_MODE_GLOBAL_AGG
            || shot_mode == GDIndex::SHOT_MODE_TRACK_AGG) {
            result_number = 
                shot_first_frames_.at(results.at(count_result).second);
        } else {
            result_number = results.at(count_result).second;
        }
        string path_to_write = db_list_.at(result_number);
        size_t strpos = path_to_write.rfind('/');
        string result_path = path_to_write.substr(0, strpos);
        // If this clip has already showed up, we skip it
        vector<string>::iterator it;
        it = find(already_scored_clips.begin(), already_scored_clips.end(), result_path);
        if (it != already_scored_clips.end()) {
            count_result++;
            continue;
        }
        // Include this clip in already_scored_clips, so that the same clip
        // will not be scored twice.
        already_scored_clips.push_back(result_path);
        float score = results.at(count_result).first;
        log_file_ << "#" << count_result + 1 << ", number " << result_number << ", file name is " << path_to_write << ", score = " << score << endl;
        file << path_to_write << endl;
        count_result++;
        count_used_result++;
    }
}


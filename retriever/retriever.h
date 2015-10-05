/*****************************************
Class that performs retrieval given a dataset and index.
 *****************************************/
#ifndef RETRIEVER_H
#define RETRIEVER_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "../indexer/global_descriptors/gdindex.h"

using namespace std;

// Default initialization for some variables
const int DEFAULT_VERBOSE_LEVEL = 1;
const uint DEFAULT_NUMBER_OUTPUT_RESULTS = 100;

class Retriever {
public:

	// Constructor, some variables are initialized.
	Retriever();

	// Destructor, clean variables.
	~Retriever();

	void retrieve_on_specific_dataset(const string gdindex_path, 
                                      const string db_list_path, 
                                      const string query_list_path, 
                                      const string output_base_path, 
                                      const string keyframe_numbers_path = "", 
                                      const string shot_list_path = "", 
                                      const int shot_mode = -1, 
                                      const uint number_scenes_to_rerank = 0, 
                                      const uint rerank_number_gaussians = 0,  
                                      const vector < vector < uint > >& group_lists = vector < vector < uint > >(), 
                                      const vector < pair < string, pair < uint, uint > > >& shot_info = vector < pair < string, pair < uint, uint > > >(), 
                                      const float word_selection_thresh_rerank = 6, 
                                      const string gdindex_path_other = "", 
                                      const bool avoid_redundant_scene_results = true);

	/**********************************************************
        Setting functions
	 *********************************************************/

	// Set verbose level.
	void set_verbose_level(uint level);

	// Set local descriptor type, which will determine its other properties
	void set_local_descriptor_mode(int mode);

	// Set number of results to output
	void set_number_output_results(uint n);

	// Set number of gaussians to use in global descriptor
	void set_number_gaussians_global_descriptor(uint n);

	// Set path where to find trained GDIndex parameters
	void set_gdindex_path(string path);

    // Set word selection mode
    void set_word_selection_mode(int mode);

    // Set word selection threshold
    void set_word_selection_thresh(float t);

    // Set min number words visited
    void set_min_num_words_visited(uint n);

private:
	/************ Variables *************/
	// GDIndex database, pointer to it
	GDIndex* gdindex_ptr_;
	// Extra GDIndex database, pointer to it; This is used in the case we
	// are retrieving in two stages (eg, first with scenes signatures, then
    // with shot signatures)
	GDIndex* gdindex_ptr_rerank_;

    // Variables related to GDIndex member variables that will be set
    int local_descriptor_mode_;
	uint number_gaussians_global_descriptor_;
    uint min_number_words_visited_;
    int word_selection_mode_;
    float word_selection_thresh_;
	string gdindex_trained_parameters_path_;

    // Number of results to output
	uint number_output_results_;

	// Verbose level
	int verbose_level_;

	// The vector db_list contains the path for each frame in the database.
	vector<string> db_list_;

	// This vector is used when we need to assign different numbers to the database instances, when doing evaluation; this has been used when dealing with different frame rates
	vector<unsigned int> keyframe_ids_for_eval_; // note: starts at 0

	// This vector is used when we need to load shot indices, to do retrieval using this information
	vector<unsigned int> shot_first_frames_;

	// Ofstream which will write log. It should be opened and closed using the same function,
	// but having it as a member variable further allows writing in any function.
	ofstream log_file_;

	/************ Functions *************/
	// Get vector of strings from a text file, one entry per line
	void get_vector_of_strings_from_file_lines(const string file_name,
                                               vector<string>& out);

	// Get vector of unsigned ints from a text file, one entry per line
	void get_vector_of_uints_from_file_lines(const string file_name,
                                             vector<uint>& out);

	// Get local descriptor filename from image path
	string get_local_descriptor_filename(const string image_path, const string ext);

	// Write query results for a given query
    void write_results(const vector< pair < float, uint > >& results,
                       const int shot_mode,
                       ofstream& file);

	// Write query results for a given query, avoiding redundant clips
    void write_results_no_redundancy(const vector< pair < float, uint > >& results,
                                     const int shot_mode,
                                     ofstream& file);
};

#endif

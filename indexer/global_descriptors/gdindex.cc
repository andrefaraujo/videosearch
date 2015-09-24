#include <algorithm>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include "../../common/feature_set/feature_set.h"
#include "gdindex.h"

extern "C" {
#include "../../common/yael_v438_modif/yael/gmm.h"
#include "../../common/yael_v438_modif/yael/matrix.h"
}

GDIndex::GDIndex() {
    // Initializing member variables to default values
    // -- Index
    index_.number_global_descriptors = 0;
    index_.word_descriptor.clear();
    index_.word_l1_norms.clear();
    index_.word_total_soft_assignment.clear();
    index_.frame_numbers_in_db.clear();
    // -- Index parameters
    index_parameters_.ld_length = LD_LENGTH_DEFAULT;
    index_parameters_.ld_frame_length = LD_FRAME_LENGTH_DEFAULT;
    index_parameters_.ld_extension = LD_EXTENSION_DEFAULT;
    index_parameters_.ld_name = LD_NAME_DEFAULT;
    index_parameters_.ld_pca_dim = LD_PCA_DIM_DEFAULT;
    index_parameters_.ld_pre_pca_power = LD_PRE_PCA_POWER_DEFAULT;
    index_parameters_.ld_mean_vector = NULL;
    index_parameters_.ld_pca_eigenvectors.clear();
    index_parameters_.gd_gmm = NULL;
    index_parameters_.gd_number_gaussians = GD_NUMBER_GAUSSIANS_DEFAULT;
    index_parameters_.gd_power = GD_POWER_DEFAULT;
    // -- Query parameters
    query_parameters_.min_number_words_visited = MIN_NUMBER_WORDS_VISITED_DEFAULT;
    query_parameters_.word_selection_mode = WORD_SELECTION_MODE_DEFAULT;
    query_parameters_.word_selection_thresh = WORD_SELECTION_THRESH_DEFAULT;
    query_parameters_.fast_corr_weights = NULL;
}

GDIndex::~GDIndex() {
    // Clearing and deleting vectors/pointers
    // -- Index
    index_.word_descriptor.clear();
    index_.word_l1_norms.clear();
    index_.word_total_soft_assignment.clear();
    index_.frame_numbers_in_db.clear();
    // -- Index parameters
    if (index_parameters_.ld_mean_vector != NULL) {
        delete[] index_parameters_.ld_mean_vector;
        index_parameters_.ld_mean_vector = NULL;
    }
    for (uint i = 0; i < index_parameters_.ld_pca_eigenvectors.size(); i++) {
        if (index_parameters_.ld_pca_eigenvectors.at(i) != NULL) {
            delete [] index_parameters_.ld_pca_eigenvectors.at(i);
            index_parameters_.ld_pca_eigenvectors.at(i) = NULL;
        }
    }
    index_parameters_.ld_pca_eigenvectors.clear();
    if (index_parameters_.gd_gmm != NULL) {
        gmm_delete(index_parameters_.gd_gmm);
    }
    // -- Query parameters
    if (query_parameters_.fast_corr_weights != NULL) {
        delete [] query_parameters_.fast_corr_weights;
        query_parameters_.fast_corr_weights = NULL;
    }
}

void GDIndex::write(const string index_path) {

}

void GDIndex::read(const string index_path) {

}

void GDIndex::write_frame_list(const string file_path) {
    
}

void GDIndex::clean_index() {

}

uint GDIndex::get_number_global_descriptors() {
    return index_.number_global_descriptors;
}

void GDIndex::generate_index(const vector<string>& feature_files, 
                             const int verbose_level) {

}

void GDIndex::generate_index_shot_based(const vector<string>& feature_files, 
                                        const vector<uint>& shot_beg_frames,
                                        const int shot_mode, const int shot_keyf, 
                                        const vector < vector < 
                                            pair < uint, uint > > >& track_lists,
                                        const int verbose_level) {

}

void GDIndex::generate_global_descriptor(const FeatureSet* feature_set, 
                                         vector<uint>& gd_word_descriptor, 
                                         vector<float>& gd_word_l1_norm, 
                                         vector<float>& gd_word_total_soft_assignment) {

}

void GDIndex::performQuery(const string local_descriptors_path, 
                           vector< pair<float,uint> >& results, 
                           const vector<uint>& indices,
                           const uint num_scenes_to_rerank,
                           const uint group_testing_number_centroids ,
                           const GDIndex* revv_other_ptr,
                           const vector < vector < uint > >& vGroupLists,
                           const vector < pair < string, pair < uint, uint > > >& shot_info,
                           const int verbose_level) {

}

void GDIndex::set_index_parameters(const uint ld_length, const uint ld_frame_length,
                                   const string ld_extension, const string ld_name,
                                   const uint ld_pca_dim, const float ld_pre_pca_power,
                                   const uint gd_number_gaussians, const float gd_power,
                                   const string trained_parameters_path,
                                   const int verbose_level) {

}

void GDIndex::set_query_parameters(const uint min_number_words_visited,
                                   const int word_selection_mode,
                                   const float word_selection_thresh,
                                   const string trained_parameters_path,
                                   const int verbose_level) {

}

void GDIndex::update_index() {

}

void GDIndex::sign_binarize(const vector<float>& gd_word_residuals, 
                            vector<uint>& gd_word_descriptor) {

}

void GDIndex::project_local_descriptr_pca(const float* desc, float* pca_desc) {

}

void GDIndex::sampleFramesFromShot(const uint number_frames_out, 
                                   const uint first_frame, 
                                   const uint number_frames_this_shot, 
                                   vector<uint>& out_frames) {

}

void GDIndex::query(const vector<uint>& query_word_descriptor,
                    const vector<float>& query_word_l1_norm,
                    const vector<float>& query_word_total_soft_assignment,
                    const vector<uint>& database_indices,
                    vector< pair<float,uint> >& database_scores_indices) {

}

void GDIndex::load_ld_mean_vector(string path) {

}

void GDIndex::load_ld_pca_eigenvectors(string path) {

}

void GDIndex::load_gd_gmm(string path) {

}

void GDIndex::load_corr_weights(string path) {

}


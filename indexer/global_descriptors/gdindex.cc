#include <algorithm>
#include <cassert>
#include <float.h>
#include <fstream>
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

/********************************
PUBLIC FUNCTIONS
********************************/

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
    query_parameters_.min_number_words_selected = MIN_NUMBER_WORDS_SELECTED_DEFAULT;
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
    int n_write; // aux. variable for writing

    // Open file for writing
    FILE* index_file = fopen(index_path.c_str(), "wb");
    if (index_file == NULL) {
        fprintf(stderr, "GDIndex::write : Cannot open: %s\n", index_path.c_str());
        exit(EXIT_FAILURE);
    }

    // Write number of global descriptors
    int number_gd_to_write = static_cast<int>(index_.number_global_descriptors);
    n_write = fwrite(&number_gd_to_write, sizeof(int), 1, index_file);
    
    // Assert that we have data in index_ and that it makes sense, 
    // and allocate helper variables
    assert(index_.word_descriptor.size() != 0);
    assert(index_.word_l1_norms.size() != 0);
    assert(index_.word_total_soft_assignment.size() != 0);
    assert(index_.word_descriptor.size() == index_.word_l1_norms.size());
    assert(index_.word_descriptor.size() == index_.word_total_soft_assignment.size());
    uint* word_descriptor_to_write = 
        new uint[index_parameters_.gd_number_gaussians];
    float* word_l1_norms_to_write = 
        new float[index_parameters_.gd_number_gaussians];
    float* word_total_soft_assignment_to_write = 
        new float[index_parameters_.gd_number_gaussians];

    // Loop over items in index and write them out
    for (int count_item = 0; count_item < number_gd_to_write; count_item++) {
        // Collect data that will be written
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians; 
             count_gaussian++) {
            word_l1_norms_to_write[count_gaussian] = 
                index_.word_l1_norms.at(count_item).at(count_gaussian);
            word_total_soft_assignment_to_write[count_gaussian] = 
                index_.word_total_soft_assignment.at(count_item).at(count_gaussian);
            word_descriptor_to_write[count_gaussian] = 
                index_.word_descriptor.at(count_item).at(count_gaussian);
        }
        // Write to file
        n_write = fwrite(word_l1_norms_to_write, sizeof(float), 
                         index_parameters_.gd_number_gaussians, index_file);
        n_write = fwrite(word_total_soft_assignment_to_write, sizeof(float), 
                         index_parameters_.gd_number_gaussians, index_file);
        n_write = fwrite(word_descriptor_to_write, sizeof(uint), 
                         index_parameters_.gd_number_gaussians, index_file);
    }

    // Clean up
    if (word_l1_norms_to_write != NULL) {
        delete [] word_l1_norms_to_write;
        word_l1_norms_to_write = NULL;
    }
    if (word_total_soft_assignment_to_write != NULL) {
        delete [] word_total_soft_assignment_to_write;
        word_total_soft_assignment_to_write = NULL;
    }
    if (word_descriptor_to_write != NULL) {
        delete [] word_descriptor_to_write;
        word_descriptor_to_write = NULL;
    }

    // Close file
    fclose(index_file);
}

void GDIndex::read(const string index_path) {
    int n_read; // aux. variable for reading

    // Open file for reading
    FILE* index_file = fopen(index_path.c_str(), "rb");
    if (index_file == NULL) {
        fprintf(stderr, "GDIndex::read : Cannot open: %s\n", index_path.c_str());
        exit(EXIT_FAILURE);
    }

    // Read number of global descriptors
    int number_gd_to_read = 0;
    n_read = fread(&number_gd_to_read, sizeof(int), 1, index_file);

    // Size of current index
    int current_db_size = index_.number_global_descriptors;

    // Allocate helper variables
    uint* word_descriptor_to_read = 
        new uint[index_parameters_.gd_number_gaussians];
    float* word_l1_norms_to_read = 
        new float[index_parameters_.gd_number_gaussians];
    float* word_total_soft_assignment_to_read = 
        new float[index_parameters_.gd_number_gaussians];

    // Depending on query_parameters_.word_selection_mode, we will load
    // either l1 norms of total soft assignment information. So the resize
    // of index_ vectors is done only to one or the other
    index_.word_descriptor.resize(current_db_size + number_gd_to_read);
    if (query_parameters_.word_selection_mode == WORD_L1_NORM) {
        index_.word_l1_norms.resize(current_db_size + number_gd_to_read);
    } else if (query_parameters_.word_selection_mode == WORD_SOFT_ASSGN) {
        index_.word_total_soft_assignment.resize(current_db_size + number_gd_to_read);
    } else {
        cout << "Error! Mode " << query_parameters_.word_selection_mode 
             << " is not allowed. Quitting..."
             << endl;
        exit(EXIT_FAILURE);
    }
    
    // Loop over items, read and insert them into index_
    for (int count_item = current_db_size; 
         count_item < current_db_size + number_gd_to_read; 
         count_item++) {
        // Read data
        n_read = fread(word_l1_norms_to_read, sizeof(float), 
                       index_parameters_.gd_number_gaussians, index_file);
        n_read = fread(word_total_soft_assignment_to_read, sizeof(float), 
                       index_parameters_.gd_number_gaussians, index_file);
        n_read = fread(word_descriptor_to_read, sizeof(uint), 
                       index_parameters_.gd_number_gaussians, index_file);

        // Insert data into index_
        if (query_parameters_.word_selection_mode == WORD_L1_NORM) {
            index_.word_l1_norms.at(count_item)
                .resize(index_parameters_.gd_number_gaussians);
            for (uint count_gaussian = 0; 
                 count_gaussian < index_parameters_.gd_number_gaussians; 
                 count_gaussian++) {
                index_.word_l1_norms.at(count_item).at(count_gaussian)
                    = word_l1_norms_to_read[count_gaussian];
            }
        } else {
            index_.word_total_soft_assignment.at(count_item)
                .resize(index_parameters_.gd_number_gaussians);
            for (uint count_gaussian = 0; 
                 count_gaussian < index_parameters_.gd_number_gaussians; 
                 count_gaussian++) {
                index_.word_total_soft_assignment.at(count_item).at(count_gaussian)
                    = word_total_soft_assignment_to_read[count_gaussian];
            }
        }
        index_.word_descriptor.at(count_item)
            .resize(index_parameters_.gd_number_gaussians);
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians; 
             count_gaussian++) {
            index_.word_descriptor.at(count_item).at(count_gaussian)
                = word_descriptor_to_read[count_gaussian];
        }
    }

    // Clean up
    if (word_l1_norms_to_read != NULL) {
        delete [] word_l1_norms_to_read;
        word_l1_norms_to_read = NULL;
    }
    if (word_total_soft_assignment_to_read != NULL) {
        delete [] word_total_soft_assignment_to_read;
        word_total_soft_assignment_to_read = NULL;
    }
    if (word_descriptor_to_read != NULL) {
        delete [] word_descriptor_to_read;
        word_descriptor_to_read = NULL;
    }

    // Close file
    fclose(index_file);

    // Update other index_ variables, since now index has changed
    update_index();
}

void GDIndex::write_frame_list(const string file_path) {
    ofstream out_file;
    out_file.open(file_path.c_str());
    uint number_frames_in_db = index_.frame_numbers_in_db.size();
    for (uint count_line = 0; count_line < number_frames_in_db; count_line++) {
        out_file << index_.frame_numbers_in_db.at(count_line) << endl;
    }
    out_file.close();
}

void GDIndex::clean_index() {
    index_.word_descriptor.clear();
    index_.word_l1_norms.clear();
    index_.word_total_soft_assignment.clear();
    index_.frame_numbers_in_db.clear();

    update_index();
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
    // Local descriptor information
    index_parameters_.ld_length = ld_length;
    index_parameters_.ld_frame_length = ld_frame_length;
    index_parameters_.ld_extension = ld_extension;
    index_parameters_.ld_name = ld_name;

    // Parameters for PCA-ing local descriptors
    index_parameters_.ld_pca_dim = ld_pca_dim;
    index_parameters_.ld_pre_pca_power = ld_pre_pca_power;
    // -- LD mean vector is stored in descriptor covariance file
    char aux_mean_vector_path[1024];
    sprintf(aux_mean_vector_path, "%s/%s.pre_alpha.%.2f.desc_covariance", 
            trained_parameters_path.c_str(), 
            index_parameters_.ld_name.c_str(), 
            index_parameters_.ld_pre_pca_power);
    string ld_mean_vector_path = aux_mean_vector_path;
    load_ld_mean_vector(ld_mean_vector_path);

    char aux_ld_pca_eigenvectors[1024];
    sprintf(aux_ld_pca_eigenvectors, "%s/%s.pre_alpha.%.2f.desc_eigenvectors", 
            trained_parameters_path.c_str(), 
            index_parameters_.ld_name.c_str(), 
            index_parameters_.ld_pre_pca_power);
    string ld_pca_eigenvectors_path = aux_ld_pca_eigenvectors;
    load_ld_pca_eigenvectors(ld_pca_eigenvectors_path);

    // Parameters used for global descriptor computation
    index_parameters_.gd_number_gaussians = gd_number_gaussians;
    index_parameters_.gd_power = gd_power;

    char aux_gd_gmm[1024];
    sprintf(aux_gd_gmm, "%s/%s.pre_alpha.%.2f.pca.%d.gmm.%d",
            trained_parameters_path.c_str(), 
            index_parameters_.ld_name.c_str(), 
            index_parameters_.ld_pre_pca_power,
            index_parameters_.ld_pca_dim,
            index_parameters_.gd_number_gaussians);
    string gd_gmm_path = aux_gd_gmm;
    load_gd_gmm(gd_gmm_path);
}

void GDIndex::set_query_parameters(const uint min_number_words_selected,
                                   const int word_selection_mode,
                                   const float word_selection_thresh,
                                   const string trained_parameters_path,
                                   const int verbose_level) {
    query_parameters_.min_number_words_selected = min_number_words_selected;
    query_parameters_.word_selection_mode = word_selection_mode;
    query_parameters_.word_selection_thresh = word_selection_thresh;

    char aux_corr_weights[1024];
    sprintf(aux_corr_weights, "%s/%s.pre_alpha.%.2f.pca.%d.gmm.%d.pre_alpha.%.2f.corr_weights",
            trained_parameters_path.c_str(), 
            index_parameters_.ld_name.c_str(), 
            index_parameters_.ld_pre_pca_power,
            index_parameters_.ld_pca_dim,
            index_parameters_.gd_number_gaussians,
            index_parameters_.gd_power);
    string corr_weights_path = aux_corr_weights;
    load_corr_weights(corr_weights_path);

    // Precompute pop-counts for fast executions
    query_parameters_.pop_count[0] = 0;
    for (int i = 1; i < 65536; i++) {
        query_parameters_.pop_count[i] = query_parameters_.pop_count[i>>1] + (i&1);
    }
}

/********************************
PRIVATE FUNCTIONS
********************************/

void GDIndex::update_index() {
    // Update number of global descriptors stored
    index_.number_global_descriptors = index_.word_descriptor.size();

    // Update number of words selected and norm. factors for each db item
    // The purpose of doing this here is to have these numbers ready when
    // querying (and not need to recalculate them at query time)
    index_.number_words_selected.resize(index_.number_global_descriptors, 
                                        index_parameters_.gd_number_gaussians);
    index_.norm_factors.resize(index_.number_global_descriptors, 
                               sqrt(index_parameters_.gd_number_gaussians
                                    * index_parameters_.ld_pca_dim));
    
    // TODO(andrefaraujo): for future work, we might want to include here
    // the possibility of skipping Gaussian residuals on the database
    // side as well
}

void GDIndex::sign_binarize(const vector<float>& gd_word_residuals, 
                            vector<uint>& gd_word_descriptor) {
    gd_word_descriptor.resize(index_parameters_.gd_number_gaussians);
    for (uint count_gaussian = 0; 
         count_gaussian < index_parameters_.gd_number_gaussians;
         count_gaussian++) {
        uint packed_block = 0;
        uint start = count_gaussian*index_parameters_.ld_pca_dim;
        uint end = (count_gaussian + 1)*index_parameters_.ld_pca_dim;
        // Sign binarize
        for (uint count_dim = start; count_dim < end; count_dim++) {
            uint bit = (gd_word_residuals.at(count_dim) > 0) ? 1 : 0;
            PUSH_BIT(packed_block, bit);
        }
        gd_word_descriptor.at(count_gaussian) = packed_block;
    }
}

void GDIndex::project_local_descriptor_pca(const float* desc, float* pca_desc) {
    vector < float > desc_pow(index_parameters_.ld_length, 0);
	// Pre power law and normalization
    float l2_norm_sq = 0;
    for (uint count_in_dim = 0; count_in_dim < index_parameters_.ld_length; count_in_dim++) {
        POWER_LAW(desc[count_in_dim], index_parameters_.ld_pre_pca_power, 
                  desc_pow.at(count_in_dim));
        l2_norm_sq += desc_pow.at(count_in_dim) * desc_pow.at(count_in_dim);
    }
    float l2_norm = sqrt(l2_norm_sq);
    for (uint count_in_dim = 0; count_in_dim < index_parameters_.ld_length; count_in_dim++) {
        desc_pow.at(count_in_dim) /= l2_norm;
    }

	// Projection onto eigenvectors
    for (uint count_out_dim = 0; count_out_dim < index_parameters_.ld_pca_dim; 
         count_out_dim++) {
        pca_desc[count_out_dim] = 0;
        for (uint count_in_dim = 0; count_in_dim < index_parameters_.ld_length; 
             count_in_dim++) {
            pca_desc[count_out_dim] +=
                (desc_pow.at(count_in_dim) - index_parameters_.ld_mean_vector[count_in_dim])
                * index_parameters_.ld_pca_eigenvectors.at(count_out_dim)[count_in_dim];
        }
    }
}

void GDIndex::sample_frames_from_shot(const uint number_frames_out, 
                                      const uint first_frame, 
                                      const uint number_frames_this_shot, 
                                      vector<uint>& out_frames) {
    // Note: number_frames_this_shot must be bigger than number_frames_out
    out_frames.clear();
    if (number_frames_out == 1) {
        uint middle_frame = first_frame 
            + static_cast<uint>(floor(static_cast<float>(number_frames_this_shot)/2));
        out_frames.push_back(middle_frame);
    } else {
        double rate_float = static_cast<double>(number_frames_this_shot)
            /static_cast<double>(number_frames_out);
        uint rate;
        if ((number_frames_this_shot % number_frames_out) == 0) {
            rate = static_cast<uint>(rate_float);
            for (uint count_sample = 0; count_sample < number_frames_this_shot; 
                 count_sample += rate) {
                out_frames.push_back(first_frame + count_sample);
            }
        } else {
            if ((number_frames_out - 1)*static_cast<uint>(ceil(rate_float))
                < number_frames_this_shot) {
                rate = static_cast<uint>(ceil(rate_float));
                for (uint count_sample = 0; count_sample < number_frames_this_shot;
                     count_sample += rate) {
                    out_frames.push_back(first_frame + count_sample);
                }
            } else {
                rate = static_cast<uint>(floor(rate_float));
                vector<uint> out_frames_aux;
                for (uint count_sample = 0; count_sample < number_frames_this_shot;
                     count_sample += rate) {
                    out_frames_aux.push_back(first_frame + count_sample);
                }
                uint start_ind, end_ind;
                uint extra_frames = (out_frames_aux.size() - number_frames_out);

                if (extra_frames) {
                    // In this case, we'll keep the middle ones, in order not to unbalance
                    // the selection too much
                    if (extra_frames % 2) {
                        // Odd: keep one more at the beginning
                        start_ind = (extra_frames - 1)/2;
                        end_ind = out_frames_aux.size() - 1 
                            - static_cast<uint>(ceil((static_cast<double>(extra_frames)/2)));
                    } else {
                        // Even
                        start_ind = extra_frames/2;
                        end_ind = out_frames_aux.size() - 1 - extra_frames/2;
                    }

                    for (uint count_sample = start_ind; count_sample < end_ind + 1; 
                         count_sample++) {
                        out_frames.push_back(out_frames_aux.at(count_sample));
                    }
                } else {
                    out_frames = out_frames_aux;
                }
            }
        }
    }
    // Make sure we've collected the correct number
    assert(number_frames_out == out_frames.size());
}

void GDIndex::query(const vector<uint>& query_word_descriptor,
                    const vector<float>& query_word_l1_norm,
                    const vector<float>& query_word_total_soft_assignment,
                    const vector<uint>& database_indices,
                    vector< pair<float,uint> >& database_scores_indices) {
    assert(index_.number_global_descriptors == database_indices.size());
    // Resize vector that will be passed back
    database_scores_indices.resize(index_.number_global_descriptors);

    // Loop over database items and get their scores
    for (uint count_elem = 0; count_elem < index_.number_global_descriptors; 
         count_elem++) {
        uint number_item = database_indices.at(count_elem);
        database_scores_indices.at(count_elem).second = number_item;

        // Check that the database item has at least the minimum number of
        // selected words that we require; if not, just set it to FLT_MAX
        if (index_.number_words_selected.at(count_elem) 
            <= query_parameters_.min_number_words_selected) {
            database_scores_indices.at(count_elem).first = FLT_MAX;
            continue;
        }

        float query_norm = 0;
        float total_correlation = 0;
        uint query_number_words_selected = 0;
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians;
             count_gaussian++) {
            float query_word_strength = 0;
            if (query_parameters_.word_selection_mode == WORD_L1_NORM) {
                query_word_strength = query_word_l1_norm.at(count_gaussian);
            } else {
                query_word_strength = query_word_total_soft_assignment.at(count_gaussian);
            }

            // Skip in case query word strength is not enough
            if (query_word_strength <= query_parameters_.word_selection_thresh) {
                continue;
            } else {
                query_number_words_selected++;
                // Compute Hamming distance
                uint aux = query_word_descriptor.at(count_gaussian)
                    ^ index_.word_descriptor.at(count_elem).at(count_gaussian);
                uint h = query_parameters_.pop_count[aux >> 16]
                    + query_parameters_.pop_count[aux & 65535];
                // Add to correlation
                total_correlation += query_parameters_.fast_corr_weights[h];
            }
        }

        // Figure out query norm to be used
        query_norm = sqrt(query_number_words_selected
                          * index_parameters_.ld_pca_dim);

        // Compute final correlation for this item
        total_correlation /= index_.norm_factors.at(count_elem) * query_norm;

        // Change sign such that smaller is better
        database_scores_indices.at(count_elem).first = -total_correlation;
    }
    // Sort scores
    sort(database_scores_indices.begin(), database_scores_indices.end(), 
         cmp_float_uint_ascend);
}

void GDIndex::load_ld_mean_vector(string path) {
    int n_read;
    FILE* mv_file = fopen(path.c_str(), "rb");
    if (mv_file == NULL) {
        printf("GDIndex::load_ld_mean_vector: Error opening: \n%s \n", path.c_str());
        exit(EXIT_FAILURE);
    }
    index_parameters_.ld_mean_vector = new float[index_parameters_.ld_length];
    n_read = fread(index_parameters_.ld_mean_vector, sizeof(float),
                   index_parameters_.ld_length, mv_file);
    fclose(mv_file);
}

void GDIndex::load_ld_pca_eigenvectors(string path) {
    int n_read;
    FILE* eig_file = fopen(path.c_str(), "rb");
    if (eig_file == NULL) {
        printf("GDIndex::load_ld_pca_eigenvectors: Error opening: \n%s \n", path.c_str());
        exit(EXIT_FAILURE);
    }
    index_parameters_.ld_pca_eigenvectors.resize(index_parameters_.ld_length, NULL);
    for (uint n = 0; n < index_parameters_.ld_length; n++) {
        index_parameters_.ld_pca_eigenvectors.at(n) 
            = new float[index_parameters_.ld_length];
        n_read = fread(index_parameters_.ld_pca_eigenvectors.at(n), sizeof(float),
                       index_parameters_.ld_length, eig_file);
    }
    fclose(eig_file);
}

void GDIndex::load_gd_gmm(string path) {
    FILE* gmm_file = fopen(path.c_str(), "rb");
    if (gmm_file == NULL) {
        printf("GDIndex::load_gd_gmm: Error opening: \n%s \n", path.c_str());
        exit(EXIT_FAILURE);
    }
    index_parameters_.gd_gmm = gmm_read(gmm_file);
    fclose(gmm_file);
}

void GDIndex::load_corr_weights(string path) {
    int n_read;
    FILE* cw_file = fopen(path.c_str(), "rb");
    if (cw_file == NULL) {
        printf("GDIndex::load_corr_weights: Error opening: \n%s \n", path.c_str());
        exit(EXIT_FAILURE);
	}
	float* weights = new float[index_parameters_.ld_pca_dim + 1];
	n_read = fread(weights, sizeof(float), index_parameters_.ld_pca_dim + 1, 
                   cw_file);
	fclose(cw_file);
    
    // Build fast_corr_weights
    query_parameters_.fast_corr_weights = new float[index_parameters_.ld_pca_dim + 1];
    // -- initialize them to zero
    for (uint i = 0; i < index_parameters_.ld_pca_dim + 1; i++) {
        query_parameters_.fast_corr_weights[i] = 0;
    }

    for (uint count_bin = 0; count_bin < index_parameters_.ld_pca_dim + 1; count_bin++) { 
        float prob = weights[index_parameters_.ld_pca_dim - count_bin]/2.0; // normalize to 1
        float corr = index_parameters_.ld_pca_dim - 2*count_bin;
        if (count_bin < CORR_WEIGHTS_CLIPPING) {
            query_parameters_.fast_corr_weights[count_bin] = prob * corr;
        }
	}

    // Clean up
    if (weights != NULL) {
        delete [] weights;
        weights = NULL;
    }
}


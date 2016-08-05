#include <algorithm>
#include <cassert>
#include <float.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>

#include "../../common/file_io/file_io.h"
#include "../../common/feature_set/feature_set.h"
#include "gdindex.h"

extern "C" {
#include "../../common/yael_v260_modif/yael/gmm.h"
#include "../../common/yael_v260_modif/yael/matrix.h"
}

/********************************
Constant strings and floats from GDIndex
********************************/
const string GDIndex::SIFT_EXTENSION = ".siftb";
const string GDIndex::SIFTGEO_EXTENSION = ".siftgeo";
const string GDIndex::SIFT_NAME = "sift";
const string GDIndex::SIFTGEO_NAME = "siftgeo";

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
    index_parameters_.ld_pca_dim = LD_PCA_DIM;
    index_parameters_.ld_pre_pca_power = LD_PRE_PCA_POWER_DEFAULT;
    index_parameters_.ld_mean_vector = NULL;
    index_parameters_.ld_pca_eigenvectors.clear();
    index_parameters_.gd_gmm = NULL;
    index_parameters_.gd_number_gaussians = GD_NUMBER_GAUSSIANS_DEFAULT;
    index_parameters_.gd_power = GD_POWER_DEFAULT;
    index_parameters_.gd_intra_normalization = GD_INTRA_NORMALIZATION_DEFAULT;
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

void GDIndex::read_all(const string index_path) {
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

    index_.word_descriptor.resize(current_db_size + number_gd_to_read);
    index_.word_l1_norms.resize(current_db_size + number_gd_to_read);
    index_.word_total_soft_assignment.resize(current_db_size + number_gd_to_read);
    
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
        index_.word_l1_norms.at(count_item)
            .resize(index_parameters_.gd_number_gaussians);
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians; 
             count_gaussian++) {
            index_.word_l1_norms.at(count_item).at(count_gaussian)
                = word_l1_norms_to_read[count_gaussian];
        }
        index_.word_total_soft_assignment.at(count_item)
            .resize(index_parameters_.gd_number_gaussians);
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians; 
             count_gaussian++) {
            index_.word_total_soft_assignment.at(count_item).at(count_gaussian)
                = word_total_soft_assignment_to_read[count_gaussian];
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
    uint number_files_to_process = feature_files.size();

    // Allocate space in index_
    index_.word_l1_norms.resize(number_files_to_process);
    index_.word_total_soft_assignment.resize(number_files_to_process);
    index_.word_descriptor.resize(number_files_to_process);

	uint number_files_processed = 0;
#pragma omp parallel for
    for (uint count_file = 0; count_file < number_files_to_process; count_file++) {
        // Load feature set
        if (!file_exists(feature_files.at(count_file))) {
            fprintf(stderr, "Missing feature file: %s\n", 
                    feature_files.at(count_file).c_str());
            exit(EXIT_FAILURE);
        }
        FeatureSet* feature_set = NULL;
        if (index_parameters_.ld_name == SIFT_NAME) {
            feature_set = readSIFTFile(feature_files.at(count_file), 
                                       index_parameters_.ld_frame_length,
                                       index_parameters_.ld_length);
        } else if (index_parameters_.ld_name == SIFTGEO_NAME) {
            feature_set = readSIFTGeoFile(feature_files.at(count_file), 
                                          index_parameters_.ld_frame_length,
                                          index_parameters_.ld_length);
        } else {
            cout << "Local feature " << index_parameters_.ld_name
                 << " is not supported" << endl;
        }

        // Generate global signature, put in index_
        index_.word_l1_norms.at(count_file).resize(index_parameters_.gd_number_gaussians);
        index_.word_total_soft_assignment.at(count_file).resize(index_parameters_.gd_number_gaussians);
        index_.word_descriptor.at(count_file).resize(index_parameters_.gd_number_gaussians);
        generate_global_descriptor(feature_set, 
                                   index_.word_descriptor.at(count_file),
                                   index_.word_l1_norms.at(count_file),
                                   index_.word_total_soft_assignment.at(count_file));

        // Report status
#pragma omp critical 
        {
            number_files_processed++;
        }
        if (verbose_level >= 2) {
            printf("[%06d, %06d] %04d features \n", 
                   count_file, number_files_processed, feature_set->m_nNumFeatures);
        }

        // Clean up feature set
        if (feature_set == NULL) {
            delete feature_set;
        }
    }
  
    update_index();
}

void GDIndex::generate_index_shot_based(const vector<string>& feature_files, 
                                        const vector<uint>& shot_beg_frames,
                                        const int shot_mode, const int shot_keyf, 
                                        const int verbose_level) {
    if (verbose_level >= 4) cout << "Starting generate_index_shot_based..." 
                                 << endl;
    uint number_shots = shot_beg_frames.size();
    if (verbose_level >= 3) cout << "Using " << number_shots << " shots" << endl;

    uint number_frames_in_db = feature_files.size();
    if (verbose_level >= 3) cout << "There are " << number_frames_in_db 
                                 << " frames in the database." << endl;

    // The number of entries in the generated index will vary depending
    // on mode and shot_keyf
    uint number_entries_in_index = 0;

    // This vector is used only if frames in shot are stored independently.
    // It contains, for each shot number (first entry), a vector with indices
    // which are the places in the index where each frame will be placed.
    vector < vector < uint > > inds_indep_mode;

    if (shot_mode == SHOT_MODE_INDEP_KEYF) {
        // The number of global signatures will vary for this mode,
        // it depends on the number of frames in each shot, so
        // we need to calculate
        for (uint count_shot = 0; count_shot < number_shots; count_shot++) {
            uint number_frames_this_shot;
            // Note: the number calculated in the following lines is a simple
            // subtraction, since the next number is not part of the shot, so 
            // there's no +1 at the end of the calculation
            if (count_shot != number_shots - 1) {
                number_frames_this_shot = shot_beg_frames.at(count_shot + 1)
                    - shot_beg_frames.at(count_shot);
            } else {
                number_frames_this_shot = number_frames_in_db
                    - shot_beg_frames.at(count_shot);
            }

            uint n;
            if (shot_keyf != -1) {
                n = min(number_frames_this_shot, static_cast<uint>(shot_keyf));
            } else {
                n = number_frames_this_shot;
            }

            vector < uint > shot_inds;
            for (uint count_f = 0; count_f < n; count_f++) {
                shot_inds.push_back(number_entries_in_index + count_f);
            }

            number_entries_in_index += n;
            inds_indep_mode.push_back(shot_inds);
        }		
    } else if (shot_mode == SHOT_MODE_SHOT_AGG) {
        // This will generate only a global signature per shot
        number_entries_in_index = number_shots;
    } else {
        cout << "Indexing for shot_mode " << shot_mode 
             << " is not currently implemented. Quitting..." << endl;
        exit(EXIT_FAILURE);
    }

    if (verbose_level >= 2) cout << "GD index will contain " 
                                 << number_entries_in_index << " entries" << endl;

    // Allocate space in index_
    index_.word_l1_norms.resize(number_entries_in_index);
    index_.word_total_soft_assignment.resize(number_entries_in_index);
    index_.word_descriptor.resize(number_entries_in_index);

    // We will loop over shots and aggregate features depending on the mode
    uint count_index = 0;
#pragma omp parallel for       
    for (uint count_shot = 0; count_shot < number_shots; count_shot++) {
        uint number_frames_this_shot;
        // Note: the number calculated in the following lines is a simple
        // subtraction, since the next number is not part of the shot, so 
        // there's no +1 at the end of the calculation
        if (count_shot != number_shots - 1) {
            number_frames_this_shot = shot_beg_frames.at(count_shot + 1)
                - shot_beg_frames.at(count_shot);
        } else {
            number_frames_this_shot = number_frames_in_db
                - shot_beg_frames.at(count_shot);
        }

        if (verbose_level >= 3) cout << "Doing shot " << count_shot 
                                     << " out of " << number_shots << endl;
        uint number_frames_to_use;
        if (shot_keyf != -1) {
            number_frames_to_use = min(number_frames_this_shot, 
                                       static_cast<uint>(shot_keyf));
        } else {
            number_frames_to_use = number_frames_this_shot;
        }

        if (verbose_level >= 3) cout << "This shot will use " 
                                     << number_frames_to_use << " frames" << endl;

        vector<uint> frames_to_use;
        sample_frames_from_shot(number_frames_to_use,
                                shot_beg_frames.at(count_shot),
                                number_frames_this_shot,
                                frames_to_use);

        if (shot_mode == SHOT_MODE_INDEP_KEYF) {
            for (uint count_ind = 0; count_ind < number_frames_to_use;
                 count_ind++) {
                uint frame_this_ind = frames_to_use.at(count_ind);
                uint ind_in_index = inds_indep_mode.at(count_shot).at(count_ind);
                // Load feature set
                if (!file_exists(feature_files.at(frame_this_ind))) {
                    fprintf(stderr, "Missing feature file: %s\n", 
                            feature_files.at(frame_this_ind).c_str());
                    exit(EXIT_FAILURE);
                }
                FeatureSet* feature_set = NULL;
                if (index_parameters_.ld_name == SIFT_NAME) {
                    feature_set = readSIFTFile(feature_files.at(frame_this_ind), 
                                               index_parameters_.ld_frame_length,
                                               index_parameters_.ld_length);
                } else if (index_parameters_.ld_name == SIFTGEO_NAME) {
                    feature_set = readSIFTGeoFile(feature_files.at(frame_this_ind), 
                                                  index_parameters_.ld_frame_length,
                                                  index_parameters_.ld_length);
                } else {
                    cout << "Local feature " << index_parameters_.ld_name
                         << " is not supported" << endl;
                }
                // Generate global signature, put in index_
                index_.word_l1_norms.at(ind_in_index).resize(index_parameters_.gd_number_gaussians);
                index_.word_total_soft_assignment.at(ind_in_index).resize(index_parameters_.gd_number_gaussians);
                index_.word_descriptor.at(ind_in_index).resize(index_parameters_.gd_number_gaussians);
                generate_global_descriptor(feature_set, 
                                           index_.word_descriptor.at(ind_in_index),
                                           index_.word_l1_norms.at(ind_in_index),
                                           index_.word_total_soft_assignment.at(ind_in_index));                
#pragma omp critical 
                {
                    count_index++;
                    index_.frame_numbers_in_db.push_back(frame_this_ind);
                }
    
                // Clean up feature set
                if (feature_set == NULL) {
                    delete feature_set;
                }
            }
        } else if (shot_mode == SHOT_MODE_SHOT_AGG) {
			// Collect all features
            FeatureSet* feature_set = new FeatureSet(index_parameters_.ld_length,
                                                     index_parameters_.ld_frame_length);
            for (uint count_ind = 0; count_ind < number_frames_to_use; count_ind++) {
                uint frame_this_ind = frames_to_use.at(count_ind);
                if (verbose_level >= 4) cout << "Doing frame " << frame_this_ind << endl;

                // Load feature set for this frame
                FeatureSet* feature_set_this_frame = NULL;
                if (index_parameters_.ld_name == SIFT_NAME) {
                    feature_set_this_frame = readSIFTFile(feature_files.at(frame_this_ind), 
                                                          index_parameters_.ld_frame_length,
                                                          index_parameters_.ld_length);
                } else if (index_parameters_.ld_name == SIFTGEO_NAME) {
                  feature_set_this_frame = readSIFTGeoFile(feature_files.at(frame_this_ind), 
                                                           index_parameters_.ld_frame_length,
                                                           index_parameters_.ld_length);
                } else {
                    cout << "Local feature " << index_parameters_.ld_name
                         << " is not supported" << endl;
                }

                if (verbose_level >= 4) cout << "Loaded features from this frame" 
                                             << endl;

                // Add this frame's features to collection of all shot's features
                uint number_features_this_frame = 
                    feature_set_this_frame->m_nNumFeatures;
                for (uint count_f = 0; count_f < number_features_this_frame; 
                     count_f++) {
                    feature_set->addFeature(feature_set_this_frame->m_vDescriptors.at(count_f),
                                            feature_set_this_frame->m_vFrames.at(count_f));
                }
                if (verbose_level >= 4) cout << "Added these features to shot's features" << endl;
                // Clean up feature set for this frame
                if (feature_set_this_frame == NULL) {
                    delete feature_set_this_frame;
                }
            }

            if (verbose_level >= 4) cout << "All features were collected, now shot contains " << feature_set->m_nNumFeatures << " features" << endl;

            // Generate global signature, put in index_
            index_.word_l1_norms.at(count_shot).resize(index_parameters_.gd_number_gaussians);
            index_.word_total_soft_assignment.at(count_shot).resize(index_parameters_.gd_number_gaussians);
            index_.word_descriptor.at(count_shot).resize(index_parameters_.gd_number_gaussians);
            generate_global_descriptor(feature_set, 
                                       index_.word_descriptor.at(count_shot),
                                       index_.word_l1_norms.at(count_shot),
                                       index_.word_total_soft_assignment.at(count_shot));
#pragma omp critical 
            {
                count_index++;
            }
    
            // Clean up feature set
            if (feature_set == NULL) {
                delete feature_set;
            }
        } 
    }

    // We need to sort frame_numbers_in_db, because it might not be in order 
    // due to parallelization. 
    // Note: this variable is used only if the mode is SHOT_MODE_INDEP_KEYF,
    // so if the mode is different it will just sort nothing, it doesnt matter
    stable_sort(index_.frame_numbers_in_db.begin(), index_.frame_numbers_in_db.end());

    // Make sure we processed the correct number of REVV signatures
    assert(count_index == number_entries_in_index);
    if (shot_mode == SHOT_MODE_INDEP_KEYF) {
        assert(count_index == index_.frame_numbers_in_db.size());
    }

    update_index();
}

void GDIndex::generate_global_descriptor(const FeatureSet* feature_set, 
                                         vector<uint>& gd_word_descriptor, 
                                         vector<float>& gd_word_l1_norm, 
                                         vector<float>& gd_word_total_soft_assignment) {
    // Resize the vectors that will be returned
    gd_word_descriptor.resize(index_parameters_.gd_number_gaussians);
    gd_word_l1_norm.resize(index_parameters_.gd_number_gaussians);
    gd_word_total_soft_assignment.resize(index_parameters_.gd_number_gaussians);

    uint unbinarized_signature_length = index_parameters_.gd_number_gaussians
        *index_parameters_.ld_pca_dim;
    float* all_pca_desc = new float[feature_set->m_nNumFeatures 
                                    * index_parameters_.ld_pca_dim];
    // Project SIFT using PCA
    for (uint count_feat = 0; count_feat < feature_set->m_nNumFeatures; count_feat++) {
        project_local_descriptor_pca(feature_set->m_vDescriptors[count_feat], 
                                     all_pca_desc + count_feat*index_parameters_.ld_pca_dim);
    }

    // Compute Fisher vector
    int gmm_flags = GMM_FLAGS_MU;
    uint fisher_output_length = gmm_fisher_sizeof(index_parameters_.gd_gmm, 
                                                 gmm_flags);
    vector<float> fisher_output(fisher_output_length, 0);
    // This will extract the FV with only the mean component and using FV normalization
    // but NOT the L2 normalization after the end nor the power normalization
    gmm_fisher_save_soft_assgn(feature_set->m_nNumFeatures, all_pca_desc, 
                               index_parameters_.gd_gmm, gmm_flags, fisher_output.data(), 
                               gd_word_total_soft_assignment.data());

    // Compute L1Norm info
    for (uint count_gaussian = 0; 
         count_gaussian < index_parameters_.gd_number_gaussians; 
         count_gaussian++) {
        double sum_abs = 0;
        for (uint count_dim = count_gaussian*index_parameters_.ld_pca_dim; 
             count_dim < (count_gaussian + 1)*index_parameters_.ld_pca_dim; 
             count_dim++) {
            sum_abs += fabs(fisher_output.at(count_dim));
        }
        gd_word_l1_norm.at(count_gaussian) = static_cast<float>(sum_abs);
    }

    // IN normalization (if selected)
    if (index_parameters_.gd_intra_normalization) {
        for (uint count_gaussian = 0; 
             count_gaussian < index_parameters_.gd_number_gaussians; 
             count_gaussian++) {
            // Compute L2 norm for this Gaussian
            float l2_norm_sq_gaussian = 0;
            for (uint count_dim = count_gaussian*index_parameters_.ld_pca_dim; 
                 count_dim < (count_gaussian + 1)*index_parameters_.ld_pca_dim; 
                 count_dim++) {
                l2_norm_sq_gaussian += fisher_output.at(count_dim) * fisher_output.at(count_dim);
            }
            // Normalize this Gaussian, if it has non-zero norm
            float l2_norm_gaussian = sqrt(l2_norm_sq_gaussian);            
            if (l2_norm_gaussian > 0) {
                for (uint count_dim = count_gaussian*index_parameters_.ld_pca_dim; 
                     count_dim < (count_gaussian + 1)*index_parameters_.ld_pca_dim; 
                     count_dim++) {
                    fisher_output.at(count_dim) /= l2_norm_gaussian;
                }
            }
        }
    }
    
    // Apply power law (if using SSR normalization), and compute L2 norm
    float l2_norm_sq = 0;
    for (uint count_dim = 0; count_dim < unbinarized_signature_length; count_dim++) {
        if (!index_parameters_.gd_intra_normalization) {
            POWER_LAW_SAME(fisher_output.at(count_dim), index_parameters_.gd_power);
        }
        l2_norm_sq += fisher_output.at(count_dim) * fisher_output.at(count_dim);
    }

    // L2 normalize
    float l2_norm = sqrt(l2_norm_sq);
    if (l2_norm > L2_NORM_SQ_THRESH) {
        for (uint count_dim = 0; count_dim < unbinarized_signature_length; count_dim++) {
            fisher_output.at(count_dim) /= l2_norm;
        }
    }
    
    // Sign binarize
    sign_binarize(fisher_output, gd_word_descriptor);
}

void GDIndex::perform_query(const string local_descriptors_path, 
                            const GDIndex* query_index_ptr,
                            const uint query_number,
                            const vector<uint>& indices,
                            vector< pair<float,uint> >& results, 
                            const uint number_2nd_stage_rerank,
                            GDIndex* gdindex_ptr_rerank,
                            const vector < vector < uint > >& group_lists_rerank,
                            const int verbose_level) {
    FeatureSet* feature_set = NULL;
    vector<uint> gd_word_descriptor;
    vector<float> gd_word_l1_norm, gd_word_total_soft_assignment;
    if (query_index_ptr != NULL) {
        // Using pre-computed global descriptor from query_index_ptr
        gd_word_descriptor = 
            query_index_ptr->index_.word_descriptor.at(query_number);
        if (query_parameters_.word_selection_mode == WORD_L1_NORM) {
            gd_word_l1_norm = 
                query_index_ptr->index_.word_l1_norms.at(query_number);
        } else {
            gd_word_total_soft_assignment = 
                query_index_ptr->index_.word_total_soft_assignment.at(query_number);
        }
    } else {
        // Computing query global descriptor

        // --> Load local descriptors
        if (index_parameters_.ld_name == SIFT_NAME) {
            feature_set = readSIFTFile(local_descriptors_path, 
                                       index_parameters_.ld_frame_length,
                                       index_parameters_.ld_length);
        } else if (index_parameters_.ld_name == SIFTGEO_NAME) {
            feature_set = readSIFTGeoFile(local_descriptors_path, 
                                          index_parameters_.ld_frame_length,
                                          index_parameters_.ld_length);
        } else {
            cout << "Local feature " << index_parameters_.ld_name
                 << " is not supported" << endl;
        }

        // --> Generate query global descriptor
        generate_global_descriptor(feature_set, 
                                   gd_word_descriptor, 
                                   gd_word_l1_norm, 
                                   gd_word_total_soft_assignment);
    }
    
    // If number_2nd_stage_rerank is 0, we're not using two-stage scoring
    if (!number_2nd_stage_rerank) {
        query(gd_word_descriptor, gd_word_l1_norm, gd_word_total_soft_assignment,
              indices, results);
    } else {
        // Using two-stage scoring
        // -- First stage
        vector<uint> indices_2_stages(index_.number_global_descriptors);
        for (uint count = 0; count < index_.number_global_descriptors; count++) {
            // Include indices corresponding to shot numbers
            indices_2_stages.at(count) = count;
        }
        vector< pair<float,uint> > results_1st_stage;
        query(gd_word_descriptor, gd_word_l1_norm, gd_word_total_soft_assignment,
              indices_2_stages, results_1st_stage);

        // -- Second stage
        vector<uint> gd_word_descriptor_rerank;
        vector<float> gd_word_l1_norm_rerank, gd_word_total_soft_assignment_rerank;
        gdindex_ptr_rerank->generate_global_descriptor(feature_set,
                                                       gd_word_descriptor_rerank,
                                                       gd_word_l1_norm_rerank,
                                                       gd_word_total_soft_assignment_rerank);
        gdindex_ptr_rerank->query_2nd_stage(gd_word_descriptor_rerank, 
                                            gd_word_l1_norm_rerank, 
                                            gd_word_total_soft_assignment_rerank, 
                                            number_2nd_stage_rerank,
                                            group_lists_rerank,
                                            results_1st_stage, results);
    }

    // Clean up
    if (feature_set != NULL) {
        delete feature_set;
        feature_set = NULL;
    }
}

void GDIndex::set_index_parameters(const uint ld_length, const uint ld_frame_length,
                                   const string ld_extension, const string ld_name,
                                   const uint ld_pca_dim, const float ld_pre_pca_power,
                                   const uint gd_number_gaussians, const float gd_power,
                                   const bool gd_intra_normalization,
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
    index_parameters_.gd_intra_normalization = gd_intra_normalization;

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

void GDIndex::score_database_item(const vector<uint>& query_word_descriptor,
                                  const vector<float>& query_word_l1_norm,
                                  const vector<float>& query_word_total_soft_assignment,
                                  const uint db_ind,
                                  float& score) {
    // Check that the database item has at least the minimum number of
    // selected words that we require; if not, just set it to FLT_MAX
    if (index_.number_words_selected.at(db_ind) 
        <= query_parameters_.min_number_words_selected) {
        score = FLT_MAX;
        return;
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
                ^ index_.word_descriptor.at(db_ind).at(count_gaussian);
            uint h = query_parameters_.pop_count[aux >> 16]
                + query_parameters_.pop_count[aux & 65535];
            // Add to correlation
            total_correlation += query_parameters_.fast_corr_weights[h];
        }
    }

    // Figure out query norm to be used
    query_norm = sqrt(query_number_words_selected
                      * index_parameters_.ld_pca_dim);

    float corr_den = index_.norm_factors.at(db_ind) * query_norm;

    // Compute final correlation for this item
    if (corr_den != 0) {
        total_correlation /= corr_den;
    } else {
        total_correlation = -FLT_MAX;
    } 

    // Change sign such that smaller is better
    score = -total_correlation;
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

        score_database_item(query_word_descriptor,
                            query_word_l1_norm,
                            query_word_total_soft_assignment,
                            count_elem,
                            database_scores_indices.at(count_elem).first);
    }
    // Sort scores
    sort(database_scores_indices.begin(), database_scores_indices.end(), 
         cmp_float_uint_ascend);
}

void GDIndex::query_2nd_stage(const vector<uint>& query_word_descriptor,
                              const vector<float>& query_word_l1_norm,
                              const vector<float>& query_word_total_soft_assignment,
                              const uint number_2nd_stage_rerank,
                              const vector < vector < uint > >& group_lists_rerank,
                              const vector< pair<float,uint> >& first_stage_scores_indices,
                              vector< pair<float,uint> >& database_scores_indices) {
    // Resize vector that will be passed back
    database_scores_indices.clear();

    // Figure out number of first stage items (groups) to re-rank
    uint number_rerank = min(number_2nd_stage_rerank,
                             static_cast<uint>(first_stage_scores_indices.size()));

    // Loop over top number_rerank items from first stage and get score
    // of the constituent second stage items
    for (uint count_top = 0; count_top < number_rerank;
         count_top++) {
        uint this_group_ind = first_stage_scores_indices.at(count_top).second;

		// Get number of 2nd stage items in this group
		uint number_items_this_group = group_lists_rerank.at(this_group_ind).size();

		// Loop over 2nd stage items in this group, get their signatures and score them
		for (uint count_item = 0; count_item < number_items_this_group;
		     count_item++) {
			uint this_item_number = group_lists_rerank.at(this_group_ind).at(count_item);
			pair < float, uint > score_this_item;
            score_this_item.second = this_item_number;
            
            score_database_item(query_word_descriptor,
                                query_word_l1_norm,
                                query_word_total_soft_assignment,
                                this_item_number,
                                score_this_item.first);
            database_scores_indices.push_back(score_this_item);
        }
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


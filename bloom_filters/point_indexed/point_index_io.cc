#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>

#include "point_index_io.h"

const int LARGE_NUM_FEATS = 10000;

void write_fv_point_index(const vector < vector < uint > >& vec_feat_assgns,
                          const vector < vector < float > >& vec_feat_assgn_weights,
                          const vector < vector < vector < float > > >& vec_feat_residuals,
                          const int residual_length,
                          const string file_path) {
    int num_items = vec_feat_assgns.size();

    // Make sure data makes sense
    assert(static_cast<size_t>(num_items) == vec_feat_assgn_weights.size());
    assert(static_cast<size_t>(num_items) == vec_feat_residuals.size());

    // Open file
    FILE* p_index_file = fopen(file_path.c_str(), "wb");
    if (p_index_file == NULL) {
        fprintf(stderr, "write_fv_point_index: Cannot open: %s\n", file_path.c_str());
        exit(EXIT_FAILURE);
    }

    // Write number of items
    int n_write = fwrite(&num_items, sizeof(int), 1, p_index_file);    

    // Write residual length
    n_write = fwrite(&residual_length, sizeof(int), 1, p_index_file);    

    // Loop and write data for each item
    for (int i = 0; i < num_items; i++) {
        int num_feats_this_item = vec_feat_assgns.at(i).size();

        // Write number of features for this item
        n_write = fwrite(&num_feats_this_item, sizeof(int), 1, p_index_file);    

        // Make sure data for this item makes sense
        assert(static_cast<size_t>(num_feats_this_item) == vec_feat_assgn_weights.at(i).size());
        assert(static_cast<size_t>(num_feats_this_item) == vec_feat_residuals.at(i).size());

        // Write feat_assgns for this item
        n_write = fwrite(vec_feat_assgns.at(i).data(), sizeof(uint), 
                         num_feats_this_item, p_index_file);
        
        // Write feat_assgn_weights for this item
        n_write = fwrite(vec_feat_assgn_weights.at(i).data(), sizeof(float),
                         num_feats_this_item, p_index_file);

        // Write feat_residuals for this item
        for (int j = 0; j < num_feats_this_item; j++) {
            assert(static_cast<size_t>(residual_length) == vec_feat_residuals.at(i).at(j).size());
            n_write = fwrite(vec_feat_residuals.at(i).at(j).data(), sizeof(float),
                             residual_length, p_index_file);            
        }
    }

    // Close file
    fclose(p_index_file);
}

void read_fv_point_index(const string file_path,
                         vector < vector < uint > >& vec_feat_assgns,
                         vector < vector < float > >& vec_feat_assgn_weights,
                         vector < vector < vector < float > > >& vec_feat_residuals,
                         const int max_feats_to_read) {
    // Clear vectors that will be returned
    vec_feat_assgns.clear();
    vec_feat_assgn_weights.clear();
    vec_feat_residuals.clear();

    // Open file
    FILE* p_index_file = fopen(file_path.c_str(), "rb");
    if (p_index_file == NULL) {
        fprintf(stderr, "read_fv_point_index: Cannot open: %s\n", file_path.c_str());
        exit(EXIT_FAILURE);
    }

    // Read number of items and resize vectors that will be output
    int num_items = 0;    
    int n_read = fread(&num_items, sizeof(int), 1, p_index_file);
    vec_feat_assgns.resize(num_items);
    vec_feat_assgn_weights.resize(num_items);
    vec_feat_residuals.resize(num_items);

    // Read residual length
    int residual_length = 0;    
    n_read = fread(&residual_length, sizeof(int), 1, p_index_file);

    // Declare aux structures that will read data
    vector<uint> vec_feat_assgns_aux(LARGE_NUM_FEATS);
    vector<float> vec_feat_assgn_weights_aux(LARGE_NUM_FEATS);
    vector<float> vec_feat_residuals_aux(residual_length);

    // Loop and read data for each item
    for (int i = 0; i < num_items; i++) {
        // Read number of features in this item
        int num_feats_this_item = 0;    
        n_read = fread(&num_feats_this_item, sizeof(int), 1, p_index_file);

        // Decide how many features will be effectively used
        int num_feats_to_use = 0;    
        if (max_feats_to_read != -1) {
          num_feats_to_use = min(max_feats_to_read, num_feats_this_item);
        } else {
          num_feats_to_use = num_feats_this_item;
        }

        // Potentially resize reading vectors if necessary
        if (num_feats_this_item > LARGE_NUM_FEATS) {
          vec_feat_assgns_aux.resize(num_feats_this_item);
          vec_feat_assgn_weights_aux.resize(num_feats_this_item);
        }

        // Read feat_assgns for this item
        vec_feat_assgns.at(i).resize(num_feats_to_use);
        n_read = fread(vec_feat_assgns_aux.data(), sizeof(uint), 
                       num_feats_this_item, p_index_file);
        copy(vec_feat_assgns_aux.begin(), 
             vec_feat_assgns_aux.begin() + num_feats_to_use, 
             vec_feat_assgns.at(i).begin());

        // Read feat_assgn_weights for this item
        vec_feat_assgn_weights.at(i).resize(num_feats_to_use);
        n_read = fread(vec_feat_assgn_weights_aux.data(), sizeof(float), 
                       num_feats_this_item, p_index_file);
        copy(vec_feat_assgn_weights_aux.begin(), 
             vec_feat_assgn_weights_aux.begin() + num_feats_to_use,
             vec_feat_assgn_weights.at(i).begin());

        // Read feat_residuals for this item
        vec_feat_residuals.at(i).resize(num_feats_to_use);
        for (int j = 0; j < num_feats_this_item; j++) {
            n_read = fread(vec_feat_residuals_aux.data(), sizeof(float),
                           residual_length, p_index_file);
            if (j < num_feats_to_use) {
              vec_feat_residuals.at(i).at(j).resize(residual_length);
              copy(vec_feat_residuals_aux.begin(), 
                   vec_feat_residuals_aux.end(),
                   vec_feat_residuals.at(i).at(j).begin());              
            }
        }
    }    
    fclose(p_index_file);
}

void write_bfv_point_index(const vector < vector < uint > >& vec_feat_assgns,
                           const vector < vector < float > >& vec_feat_assgn_weights,
                           const vector < vector < uint > >& vec_feat_residuals_binarized,
                           const string file_path) {
    int num_items = vec_feat_assgns.size();

    // Make sure data makes sense
    assert(static_cast<size_t>(num_items) == vec_feat_assgn_weights.size());
    assert(static_cast<size_t>(num_items) == vec_feat_residuals_binarized.size());

    // Open file
    FILE* p_index_file = fopen(file_path.c_str(), "wb");
    if (p_index_file == NULL) {
        fprintf(stderr, "write_bfv_point_index: Cannot open: %s\n", file_path.c_str());
        exit(EXIT_FAILURE);
    }

    // Write number of items
    int n_write = fwrite(&num_items, sizeof(int), 1, p_index_file);    

    // Loop and write data for each item
    for (int i = 0; i < num_items; i++) {
        int num_feats_this_item = vec_feat_assgns.at(i).size();

        // Make sure data for this item makes sense
        assert(static_cast<size_t>(num_feats_this_item) == vec_feat_assgn_weights.at(i).size());
        assert(static_cast<size_t>(num_feats_this_item) == vec_feat_residuals_binarized.at(i).size());

        // Write number of features for this item
        n_write = fwrite(&num_feats_this_item, sizeof(int), 1, p_index_file);    

        // Write feat_assgns for this item
        n_write = fwrite(vec_feat_assgns.at(i).data(), sizeof(uint), 
                         num_feats_this_item, p_index_file);
        
        // Write feat_assgn_weights for this item
        n_write = fwrite(vec_feat_assgn_weights.at(i).data(), sizeof(float),
                         num_feats_this_item, p_index_file);

        // Write feat_residuals_binarized for this item
        n_write = fwrite(vec_feat_residuals_binarized.at(i).data(), sizeof(uint),
                         num_feats_this_item, p_index_file);

    }

    // Close file
    fclose(p_index_file);
}

void read_bfv_point_index(const string file_path,
                          vector < vector < uint > >& vec_feat_assgns,
                          vector < vector < float > >& vec_feat_assgn_weights,
                          vector < vector < uint > >& vec_feat_residuals_binarized,
                          const int max_feats_to_read) {
    // Clear vectors that will be returned
    vec_feat_assgns.clear();
    vec_feat_assgn_weights.clear();
    vec_feat_residuals_binarized.clear();

    // Open file
    FILE* p_index_file = fopen(file_path.c_str(), "rb");
    if (p_index_file == NULL) {
      fprintf(stderr, "read_bfv_point_index: Cannot open: %s\n", file_path.c_str());
      exit(EXIT_FAILURE);
    }

    // Read number of items and resize vectors that will be output
    int num_items = 0;    
    int n_read = fread(&num_items, sizeof(int), 1, p_index_file);

    vec_feat_assgns.resize(num_items);
    vec_feat_assgn_weights.resize(num_items);
    vec_feat_residuals_binarized.resize(num_items);

    // Declare aux structures that will read data
    vector<uint> vec_feat_assgns_aux(LARGE_NUM_FEATS);
    vector<float> vec_feat_assgn_weights_aux(LARGE_NUM_FEATS);
    vector<uint> vec_feat_residuals_binarized_aux(LARGE_NUM_FEATS);

    // Loop and read data for each item
    for (int i = 0; i < num_items; i++) {
        // Read number of features in this item
        int num_feats_this_item = 0;    
        n_read = fread(&num_feats_this_item, sizeof(int), 1, p_index_file);

        // Decide how many features will be effectively used
        int num_feats_to_use = 0;    
        if (max_feats_to_read != -1) {
          num_feats_to_use = min(max_feats_to_read, num_feats_this_item);
        } else {
          num_feats_to_use = num_feats_this_item;
        }

        // Potentially resize reading vectors if necessary
        if (num_feats_this_item > LARGE_NUM_FEATS) {
          vec_feat_assgns_aux.resize(num_feats_this_item);
          vec_feat_assgn_weights_aux.resize(num_feats_this_item);
          vec_feat_residuals_binarized_aux.resize(num_feats_this_item);          
        }

        // Read feat_assgns for this item
        vec_feat_assgns.at(i).resize(num_feats_to_use);
        n_read = fread(vec_feat_assgns_aux.data(), sizeof(uint), 
                       num_feats_this_item, p_index_file);
        copy(vec_feat_assgns_aux.begin(), 
             vec_feat_assgns_aux.begin() + num_feats_to_use, 
             vec_feat_assgns.at(i).begin());

        // Read feat_assgn_weights for this item
        vec_feat_assgn_weights.at(i).resize(num_feats_to_use);
        n_read = fread(vec_feat_assgn_weights_aux.data(), sizeof(float), 
                       num_feats_this_item, p_index_file);
        copy(vec_feat_assgn_weights_aux.begin(), 
             vec_feat_assgn_weights_aux.begin() + num_feats_to_use,
             vec_feat_assgn_weights.at(i).begin());

        // Read feat_residuals_binarized for this item
        vec_feat_residuals_binarized.at(i).resize(num_feats_to_use);
        n_read = fread(vec_feat_residuals_binarized_aux.data(), sizeof(uint), 
                       num_feats_this_item, p_index_file);
        copy(vec_feat_residuals_binarized_aux.begin(), 
             vec_feat_residuals_binarized_aux.begin() + num_feats_to_use,
             vec_feat_residuals_binarized.at(i).begin());
    }    

    fclose(p_index_file);
}


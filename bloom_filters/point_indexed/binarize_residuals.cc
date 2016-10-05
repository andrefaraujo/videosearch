#include <cassert>
#include <vector>

#include "binarize_residuals.h"

// This is hard-coded, as we assume throughout this codebase that residuals
// contain 32 floating point dimensions; this allows us to perform
// simple binarization that fits within 4-byte ints
const uint RESIDUAL_LENGTH = 32;

#ifndef PUSH_BIT
#define PUSH_BIT(packed, bit) \
  packed = packed << 1; \
  packed += bit;
#endif

// Function to binarize only one item
void binarize_residual(const vector < vector < float > >& feat_residuals,
                       vector<uint>& feat_residuals_binarized) {
    uint num_features = feat_residuals.size();
    feat_residuals_binarized.resize(num_features);
    for (uint count_feat = 0; count_feat < num_features; count_feat++) {
        assert(feat_residuals.at(count_feat).size()
               == RESIDUAL_LENGTH);
        uint packed_block = 0;
        for (uint count_dim = 0; count_dim < RESIDUAL_LENGTH; count_dim++) {
            uint bit_uint = 
                (feat_residuals.at(count_feat).at(count_dim) > 0) 
                ? 1 : 0;
            PUSH_BIT(packed_block, bit_uint);
        }
        feat_residuals_binarized.at(count_feat) = packed_block;
    }
}

// Function to binarize several items
void binarize_residuals(const vector < vector < vector < float > > >& vec_feat_residuals, 
                       vector < vector < uint > >& vec_feat_residuals_binarized) {
    uint num_items = vec_feat_residuals.size();    
    vec_feat_residuals_binarized.resize(num_items);

    for (uint count_item = 0; count_item < num_items; count_item++) {
        binarize_residual(vec_feat_residuals.at(count_item),
                          vec_feat_residuals_binarized.at(count_item));
    }
}

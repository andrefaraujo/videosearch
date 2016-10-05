#include <algorithm>
#include <vector>

#include "hasher.h"

using namespace std;

const uint BITS_IN_BYTE = 8;

VBHHasher::VBHHasher(size_t num_bits, size_t item_length) {
    vector<size_t> indices;
    for (size_t i = 0; i < item_length; i++) 
        indices.push_back(i);
    // Shuffle indices (to avoid getting repeated ones)
    random_shuffle(indices.begin(), indices.end());
    for (size_t i = 0; i < num_bits; i++) {
        indices_to_sample_.push_back(indices.at(i));
    }
    
    // Number of bits in uint, used during hashing
    bits_in_uint_ = sizeof(uint) * BITS_IN_BYTE;
}

size_t VBHHasher::hash(const vector<uint>& item) {
    size_t output = 0;
    for (size_t i = 0; i < indices_to_sample_.size(); i++) {
        size_t idx = indices_to_sample_.at(i);
        uint val = item.at(idx / bits_in_uint_);
        output = (output << 1) | ( (val >> (idx%bits_in_uint_) ) & 1 );
    }
    return output;
}

GBHHasher::GBHHasher(size_t num_bits, size_t residual_length,
                     size_t first_residual_ind) {
    vector<size_t> indices;
    for (size_t i = 0; i < residual_length; i++) 
        indices.push_back(i);
    // Shuffle indices (to avoid getting repeated ones)
    random_shuffle(indices.begin(), indices.end());
    for (size_t i = 0; i < num_bits; i++) {
        indices_to_sample_.push_back(first_residual_ind + indices.at(i));
    }
    first_ind_ = first_residual_ind;
    
    // Number of bits in uint, used during hashing
    bits_in_uint_ = sizeof(uint) * BITS_IN_BYTE;
}

size_t GBHHasher::hash(const vector<uint>& item) {
    return hash(item, first_ind_/bits_in_uint_);
}

size_t GBHHasher::hash(const vector<uint>& item, uint ind) {
    size_t output = 0;
    uint val = item.at(ind);
    for (size_t i = 0; i < indices_to_sample_.size(); i++)
        output = (output << 1) | ( (val >> (indices_to_sample_.at(i) - first_ind_)) & 1 );

    return output;
}

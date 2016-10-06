#include <algorithm>
#include <vector>

#include "hasher.h"

using namespace std;

const uint BITS_IN_BYTE = 8;

GBHHasher::GBHHasher(size_t num_bits, size_t residual_length) {
    vector<size_t> indices;
    for (size_t i = 0; i < residual_length; i++) {
        indices.push_back(i);
    }
    // Shuffle indices (to avoid getting repeated ones)
    random_shuffle(indices.begin(), indices.end());
    for (size_t i = 0; i < num_bits; i++) {
        indices_to_sample_.push_back(indices.at(i));
    }
    
    // Number of bits in uint, used during hashing
    bits_in_uint_ = sizeof(uint) * BITS_IN_BYTE;
}

uint GBHHasher::hash(const uint item) {
    uint output = 0;
    for (size_t i = 0; i < indices_to_sample_.size(); i++)
        output = (output << 1) | ( (item >> indices_to_sample_.at(i)) & 1 );
    return output;
}

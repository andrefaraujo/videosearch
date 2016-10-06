#ifndef HASHER_H
#define HASHER_H

#include <vector>

using namespace std;

typedef unsigned int uint;

class GBHHasher {
 public:
    GBHHasher(size_t num_bits, size_t residual_length);
    uint hash(const uint item);

 private:
    vector<size_t> indices_to_sample_;
    uint bits_in_uint_;
};

#endif

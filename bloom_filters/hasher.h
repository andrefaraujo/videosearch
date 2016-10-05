#ifndef HASHER_H
#define HASHER_H

#include <vector>

using namespace std;

typedef unsigned int uint;

template <class T>
class HashFn {
 public:
    virtual ~HashFn() {}
    virtual size_t hash(const T& item) = 0;
};

class VBHHasher: public HashFn< vector<uint> > {
 public:
    VBHHasher(size_t num_bits, size_t item_length);
    size_t hash(const vector<uint>& item);

 private:
    vector<size_t> indices_to_sample_;
    uint bits_in_uint_;
};

class GBHHasher: public HashFn< vector<uint> > {
 public:
    GBHHasher(size_t num_bits, size_t residual_length,
              size_t first_residual_ind);
    size_t hash(const vector<uint>& item);
    
 private:
    size_t hash(const vector<uint>& item, uint ind);

    vector<size_t> indices_to_sample_;
    size_t first_ind_;
    uint bits_in_uint_;
};

#endif

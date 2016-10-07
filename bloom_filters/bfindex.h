/*
Index of Bloom filters, using point-indexed features
*/

#ifndef BFINDEX_H
#define BFINDEX_H

#include <vector>

#include "bloom.h"
#include "hasher.h"
#include "point_indexed/point_index_io.h"

using namespace std;

typedef unsigned int uint;

class BFIndex {
 public:
    // Constructor/Destructor
    BFIndex();
    ~BFIndex();

    // Initialization function, should be called after constructor
    void initialize(size_t n_bits, size_t n_hashers,
                    size_t residual_len = 32, float alpha = 1,
                    int verbose = 0);

    // BF construction for each video clip.
    // This must be called before perform_query()
    void insert_from_indexes(const vector<string>& index_paths);
    
    // Function to perform querying of an image against the BF database
    void perform_query(const vector<uint>& query_emb_features,
                       const vector<uint>& query_hash_numbers,
                       vector< pair<float,uint> >& results);

 private:
    InvertedIndexBloom* inverted_index_bloom_ptr_;
    size_t num_bloom_filters_;
    size_t residual_length_;
    vector<GBHHasher*> hash_functions_;
    size_t num_hashers_;
    size_t num_bits_;
    int verbose_;
    float alpha_;
    vector < vector < float > > idfs_;
    vector < float > idf_norms_;

    // Function to get number of items/indices
    void get_number_indices(size_t& number_indices);

    // Function to pre-process idf weights
    void idf_processing();
};

#endif

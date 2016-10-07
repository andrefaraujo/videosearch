#include <cmath>
#include <vector>

#include "bloom.h"

using namespace std;

InvertedIndexBloom::InvertedIndexBloom(const size_t num_bits,
                                       const size_t num_hashers)
    : num_bits_(num_bits), num_hashers_(num_hashers), num_items_(0) {
    size_t max_num_buckets = pow(2, num_bits);
    partitioned_index_.resize(num_hashers_);
    for (size_t i = 0; i < num_hashers_; i++) {
        partitioned_index_.at(i).resize(max_num_buckets);
    }
}

InvertedIndexBloom::~InvertedIndexBloom() {
    partitioned_index_.clear();
}

void InvertedIndexBloom::insert(const uint hash_ind, const uint hash_number,
                                const uint item_number) {
}

void InvertedIndexBloom::get_properties(size_t& number_items,
                                        size_t& number_indices) {
}

void InvertedIndexBloom::get_counts_per_hash(const uint hash_number,
                                             vector<uint>& counts) {
}

void InvertedIndexBloom::get_idf_norm(const vector< vector<float> >& idfs,
                                      const float alpha,
                                      const size_t num_bloom_filters,
                                      vector<float>& idf_norms) {
}

void InvertedIndexBloom::perform_query(const vector<uint>& query_hash_indices,
                                       const vector<uint>& query_hash_numbers,
                                       const vector< vector<float> >& idfs,
                                       const vector<float>& idf_norms,
                                       vector< pair<float,uint> >& results) {
}

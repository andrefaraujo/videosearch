#include <cmath>
#include <vector>

#include "bloom.h"

using namespace std;

InvertedIndexBloom::InvertedIndexBloom(const size_t num_bits,
                                       const size_t num_hashers)
    : num_bits_(num_bits), num_hashers_(num_hashers) {
    size_t max_num_buckets = pow(2, num_bits_);
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
    if (!binary_search(partitioned_index_.at(hash_number).at(hash_ind).begin(),
                       partitioned_index_.at(hash_number).at(hash_ind).end(),
                       item_number)) {
        // Find position where to insert
        vector<uint>::iterator position_to_insert =
            partitioned_index_.at(hash_number).at(hash_ind).begin();
        if (position_to_insert != partitioned_index_.at(hash_number)
            .at(hash_ind).end()) {
            while (*position_to_insert < item_number) {
                position_to_insert++;
                if (position_to_insert ==
                    partitioned_index_.at(hash_number).at(hash_ind).end()) {
                    break;
                }
            }
        }
        // Insert element
        partitioned_index_.at(hash_number).at(hash_ind)
            .insert(position_to_insert, item_number);
    }
}

void InvertedIndexBloom::get_number_indices(size_t& number_indices) {
    number_indices = 0;
    for (size_t i = 0; i < num_hashers_; i++) {
        for (size_t j = 0; j < partitioned_index_.at(i).size(); j++) {
            number_indices += partitioned_index_.at(i).at(j).size();
        }
    }
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

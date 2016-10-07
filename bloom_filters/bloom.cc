#include <algorithm>
#include <cassert>
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
    for (uint count_ind = 0;
         count_ind < partitioned_index_.at(hash_number).size();
         count_ind++) {
        counts.at(count_ind) = partitioned_index_.at(hash_number)
            .at(count_ind).size();
    }
}

void InvertedIndexBloom::get_idf_norm(const vector< vector<float> >& idfs,
                                      const float alpha,
                                      const size_t num_bloom_filters,
                                      vector<float>& idf_norms) {
    // Initialize norms to 0
    for (uint count_bf = 0; count_bf < num_bloom_filters; count_bf++) {
        idf_norms.at(count_bf) = 0;
    }

    // Loop over index and construct partial BFs, so that we'll be able to
    // get their IDF norms
    for (uint count_hasher = 0; count_hasher < num_hashers_;
         count_hasher++) {
        vector < vector < uint > >
            bf_partitioned_indices_per_hash(num_bloom_filters);
        for (uint count_ind = 0;
             count_ind < partitioned_index_.at(count_hasher).size();
             count_ind++) {
            size_t list_size = partitioned_index_.at(count_hasher).at(count_ind).size();
            for (size_t count_l = 0; count_l < list_size; count_l++) {
                bf_partitioned_indices_per_hash
                    .at(partitioned_index_.at(count_hasher).at(count_ind).at(count_l))
                    .push_back(count_ind);
            }
        }

        for (uint count_bf = 0; count_bf < num_bloom_filters; count_bf++) {
            for (uint count_ind = 0;
                 count_ind < bf_partitioned_indices_per_hash.at(count_bf).size();
                 count_ind++) {
                idf_norms.at(count_bf) +=
                    pow(idfs.at(count_hasher).
                        at(bf_partitioned_indices_per_hash.at(count_bf).at(count_ind)), 2);
            }
        }
    }

    for (uint count_bf = 0; count_bf < num_bloom_filters; count_bf++) {
        idf_norms.at(count_bf) = pow(idf_norms.at(count_bf), alpha);
    }
}

void InvertedIndexBloom::perform_query(const vector<uint>& query_hash_indices,
                                       const vector<uint>& query_hash_numbers,
                                       const vector< vector<float> >& idfs,
                                       const vector<float>& idf_norms,
                                       const size_t num_bloom_filters,
                                       vector< pair<float,uint> >& results) {
    // results: first item is score (smaller is better)
    //          second item is group number

    // Initialize results
    results.clear();
    results.resize(num_bloom_filters);
    for (size_t n = 0; n < num_bloom_filters; n++) {
        results.at(n).first = 0;
        results.at(n).second = n;
    }

    // Loop over obtained hashes, increment score for each encountered BF
    uint num_hash_inds = query_hash_indices.size();
    assert(query_hash_numbers.size() == num_hash_inds);
    size_t i = 0; // Hash number (ie, which hash is being used)
    for (size_t count_hash_ind = 0; count_hash_ind < num_hash_inds;
         count_hash_ind++) {
        i = query_hash_numbers.at(count_hash_ind);
        uint hash_ind = query_hash_indices.at(count_hash_ind);
        size_t list_size = partitioned_index_.at(i).at(hash_ind).size();
        float partial_idf_score = pow(idfs.at(i).at(hash_ind), 2);
        for (size_t j = 0; j < list_size; j++) {
            // Negative to follow score ordering convention
            results.at(partitioned_index_.at(i).at(hash_ind).at(j)).first
                -= partial_idf_score
                /idf_norms.at(partitioned_index_.at(i).at(hash_ind).at(j));
        }
    }
}

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "bfindex.h"
#include "bloom.h"
#include "point_indexed/point_index_io.h"

using namespace std;

/************************************************
PUBLIC FUNCTIONS
 ***********************************************/


BFIndex::BFIndex() {
    inverted_index_bloom_ptr_ = NULL;
    num_bloom_filters_ = 0;
    residual_length_ = 0;
    hash_functions_.clear();
    num_hashers_ = 0;
    num_bits_ = 0;
    verbose_ = 0;    
    alpha_ = 1;
    idfs_.clear();
    idf_norms_.clear();
}

BFIndex::~BFIndex() {
    if (inverted_index_bloom_ptr_ != NULL) {
        delete inverted_index_bloom_ptr_;
    }
    inverted_index_bloom_ptr_ = NULL;
    for (size_t i = 0; i < num_hashers_; i++) {
        if (hash_functions_.at(i) != NULL) {
            delete hash_functions_.at(i);
        }
        hash_functions_.at(i) = NULL;
    }
    hash_functions_.clear();
    idfs_.clear();
    idf_norms_.clear();
}

void BFIndex::initialize(size_t n_bits, size_t n_hashers,
                         size_t residual_len, float alpha,
                         int verbose) {
    if (verbose) cout << "Initializing inverted index Bloom filter..." << endl;
    num_hashers_ = n_hashers;
    num_bits_ = n_bits;
    residual_length_ = residual_len;
    alpha_ = alpha;
    verbose_ = verbose;

    // Generate hash functions
    for (size_t i = 0; i < num_hashers_; i++) {
        hash_functions_.push_back(new GBHHasher(num_bits_,
                                                residual_length_));        
    }

    // Initialize InvertedIndexBloom
    inverted_index_bloom_ptr_ =
        new InvertedIndexBloom(num_bits_, num_hashers_);
    if (verbose_) cout << "Done initializing inverted index Bloom filters." << endl;
}

void BFIndex::insert_from_indexes(const vector<string>& index_paths) {
    if (verbose_) cout << "Inserting point-indexed FVs into BFs" << endl;
    size_t num_bloom_filters_ = index_paths.size();

    // Vectors that collect Fisher-embedded features for each clip 
    vector< vector< vector<uint> > > emb_features(num_bloom_filters_);
    vector< vector< vector<uint> > > hash_numbers(num_bloom_filters_);

    for (size_t c = 0; c < num_bloom_filters_; c++) {
        vector < vector < float > > dummy_vec;
        read_bfv_point_index(index_paths.at(c), hash_numbers.at(c),
                             dummy_vec, emb_features.at(c));
    }
    if (verbose_) cout << "All features were read..." << endl;

    size_t number_items = 0;
    for (size_t c = 0; c < num_bloom_filters_; c++) {
        if (verbose_ >= 2) {
            if (c % 100 == 0) {
                cout << "Clip " << c << " out of " << num_bloom_filters_ << endl;
            }
        }
        for (size_t f = 0; f < emb_features.at(c).size(); f++) {
            for (size_t i = 0; i < emb_features.at(c).at(f).size(); i++) {
                uint hash_ind = hash_functions_.at(hash_numbers.at(c).at(f).at(i))
                    ->hash(emb_features.at(c).at(f).at(i));
                inverted_index_bloom_ptr_->insert(hash_ind,
                                                  hash_numbers.at(c).at(f).at(i),
                                                  c);
            }
            number_items++;
        }
        emb_features.at(c).clear();
    }
    if (verbose_) cout << "done inserting!" << endl;

    size_t number_indices;
    get_number_indices(number_indices);
    if (verbose_) cout << "Total number_items = " << number_items
                       << ", number_indices = "
                       << number_indices << endl;

    if (verbose_) cout << "Now, processing idf..." << endl;    
    idf_processing();
    if (verbose_) cout << "done!" << endl;
}

static bool cmp_float_uint_ascend(const pair<float,uint>& pair1,
                                  const pair<float,uint>& pair2) {
    return pair1.first < pair2.first;
}

void BFIndex::perform_query(const vector<uint>& query_emb_features,
                            const vector<uint>& query_hash_numbers,
                            vector< pair<float,uint> >& results) {
    size_t num_features = query_hash_numbers.size();
        
    // First, get hash indices for query_emb_features
    vector<uint> query_hash_indices(num_features, 0);
    for (size_t f = 0; f < num_features; f++) {
        query_hash_indices.at(f) =
            hash_functions_.at(query_hash_numbers.at(f))->
            hash(query_emb_features.at(f));
    }
    
    // Now, query index
    inverted_index_bloom_ptr_->perform_query(query_hash_indices,
                                             query_hash_numbers,
                                             idfs_, idf_norms_,
                                             results);
    // Sorting results
    sort(results.begin(), results.end(), cmp_float_uint_ascend); 
}

/************************************************
PRIVATE FUNCTIONS
 ***********************************************/

void BFIndex::get_number_indices(size_t& number_indices) {
    number_indices = 0;
    inverted_index_bloom_ptr_->get_number_indices(number_indices);
}

void BFIndex::idf_processing() {
    // First, obtain idfs
    idfs_.clear();
    idfs_.resize(num_hashers_);
    size_t num_buckets = pow(2, num_bits_);
    for (uint count_hasher = 0; count_hasher < num_hashers_;
         count_hasher++) {
        idfs_.at(count_hasher).resize(num_buckets, 0);
        vector < uint > counts(num_buckets, 0);
        inverted_index_bloom_ptr_->get_counts_per_hash(count_hasher,
                                                       counts);
        for (uint count_bucket = 0; count_bucket < num_buckets;
             count_bucket++) {
            if (counts.at(count_bucket)) {
                idfs_.at(count_hasher).at(count_bucket)
                    = log(num_bloom_filters_/counts.at(count_bucket));
            }
        }
    }

    // Then, obtain idf_norms
    idf_norms_.clear();
    idf_norms_.resize(num_bloom_filters_);
    inverted_index_bloom_ptr_->get_idf_norm(idfs_, alpha_,
                                            num_bloom_filters_, idf_norms_);
}

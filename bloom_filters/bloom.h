/*
Bloom Filter class
*/

#ifndef BLOOM_H
#define BLOOM_H

using namespace std;

typedef unsigned int uint;

class InvertedIndexBloom {
 public:
    // Constructor/Destructor
    InvertedIndexBloom(const size_t num_bits, const size_t num_hashers);
    ~InvertedIndexBloom();

    // Insert into index
    void insert(const uint hash_ind, const uint hash_number,
                const uint item_number);

    // Function to get number of indices stored in the database
    void get_number_indices(size_t& number_indices);

    // Function to get number of times each bucket is set over the index,
    // for a given hash function
    void get_counts_per_hash(const uint hash_number,
                             vector<uint>& counts);

    void get_idf_norm(const vector< vector<float> >& idfs,
                      const float alpha,
                      const size_t num_bloom_filters,
                      vector<float>& idf_norms);
    
    // Querying the database
    void perform_query(const vector<uint>& query_hash_indices,
                       const vector<uint>& query_hash_numbers,
                       const vector< vector<float> >& idfs,
                       const vector<float>& idf_norms,
                       const size_t num_bloom_filters,
                       vector< pair<float,uint> >& results);
    
 private:
    // Index
    vector< vector< vector<uint> > > partitioned_index_;

    // Parameters
    size_t num_bits_;
    size_t num_hashers_;
};

#endif

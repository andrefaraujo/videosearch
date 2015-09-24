#ifndef GDINDEX_H
#define GDINDEX_H

#include <string>
#include <vector>

#include "../../common/feature_set/feature_set.h"

extern "C" {
#include "../../common/yael_v438_modif/yael/gmm.h"
}

using namespace std;

#ifndef PUSH_BIT
#define PUSH_BIT(packed, bit) \
  packed = packed << 1; \
  packed += bit;
#endif

#ifndef POWER_LAW
#define POWER_LAW(v, a, w)                        \
  w = (v >= 0) ? pow(v,a) : -1 * pow(-v,a);
#endif

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;

// Default initializations for private variables
const uint LD_LENGTH_DEFAULT = 128;
const uint LD_FRAME_LENGTH_DEFAULT = 4;
const string LD_EXTENSION_DEFAULT = ".siftb";
const string LD_NAME_DEFAULT = "sift";
const uint LD_PCA_DIM_DEFAULT = 32;
const float LD_PRE_PCA_POWER_DEFAULT = 0.5;
const uint GD_NUMBER_GAUSSIANS_DEFAULT = 512;
const float GD_POWER_DEFAULT = 0.5;
const uint MIN_NUMBER_WORDS_SELECTED_DEFAULT = 20;
const int WORD_SELECTION_MODE_DEFAULT = 0;
const float WORD_SELECTION_THRESH_DEFAULT = 7;

// Hamming distance above which score is just set to zero
const uint CORR_WEIGHTS_CLIPPING = 16;

class GDIndex
{
 public:

  // Constructor
  GDIndex();

  // Destructor
  ~GDIndex();

  // Index I/O
  // -- write function that is used after indexing, writing 
  //    descriptors, l1 norms and total soft assignment information
  //    per Gaussian
  void write(const string index_path);
  // -- read function that will load index into index_ variables
  //    Note: this will load either L1 norms OR Total Soft Assignment
  //    information, depending on query_parameters_.word_selection_mode
  void read(const string index_path);
  // -- This function writes a file useful in the case of using shots 
  //    with independent keyframes indexing
  void write_frame_list(const string file_path);

  // Clean variables in index_
  void clean_index();

  // Get number of stored global descriptors in index
  uint get_number_global_descriptors();

  // Generate index: these will populate index_
  // -- frame-based signatures
  void generate_index(const vector<string>& feature_files, 
                      const int verbose_level = 1);
  // -- shot/scene-based signatures
  void generate_index_shot_based(const vector<string>& feature_files, 
                                 const vector<uint>& shot_beg_frames,
                                 const int shot_mode, const int shot_keyf, 
                                 const vector < vector < 
                                   pair < uint, uint > > >& track_lists,
                                 const int verbose_level = 1);
  // -- generate one global descriptor from features
  void generate_global_descriptor(const FeatureSet* feature_set, 
                                  vector<uint>& gd_word_descriptor, 
                                  vector<float>& gd_word_l1_norm, 
                                  vector<float>& gd_word_total_soft_assignment);

  // Query index from query's local descriptor path
  void performQuery(const string local_descriptors_path, 
                    vector< pair<float,uint> >& results, 
                    const vector<uint>& indices = vector<uint>(),
                    const uint num_scenes_to_rerank = 0,
                    const uint group_testing_number_centroids = 0,
                    const GDIndex* revv_other_ptr = NULL,
                    const vector < vector < uint > >& vGroupLists 
                      = vector < vector < uint > >(), 
                    const vector < pair < string, pair < uint, uint > > >& shot_info
                      = vector < pair < string, pair < uint, uint > > >(),
                    const int verbose_level = 1);

  // Function to set index_parameters_
  void set_index_parameters(const uint ld_length, const uint ld_frame_length,
                            const string ld_extension, const string ld_name,
                            const uint ld_pca_dim, const float ld_pre_pca_power,
                            const uint gd_number_gaussians, const float gd_power,
                            const string trained_parameters_path,
                            const int verbose_level = 1);

  // Function to set query_parameters_
  // -- Note that the correlation weights loaded here require that
  //    some of the index_parameters_ variables be set. So, ALWAYS
  //    load index_parameters_ BEFORE loading query_parameters_
  void set_query_parameters(const uint min_number_words_selected,
                            const int word_selection_mode,
                            const float word_selection_thresh,
                            const string trained_parameters_path,
                            const int verbose_level = 1);
 private:
  /************ Constants *************/
  // Modes for shot detection
  enum {SHOT_MODE_INDEP_KEYF = 0, SHOT_MODE_SHOT_AGG = 1, SHOT_MODE_GLOBAL_AGG = 2, SHOT_MODE_TRACK_AGG = 3};
  // Mode used in word selection for querying
  enum {WORD_L1_NORM = 0, WORD_SOFT_ASSGN = 1};
  /************ End of Constants *************/

  /************ Variable  *************/  
  // Variables that will hold the index and number of signatures stored
  struct struct_index {
      vector < vector < uint > > word_descriptor;
      vector < vector < float > > word_l1_norms;
      vector < vector < float > > word_total_soft_assignment;

      // Vector that keeps frame numbers that are actually indexed in the db,
      // used only when SHOT_MODE_INDEP_KEYF mode is used
      vector<uint> frame_numbers_in_db;

      // Number of global descriptors in database;
      // this variable is always updated in function update_index
      uint number_global_descriptors;

      // Variables that are used when scoring; these hold values
      // for each database item; these are always updated in
      // function update_index
      vector < uint > number_words_selected;
      vector < float > norm_factors;
  };
  struct_index index_;

  // Index parameters
  struct struct_index_parameters {
      // Local descriptor information
      uint ld_length;
      uint ld_frame_length;
      string ld_extension;
      string ld_name;

      // Parameters for PCA-ing local descriptors
      uint ld_pca_dim;
      float ld_pre_pca_power;
      float* ld_mean_vector;
      vector<float*> ld_pca_eigenvectors;      

      // Parameters used for global descriptor computation
      gmm_t* gd_gmm;
      uint gd_number_gaussians;
      float gd_power;
  };
  struct_index_parameters index_parameters_;

  // Variables relevant for query time
  struct struct_query_parameters {
      // -- Number of minimum words to require for matching
      uint min_number_words_selected;
      // -- Type of word selection mode in use
      // WORD_L1_NORM: only globalWordL1Norm is used
      // WORD_SOFT_ASSGN: only globalWordTotalSoftAssignment is used
      int word_selection_mode;
      // -- Threshold to use in visual word selection (used in asymmetric mode)
      float word_selection_thresh;
      // -- Parameters used in scoring
      float* fast_corr_weights;
      int pop_count[65536];
  };
  struct_query_parameters query_parameters_;

  /************ End of Variables  *************/

  /************ Functions *************/
  // Update index after additions or deletions to it
  void update_index();

  // Sign binarize residuals
  void sign_binarize(const vector<float>& gd_word_residuals, 
                     vector<uint>& gd_word_descriptor);

  // PCA projection for local descriptors
  void project_local_descriptor_pca(const float* desc, float* pca_desc);

  // This function samples the number_frames_out from the shot that begins at
  // frame first_frame and contains number_frames_this_shot, and returns the frame
  // indices in a sorted vector.
  // If one frame is requested (number_frames_out = 1), it will return the center one;
  // otherwise, it will try to take equally spaced frames; if this results in frames 
  // too concentrated at the beginning or the end, it will take only the center ones.
  void sample_frames_from_shot(const uint number_frames_out, 
                               const uint first_frame, 
                               const uint number_frames_this_shot, 
                               vector<uint>& out_frames);

  // Query index with query global descriptor
  void query(const vector<uint>& query_word_descriptor,
             const vector<float>& query_word_l1_norm,
             const vector<float>& query_word_total_soft_assignment,
             const vector<uint>& database_indices,
             vector< pair<float,uint> >& database_scores_indices);

  // Functions to load trained parameters from index_parameters_
  void load_ld_mean_vector(string path);
  void load_ld_pca_eigenvectors(string path);
  void load_gd_gmm(string path);

  // Functions to load trained parameters from query_parameters_
  void load_corr_weights(string path);

  // Helper function to compare pairs
  static bool cmp_float_uint_ascend(const pair<float,uint> pair1, 
                                    const pair<float,uint> pair2) {
      return pair1.first < pair2.first;
  }
  /************ End of Functions *************/


};

#endif

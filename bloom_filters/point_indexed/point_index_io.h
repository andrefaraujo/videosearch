#ifndef POINT_INDEX_IO_H
#define POINT_INDEX_IO_H

#include <vector>

using namespace std;

typedef unsigned int uint;

/* I/O in floating-point format */
void write_fv_point_index(const vector < vector < uint > >& vec_feat_assgns,
                          const vector < vector < float > >& vec_feat_assgn_weights,
                          const vector < vector < vector < float > > >& vec_feat_residuals,
                          const int residual_length,
                          const string file_path);

void read_fv_point_index(const string file_path,
                         vector < vector < uint > >& vec_feat_assgns,
                         vector < vector < float > >& vec_feat_assgn_weights,
                         vector < vector < vector < float > > >& vec_feat_residuals,
                         const int max_feats_to_read = -1);

/* I/O in binary format */
void write_bfv_point_index(const vector < vector < uint > >& vec_feat_assgns,
                           const vector < vector < float > >& vec_feat_assgn_weights,
                           const vector < vector < uint > >& vec_feat_residuals_binarized,
                           const string file_path);

void read_bfv_point_index(const string file_path,
                          vector < vector < uint > >& vec_feat_assgns,
                          vector < vector < float > >& vec_feat_assgn_weights,
                          vector < vector < uint > >& vec_feat_residuals_binarized,
                          const int max_feats_to_read = -1);

#endif

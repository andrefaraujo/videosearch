#ifndef BINARIZE_RESIDUALS_H
#define BINARIZE_RESIDUALS_H

#include <vector>

typedef unsigned int uint;

using namespace std;

// Functions to perform binarization of point-indexed residuals
void binarize_residual(const vector < vector < float > >& feat_residuals,
                       vector<uint>& feat_residuals_binarized);

void binarize_residuals(const vector < vector < vector < float > > >& vec_feat_residuals, 
                        vector < vector < uint > >& vec_feat_residuals_binarized);

#endif


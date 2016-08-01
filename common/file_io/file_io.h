#ifndef FILE_IO_H
#define FILE_IO_H

#include <algorithm>
#include <string>

#include "../feature_set/feature_set.h"

using namespace std;

const int MAX_NUM_FEATURES = 5000;
const float DESC_DIVISOR = 512.0;
const int SIFTGEO_DESC_SIZE = 168;

FeatureSet* readSIFTFile(string sFilename, uint nFrameLength, uint nDescLength);
FeatureSet* readSIFTGeoFile(string sFilename, uint nFrameLength, uint nDescLength);

#endif

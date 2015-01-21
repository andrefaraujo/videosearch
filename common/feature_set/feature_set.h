#ifndef _FEATURE_SET_H_
#define _FEATURE_SET_H_

#include <vector>
#include <cstring>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
using namespace std;

typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short ushort;

/* ----------------------------------------------------------------
 * class: FeatureSet
 * ----------------------------------------------------------------
 * Contains a set of features.
 * ----------------------------------------------------------------
 */
class FeatureSet
{
 public:
  
  /* ----------------------------------------------------------------
   * Constructor: FeatureSet
   * ----------------------------------------------------------------
   * Creates a FeatureSet object.
   * ----------------------------------------------------------------
   */
  FeatureSet(unsigned int nDescriptorLength, unsigned int nFrameLength);
  FeatureSet(FeatureSet* pOther);

  /* ----------------------------------------------------------------
   * Destructor: ~FeatureSet
   * ----------------------------------------------------------------
   * Destroys a FeatureSet object.
   * ----------------------------------------------------------------
   */
  ~FeatureSet();

  /* ----------------------------------------------------------------
   * function: addFeature
   * ----------------------------------------------------------------
   * Adds a new feature to the set.
   * ----------------------------------------------------------------
   */
  void addFeature(float* pDescriptor, float* pFrame);

  /* ----------------------------------------------------------------
   * Method: print
   * ----------------------------------------------------------------
   * Print a summary of the feature set.
   * ----------------------------------------------------------------
   */
  void print();

  // Member variables
  vector<float*> m_vDescriptors;
  vector<float*> m_vFrames;
  unsigned int m_nDescriptorLength;
  unsigned int m_nFrameLength;
  unsigned int m_nNumFeatures;

};

#endif

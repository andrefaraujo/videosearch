#include "feature_set.h"

/* ----------------------------------------------------------------
 * ----------------------------------------------------------------
 * PUBLIC FUNCTIONS
 * ----------------------------------------------------------------
 * ----------------------------------------------------------------
 */

FeatureSet::FeatureSet(unsigned int nDescriptorLength, unsigned int nFrameLength)
{
  m_nDescriptorLength = nDescriptorLength;
  m_nFrameLength = nFrameLength;
  m_nNumFeatures = 0;
}

FeatureSet::FeatureSet(FeatureSet* pOther)
{
  if (pOther == NULL)
    {
      return;
    }
  m_nDescriptorLength = pOther->m_nDescriptorLength;
  m_nFrameLength = pOther->m_nFrameLength;
  m_nNumFeatures = 0;
  for (int nFeature = 0; nFeature < pOther->m_nNumFeatures; nFeature++)
    {
      float* pDescCopy = new float[m_nDescriptorLength];
      float* pFrameCopy = new float[m_nFrameLength];
      memcpy(pDescCopy, pOther->m_vDescriptors[nFeature], m_nDescriptorLength * sizeof(float));
      memcpy(pFrameCopy, pOther->m_vFrames[nFeature], m_nFrameLength * sizeof(float));
      addFeature(pDescCopy, pFrameCopy);
    } // nFeature
}

FeatureSet::~FeatureSet()
{
  vector<float*>::iterator iElem;
  for (iElem = m_vDescriptors.begin(); iElem != m_vDescriptors.end(); iElem++)
    {
      float* pDescriptor = *iElem;
      if (pDescriptor != NULL)
	{
	  delete [] pDescriptor;
	  *iElem = NULL;
	}
    }

  for (iElem = m_vFrames.begin(); iElem != m_vFrames.end(); iElem++)
    {
      float* pFrame = *iElem;
      if (pFrame != NULL)
	{
	  delete [] pFrame;
	  *iElem = NULL;
	}
    }
  m_nNumFeatures = 0;
}

void FeatureSet::addFeature(float* pDescriptor, float* pFrame)
{
  m_vDescriptors.push_back(pDescriptor);
  m_vFrames.push_back(pFrame);
  m_nNumFeatures++;
}

void FeatureSet::print()
{
	for (int nElem = 0; nElem < m_vDescriptors.size(); nElem++) {
		if (m_nFrameLength > 3) {
		  printf("(%d) frames [%.2f %.2f %.2f %.2f], descriptors [%.2f %.2f %.2f %.2f ... %.2f] \n",
			 nElem+1, m_vFrames[nElem][0], m_vFrames[nElem][1], m_vFrames[nElem][2], m_vFrames[nElem][3],
			 m_vDescriptors[nElem][0], m_vDescriptors[nElem][1],
			 m_vDescriptors[nElem][2], m_vDescriptors[nElem][3],
			 m_vDescriptors[nElem][m_nDescriptorLength-1]);
		} else {
		  printf("(%d) frames [%.2f %.2f], descriptors [%.2f %.2f %.2f %.2f ... %.2f] \n",
			 nElem+1, m_vFrames[nElem][0], m_vFrames[nElem][1],
			 m_vDescriptors[nElem][0], m_vDescriptors[nElem][1],
			 m_vDescriptors[nElem][2], m_vDescriptors[nElem][3],
			 m_vDescriptors[nElem][m_nDescriptorLength-1]);
		}
	}
}


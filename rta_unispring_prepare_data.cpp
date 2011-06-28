#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;


/** 
 Pre-uniformize coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void UniSpring::preUniformize(){
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back(mPoints[i*DIM]);
		mPointsY.push_back(mPoints[i*DIM+1]);
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), mPoints[i*DIM]);
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), mPoints[i*DIM+1]);
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		
		// Scale coordinate
		mPoints[i*DIM] = 1 + RECT_SCALE * index_x / (mNpoints - 1) - RECT_SCALE/2;
		mPoints[i*DIM+1] = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		// Assign values so that coordinates won't be found again
		*it_preuni_x = -1; // Find better option... This forces input coordinates to be > 0 (are all descriptor values > 0 ?)
		*it_preuni_y = -1;
		
	}	
	
}

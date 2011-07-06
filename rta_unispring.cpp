#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;

Disk::Disk (float r, float cx, float cy) { 
	
	type = shape_disk;
	ratio = 1;
	scale_factor = r;
	shift_scaled_x = cx - r;
	shift_scaled_y = cy - r;
	
}

double Disk::fd_compute2 (double px, double py) {
	
	double mindist = UniSpring::fd_disk(px, py, 1, 1, 1);
	return mindist;
	
}

Square::Square (float s, float llx, float lly) {
	
	type = shape_square;
	ratio = 1;
	scale_factor = s/2;
	shift_scaled_x = llx;
	shift_scaled_y = lly;
	
}

double Square::fd_compute2 (double px, double py) {
	
	double mindist = UniSpring::fd_rect(px, py, 0, 0, 2, 2);
	return mindist;
	
}


Rectangle::Rectangle (float llx, float lly, float urx, float ury) {
	
	type = shape_rect; 
	ratio = (urx-llx)/(ury-lly);
	scale_factor = (ury-lly)/2;
	shift_scaled_x = llx;
	shift_scaled_y = lly;
	
}

double Rectangle::fd_compute2 (double px, double py) {

	double mindist = UniSpring::fd_rect(px, py, 0, 0, 2*ratio, 2);
	return mindist;

}


UniSpring::UniSpring() {
	
	max_displ_old = 0;
	max_displ_prev = 0;
	dptol = 0.0016;
	stop = 0;
	
};
		

void UniSpring::set_points(int n, float *points, Shape *shape) {
	
	// Get pointer to Shape. Use of pointer allows to use preserve the dynamic type of shape, as defined elsewhere (herited classes Disk, Rectangle, Square)
	mShape = shape;	
	
	// Copy point data
	mNpoints = n;
	mPoints = new coordT[mNpoints*(DIM+1)];	//+1 for convex hull
	mPointsOld = new coordT[mNpoints*(DIM+1)];
	
	for (int i=0; i<mNpoints; i++) {
		
		mPoints[i*DIM] = points[i*DIM];
		mPoints[i*DIM+1] = points[i*DIM+1];
		mPointsOld[i*DIM] = points[i*DIM];
		mPointsOld[i*DIM+1] = points[i*DIM+1];		
		
	}
	
	// Init total force vectors
	std::vector<double> F_temp(2,0); 
	Ftot.resize(mNpoints, F_temp);

	preUniformize();
	setupQhull();
	triangulate();
	getEdgeVector();

}

/** 
 Pre-uniformize x-coordinates between 1 - RECT_SCALE*mShape.ratio/2 and 1 + RECT_SCALE*mShape.ratio/2, y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
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
		
		// Scale coordinates
		mPoints[i*DIM] = mShape->ratio + RECT_SCALE * mShape->ratio * index_x / (mNpoints - 1) - RECT_SCALE/2 * mShape->ratio;
		mPoints[i*DIM+1] = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		// Assign values so that coordinates won't be found again
		*it_preuni_x = -1; // Find better option... This forces input coordinates to be > 0 (are all descriptor values > 0 ?) // TODO: check this problem
		*it_preuni_y = -1;
		
	}	
	
}


void UniSpring::get_points_scaled(float *points) {
	
	//printf("%f",mShape->scale_factor);
	double stop;
		
	for (int i=0; i<mNpoints; i++) {
		
		points[i*DIM] = mPoints[i*DIM] * mShape->scale_factor + mShape->shift_scaled_x;
		points[i*DIM+1] = mPoints[i*DIM+1] * mShape->scale_factor + mShape->shift_scaled_y;	
		
	}
	
}

int UniSpring::update() {
			
	updatePositions();

	if (max_displ_prev / H0 < dptol) stop = 1;
	
	if (stop==0 && max_displ_old / H0 > TTOL) { // Retriangulate
		
		mEdges.clear(); // Reset		
		memcpy(mPointsOld,mPoints,mNpoints*(DIM+1)*sizeof(coordT)); // Copy old points positions
		retriangulate();
		getEdgeVector();
		freeQhullMemory(); // Free memory
		
	}	

	resetPhysicalModel();

	return stop;
	
}

void UniSpring::set_tolerance (float tol) {

	dptol = tol;
	
}

#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

#include <iostream>

using namespace UniSpringSpace ;

Disk::Disk (float r, float cx, float cy) { 
	
	type = shape_disk;
	ratio = 1;
	scale_factor = r;
	shift_scaled_x = cx - r;
	shift_scaled_y = cy - r;
	
}

double Disk::fd_compute (double px, double py) {
	
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

double Square::fd_compute (double px, double py) {
	
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

double Rectangle::fd_compute (double px, double py) {

	double mindist = UniSpring::fd_rect(px, py, 0, 0, 2*ratio, 2);
	return mindist;

}

Sphere::Sphere (float r, float cx, float cy, float cz) {
	
	type = shape_3D_sphere; 
	ratio_1 = 1;
	ratio_2 = 1;
	scale_factor = r;
	shift_scaled_x = cx - r; // TODO: check. should be ok
	shift_scaled_y = cy - r;
	shift_scaled_y = cz - r;
	
}

double Sphere::fd_compute (double px, double py, double pz) {
	
	double mindist = UniSpring::fd_sphere(px, py, pz, 1, 1, 1, 1); // error with 1,1,1,1, 2,1,1,1 seems ok. why ?
	return mindist;
	
}

Cube::Cube (float s, float llbx, float llby, float llbz) {
	
	type = shape_3D_cube; 
	ratio_1 = 1;
	ratio_2 = 1;
	scale_factor = s/2;
	shift_scaled_x = llbx;
	shift_scaled_y = llby;
	shift_scaled_z = llbz;
	
}

double Cube::fd_compute (double px, double py, double pz) {
	
	double mindist = UniSpring::fd_rparallel(px, py, pz, 0, 0, 0, 2, 2, 2);
	return mindist;
	
}

RParallel::RParallel (float llbx, float llby, float llbz, float urtx, float urty, float urtz) {
	
	type = shape_3D_cube; 
	ratio_1 = (urtx-llbx)/(urtz-llbz);
	ratio_2 = (urty-llby)/(urtz-llbz);
	scale_factor = (urtz-llbz)/2;
	shift_scaled_x = llbx;
	shift_scaled_y = llby;
	shift_scaled_z = llbz;
	
}

double RParallel::fd_compute (double px, double py, double pz) {
	
	double mindist = UniSpring::fd_rparallel(px, py, pz, 0, 0, 0, 2*ratio_1, 2*ratio_2, 2);
	return mindist;
	
}

UniSpring::UniSpring() {
	
	max_displ_old = 0;
	max_displ_prev = 0;
	dptol = 0.0016;
	stop = 0;
	
};
		
// TODO : check if mNpoints is 0. Also check later
void UniSpring::set_points(int n, float *points, Shape *shape, bool preUni) {
	
	// Get pointer to Shape. Use of pointer allows to use preserve the dynamic type of shape, as defined elsewhere (herited classes Disk, Rectangle, Square)
	mShape = shape;	
	
	// Copy point data
	mNpoints = n;
	//mPoints = new coordT[mNpoints*(DIM+1)];	//+1 for convex hull
	
	// Allocate space
	mPoints = new hed::Node[mNpoints];
	mPointsOld = new hed::Node[mNpoints];
	

	//nodes.push_back(&mInputPoints[mNInputPoints]);

	
	for (int i=0; i<mNpoints; i++) {
		
		mPoints[i].setPosition(points[i*DIM],points[i*DIM+1]);
		mPointsOld[i].setPosition(points[i*DIM],points[i*DIM+1]);
		nodes.push_back(&mPoints[i]);	
		
	}
	
	// Init total force vectors
	std::vector<double> F_temp(2,0); 
	Ftot.resize(mNpoints, F_temp);

	// Scaling, optional pre-uniformisation
	if (preUni == true) preUniformize();
	else scale();	
	
	//setupQhull();
	//triangulate();
	
	// Triangulate
	triang.createDelaunay(nodes.begin(), nodes.end()); //new_end if using unique
	
	getEdgeVector();

}

void UniSpring::set_points_3D(int n, float *points, Shape_3D *shape) {
	
	// Get pointer to Shape_3D. Use of pointer allows to use preserve the dynamic type of shape, as defined elsewhere (herited classes Cube, Sphere, RParrallel)
	mShape_3D = shape;	
	
	// Copy point data
	mNpoints = n;
	mPoints = (hed::Node*) malloc(sizeof(hed::Node)*mNpoints);
	mPointsOld = (hed::Node*) malloc(sizeof(hed::Node)*mNpoints);
	
	for (int i=0; i<mNpoints; i++) {
		
		mPoints[i].setPosition(points[i*DIM],points[i*DIM+1],points[i*DIM+2]);
		mPointsOld[i].setPosition(points[i*DIM],points[i*DIM+1],points[i*DIM+2]);
		nodes.push_back(&mPoints[i]);
		
	}
	
	// Init total force vectors
	std::vector<double> F_temp(3,0); 
	Ftot.resize(mNpoints, F_temp);
	
	preUniformize_3D();
	//setupQhull();
	//triangulate();
	
	triang.createDelaunay(nodes.begin(), nodes.end()); //new_end if using unique

	getEdgeVector_3D();
	
}

/** 
 Pre-uniformize x-coordinates between 1 - RECT_SCALE*mShape.ratio/2 and 1 + RECT_SCALE*mShape.ratio/2, y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void UniSpring::preUniformize(){
		
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back(mPoints[i].x());
		mPointsY.push_back(mPoints[i].y());
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), mPoints[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), mPoints[i].y());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		
		// Scale coordinates
		
		double scaledX = mShape->ratio + RECT_SCALE * mShape->ratio * index_x / (mNpoints - 1) - RECT_SCALE/2 * mShape->ratio;
		double scaledY = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		mPoints[i].setPosition(scaledX, scaledY);
		
		// Assign values so that coordinates won't be found again
		*it_preuni_x = -1; // Find better option... This forces input coordinates to be > 0 (are all descriptor values > 0 ?) // TODO: check this problem
		*it_preuni_y = -1;
		
	}	
	
}

/** 
 Pre-uniformize x-coordinates between 1 - RECT_SCALE*mShape.ratio/2 and 1 + RECT_SCALE*mShape.ratio/2, y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void UniSpring::preUniformize_3D(){
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back(mPoints[i].x());
		mPointsY.push_back(mPoints[i].y());
		mPointsZ.push_back(mPoints[i].z());
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());
	sort (mPointsZ.begin(), mPointsZ.end());
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	std::vector<double>::iterator it_preuni_z;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), mPoints[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), mPoints[i].y());
		it_preuni_z = find (mPointsZ.begin(), mPointsZ.end(), mPoints[i].z());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		int index_z = it_preuni_z - mPointsZ.begin();
		
		// Scale coordinates
		
		double scaledX = mShape_3D->ratio_1 + RECT_SCALE * mShape_3D->ratio_1 * index_x / (mNpoints - 1) - RECT_SCALE/2 * mShape_3D->ratio_1;
		double scaledY = mShape_3D->ratio_2 + RECT_SCALE * mShape_3D->ratio_2 * index_y / (mNpoints - 1) - RECT_SCALE/2 * mShape_3D->ratio_2;
		double scaledZ = 1 + RECT_SCALE * index_z / (mNpoints - 1) - RECT_SCALE/2;
		
		mPoints[i].setPosition(scaledX, scaledY, scaledZ);
		
		// Assign values so that coordinates won't be found again
		*it_preuni_x = -1; // Find better option... This forces input coordinates to be > 0 (are all descriptor values > 0 ?) // TODO: check this problem
		*it_preuni_y = -1;
		*it_preuni_z = -1;
		
	}	
	
}

/** 
 Scale x-coordinates between mShape->ratio - RECT_SCALE*mShape.ratio/2 and mShape->ratio + RECT_SCALE*mShape.ratio/2, y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2. Alternative to preUniformise().
 */
void UniSpring::scale(){
	
	double minX = mPoints[0].x();
	double maxX = mPoints[0].x();
	double minY = mPoints[0].y();
	double maxY = mPoints[0].y();
		
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if (mPoints[i].x() < minX) minX = mPoints[i].x();	
		if (mPoints[i].x() > maxX) maxX = mPoints[i].x();
		if (mPoints[i].y() < minY) minY = mPoints[i].y();
		if (mPoints[i].y() > maxY) maxY = mPoints[i].y();
		
	}
	
	std::cout << "minX " << minX << " maxX " << maxX << "\n";
	std::cout << "minY " << minY << " maxY " << maxY << "\n";

	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = mShape->ratio + RECT_SCALE * mShape->ratio * (mPoints[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * mShape->ratio;
		double scaledY = 1 + RECT_SCALE * (mPoints[i].y()-minY)/(maxY-minY) - RECT_SCALE/2;
		
		mPoints[i].setPosition(scaledX, scaledY);
		
		
	}	
	
}

/** 
 Scale x-coordinates between mShape->ratio - RECT_SCALE*mShape.ratio/2 and mShape->ratio + RECT_SCALE*mShape.ratio/2, y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2. Alternative to preUniformise().
 */
//void UniSpring::scale_3D(){ // TODO: update for 3D


void UniSpring::get_points_scaled(float *points) {
			
	for (int i=0; i<mNpoints; i++) {
		
		points[i*DIM] = mPoints[i].x() * mShape->scale_factor + mShape->shift_scaled_x;
		points[i*DIM+1] = mPoints[i].y() * mShape->scale_factor + mShape->shift_scaled_y;	
		
	}
	
}

void UniSpring::get_points_scaled_3D(float *points) {
			
	for (int i=0; i<mNpoints; i++) {
		
		points[i*DIM] = mPoints[i].x() * mShape_3D->scale_factor + mShape_3D->shift_scaled_x;
		points[i*DIM+1] = mPoints[i].y() * mShape_3D->scale_factor + mShape_3D->shift_scaled_y;
		points[i*DIM+2] = mPoints[i].z() * mShape_3D->scale_factor + mShape_3D->shift_scaled_z;
		
	}
	
}

int UniSpring::update() {
				
	if (max_displ_old / H0 > TTOL) { // Retriangulate
		
		mEdges.clear(); // Reset		
		memcpy(mPointsOld,mPoints,mNpoints*sizeof(hed::Node)); // Copy old points positions
		triang.createDelaunay(nodes.begin(), nodes.end());
		getEdgeVector();
		//freeQhullMemory(); // Free memory
		
	}	
	
	resetPhysicalModel();
	updatePositions();

	if (max_displ_prev / H0 < dptol) stop = 1;

	return stop;
	
}

int UniSpring::update_3D() {
	
	if (max_displ_old / H0 > TTOL) { // Retriangulate
		
		mEdges.clear(); // Reset		
		memcpy(mPointsOld,mPoints,mNpoints*(DIM+1)*sizeof(hed::Node)); // Copy old points positions
		triang.createDelaunay(nodes.begin(), nodes.end());
		getEdgeVector_3D();
		
	}	
	
	resetPhysicalModel_3D();
	updatePositions_3D();
	
	if (max_displ_prev / H0 < dptol) stop = 1;
	
	return stop;
	
}

void UniSpring::set_tolerance (float tol) {

	dptol = tol;
	
}

std::vector< std::vector<int> > UniSpring::get_edges() {
	
	return mEdges;
	
}

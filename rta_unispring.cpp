#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;

UniSpring::UniSpring() {

	max_displ_old = 0;
	max_displ_prev = 0;
	dptol = 0.0015;
	stop = 0;
	
	
}

void UniSpring::set_points(int n, double *points, Shape shape) { // TODO: points as float

	// Copy point data
	mNpoints = n;
	mPoints = new coordT[mNpoints*(DIM+1)];	
	mPointsOld = new coordT[mNpoints*(DIM+1)];	
	memcpy(mPoints,points,mNpoints*(DIM+1)*sizeof(coordT));
	memcpy(mPointsOld,points,mNpoints*(DIM+1)*sizeof(coordT));
	
	// Init total force vectors
	std::vector<double> F_temp(2,0); 
	Ftot.resize(mNpoints, F_temp);

	preUniformize();
	setupQhull();
	triangulate();
	getEdgeVector();

}

void UniSpring::get_points_scaled(double *points) {
	
	for (int i=0; i<mNpoints; i++) {
		
		points[i*DIM]=mPoints[i*DIM];
		points[i*DIM+1]=mPoints[i*DIM+1];
		
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

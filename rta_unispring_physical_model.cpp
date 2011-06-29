#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;


/** 
 Compute euclidean distance between points i1 and i2.
 */
double UniSpring::euclDistance(int i1, int i2){
	
	double diffx = mPoints[i1*DIM] - mPoints[i2*DIM];
	double diffy = mPoints[i1*DIM+1] - mPoints[i2*DIM+1];
	double dist = pow(diffx,2) + pow(diffy,2);
	
    return sqrt(dist);
	
}


/** 
 Compute displacement of point i1 from its previous position.
 */
double UniSpring::euclDispl(int i1){
	
	double diffx = mPoints[i1*DIM] - mPointsOld[i1*DIM];
	double diffy = mPoints[i1*DIM+1] - mPointsOld[i1*DIM+1];
	double dist = pow(diffx,2) + pow(diffy,2);
	
    return sqrt(dist);
	
}


/** 
 Sum all elements of a vector.
 */
double UniSpring::sum(std::vector<double> v){
	
	double sum = 0;
	
	for(std::vector<double>::iterator it = v.begin(); it!=v.end(); it++){
		
		sum += *it;
	}
	
	return sum;
	
}


/** 
 Definition of the signed distance function (square)
 */
double UniSpring::fd_rect(double px, double py, double llx, double lly, double urx, double ury){
	
	double mindist = -std::min(std::min(std::min(-lly+py,ury-py),-llx+px),urx-px);
	return mindist;
	
}
						  
/** 
 Definition of the signed distance function (disk)
 */
double UniSpring::fd_disk(double px, double py, double r, double cx, double cy){
							  
	double mindist = sqrt(pow(px-cx,2)+pow(py-cy,2)) - r;
	return mindist;
																										
}

double UniSpring::fd_compute(double px, double py){
	
	if (mShape.type == shape_square) {
		
		double mindist = fd_rect(px, py, 0, 2, 0, 2);
		return mindist;
		
	}
	
	if (mShape.type == shape_disk) {
		
		double mindist = fd_disk(px, py, 1, 1, 1);
		return mindist;
		
	}
	
	if (mShape.type == shape_rect) {
		
		double mindist = fd_rect(px, py, 0, 0, 2*mShape.ratio, 2);
		return mindist;
		
	}

}



/** 
 Definition of the desired length function.
 */
double UniSpring::fh(double px, double py){
	
	double desdist = 1; // Uniform
	return desdist;
	
}

/** 
 Update point positions.
 */
void UniSpring::updatePositions(){
	
	double deps = sqrt(EPS)*H0;
	
	// Compute edge lengths and initial target distances
	for (int i = 0; i < mEdges.size(); i++) {
		
		double middlex;
		double middley;
		double length;
		
		middlex = (mPoints[mEdges[i][0]*DIM]+mPoints[mEdges[i][1]*DIM])/2;
		middley = (mPoints[mEdges[i][0]*DIM+1]+mPoints[mEdges[i][1]*DIM+1])/2;
		length = euclDistance(mEdges[i][0],mEdges[i][1]);
		
		double length2 = pow(length,2);
		
		hbars2.push_back(pow(fh(middlex,middley),2)); // compute squared value of target distance on edge's middle
		L.push_back(length);
		L2.push_back(pow(length,2));
		
	}
	
	// Compute total edge length and initial target distance	
	hbars2_sum = sum(hbars2);
	L2_sum = sum(L2);
	
	// Compute total force components
	for (int i = 0; i < mEdges.size(); i++) {
		
		std::vector<double> Fvec(2);
		std::vector<double> barvec(2);
		
		// Compute edge force
		double L0 = sqrt(hbars2[i]) * FSCALE * sqrt(L2_sum/hbars2_sum);
		double F = std::max(L0-L[i],0.);
		
		// Get edge unit vector
		barvec[0] = ( mPoints[mEdges[i][0]*DIM] - mPoints[mEdges[i][1]*DIM] ) / L[i];
		barvec[1] = ( mPoints[mEdges[i][0]*DIM+1] - mPoints[mEdges[i][1]*DIM+1] ) / L[i];		
		
		// Project force on edge unit vector
		Fvec[0] = F * barvec[0];
		Fvec[1] = F * barvec[1];
		
		// Assign force due to current edge to point indices in Ftot		
		int currentIndex = mEdges[i][0];
		
		Ftot[mEdges[i][0]][0] += Fvec[0];
		Ftot[mEdges[i][0]][1] += Fvec[1];
		Ftot[mEdges[i][1]][0] += -Fvec[0];
		Ftot[mEdges[i][1]][1] += -Fvec[1];
		
	}
	
	// Move points & bring outside points back to boundary
	std::vector<double> displ_temp(2);
	
	for (int i=0; i<mNpoints; i++) {
		
		displ_temp[0] = DELTAT * Ftot[i][0];
		displ_temp[1] = DELTAT * Ftot[i][1];
		
		mPoints[i*DIM] = mPoints[i*DIM] + displ_temp[0];
		mPoints[i*DIM+1] = mPoints[i*DIM+1] + displ_temp[1];
				
		// Check if point has moved outside
		double d = fd_compute(mPoints[i*DIM], mPoints[i*DIM+1]);
		
		if (d > 0) {
			
			//Bring it back to boundary
			double dgradx = ( fd_compute(mPoints[i*DIM] + deps, mPoints[i*DIM+1]) - d ) / deps;
			double dgrady = ( fd_compute(mPoints[i*DIM], mPoints[i*DIM+1] + deps) - d ) / deps;
			mPoints[i*DIM] = mPoints[i*DIM] - d * dgradx; 
			mPoints[i*DIM+1] = mPoints[i*DIM+1] - d * dgrady;			
			
		}
				
		// Compute displacement relative to previous triangulation (_old) and previous step (_prev)
		double displ_old = euclDispl(i);
		double displ_prev = sqrt(pow(displ_temp[0],2) + pow(displ_temp[1],2));
		
		// For interior point, update max displacement relative to previous step
		if (d < -GEPS && displ_prev > max_displ_prev) {
			
			max_displ_prev = displ_prev;
			
		}
		
		// For all points, update max displacement relative to previous triangulation
		if (displ_old > max_displ_old) {
			
			max_displ_old = displ_old;
			
		}
		
	}
	
}


/** 
 Prepare physical model for next iteration.
 */
void UniSpring::resetPhysicalModel(){
	
	max_displ_old = 0;
	max_displ_prev = 0;
	std::vector<double> F_temp(2,0); 
	for (int i=0; i<mNpoints; i++) {
		
		Ftot[i]=F_temp;
		
	}
	hbars2.clear();
	L.clear();
	L2.clear();
	
}
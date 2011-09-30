#include "rta_unispring.h"
#define DIM RTA_UNISPRING_NDIM

using namespace UniSpringSpace ;


/** 
 Compute euclidean distance between points i1 and i2.
 */
double UniSpring::euclDistance(int i1, int i2){
	
	double diffx = mPoints[i1].x() - mPoints[i2].x();
	double diffy = mPoints[i1].y() - mPoints[i2].y();
	double dist = pow(diffx,2) + pow(diffy,2);
	
    return sqrt(dist);
	
}

double UniSpring::euclDistance_3D(int i1, int i2){
	
	double diffx = mPoints[i1].x() - mPoints[i2].x();
	double diffy = mPoints[i1].y() - mPoints[i2].y();
	double diffz = mPoints[i1].z() - mPoints[i2].z();
	double dist = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
	
    return sqrt(dist);
	
}


/** 
 Compute displacement of point i1 from its previous position.
 */
double UniSpring::euclDispl(int i1){
	
	double diffx = mPoints[i1].x() - mPointsOld[i1].x();
	double diffy = mPoints[i1].y() - mPointsOld[i1].y();
	double dist = pow(diffx,2) + pow(diffy,2);
	
    return sqrt(dist);
	
}

double UniSpring::euclDispl_3D(int i1){
	
	double diffx = mPoints[i1].x() - mPointsOld[i1].x();
	double diffy = mPoints[i1].y() - mPointsOld[i1].y();
	double diffz = mPoints[i1].z() - mPointsOld[i1].z();
	double dist = pow(diffx,2) + pow(diffy,2) + pow(diffz,2);
	
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
 Definition of the signed distance function (rectangle)
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

/** 
 Definition of the signed distance function (right parallelepiped)
 */
double UniSpring::fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz){
	
	double mindist = -std::min(std::min(std::min(std::min(std::min(-llbz+pz,urtz-pz),-llby+py),urty-py),-llbx+px),urtx-px);
	return mindist;
	//    The formula used here is not quite correct.  In particular, it is wrong
	//    for points exterior to the cube whose nearest point on the cube is at a corner.
	//	
	//   For DISTMESH_3D's purposes, though, this computation is accurate enough.
	
	//d = - min ( min ( min ( min ( min ( -0.0+p(:,3), 1.0-p(:,3) ), ...
	//							 -0.0+p(:,2) ), ...
	//					   1.0-p(:,2) ), ...
	//				 -0.0+p(:,1) ), ...
	//		   1.0-p(:,1) );
	
}

/** 
 Definition of the signed distance function (sphere)
 */
double UniSpring::fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz){
	
	double mindist = sqrt(pow(px-cx,2)+pow(py-cy,2)+pow(pz-cz,2)) - r;
	return mindist;
	
}

//UNUSED
/*
double UniSpring::fd_compute(double px, double py){
	
	if (mShape->type == shape_square) {
		
		double mindist = fd_rect(px, py, 0, 2, 0, 2);
		return mindist;
		
	}
	
	if (mShape->type == shape_disk) {
		
		double mindist = fd_disk(px, py, 1, 1, 1);
		return mindist;
		
	}
	
	if (mShape->type == shape_rect) {
		
		double mindist = fd_rect(px, py, 0, 0, 2*mShape->ratio, 2);
		return mindist;
		
	}

}
 */



/** 
 Definition of the desired length function.
 */
double UniSpring::fh(double px, double py){
	
	double desdist = 1; // Uniform
	return desdist;
	
}

double UniSpring::fh_3D(double px, double py, double pz){
	
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
		
		middlex = (mPoints[mEdges[i][0]].x()+mPoints[mEdges[i][1]].x())/2;
		middley = (mPoints[mEdges[i][0]].y()+mPoints[mEdges[i][1]].y())/2;
		length = euclDistance(mEdges[i][0],mEdges[i][1]);
				
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
		barvec[0] = ( mPoints[mEdges[i][0]].x() - mPoints[mEdges[i][1]].x() ) / L[i];
		barvec[1] = ( mPoints[mEdges[i][0]].y() - mPoints[mEdges[i][1]].y() ) / L[i];		
		
		// Project force on edge unit vector
		Fvec[0] = F * barvec[0];
		Fvec[1] = F * barvec[1];
		
		// Assign force due to current edge to point indices in Ftot		
		//int currentIndex = mEdges[i][0]; //debug
		
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
		
		mPoints[i].setPosition( mPoints[i].x() + displ_temp[0], mPoints[i].y() + displ_temp[1]);	
		
		// Check if point has moved outside
		double d = mShape->fd_compute(mPoints[i].x(), mPoints[i].y());
		
		if (d > 0) {
			
			//Bring it back to boundary
			double dgradx = ( mShape->fd_compute(mPoints[i].x() + deps, mPoints[i].y()) - d ) / deps;
			double dgrady = ( mShape->fd_compute(mPoints[i].x(), mPoints[i].y() + deps) - d ) / deps;
			
			mPoints[i].setPosition( mPoints[i].x() - d * dgradx, mPoints[i].y() - d * dgrady);
			
			
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

void UniSpring::updatePositions_3D(){
	
	double deps = sqrt(EPS)*H0;
	
	// Compute edge lengths and initial target distances
	for (int i = 0; i < mEdges.size(); i++) {
		
		double middlex;
		double middley;
		double middlez;
		double length;
		
		middlex = (mPoints[mEdges[i][0]].x()+mPoints[mEdges[i][1]].x())/2;
		middley = (mPoints[mEdges[i][0]].y()+mPoints[mEdges[i][1]].y())/2;
		middlez = (mPoints[mEdges[i][0]].z()+mPoints[mEdges[i][1]].z())/2;
		length = euclDistance_3D(mEdges[i][0],mEdges[i][1]);
		
		//double length3 = pow(length,3);
		
		hbars3.push_back(pow(fh_3D(middlex,middley,middlez),3)); // compute squared value of target distance on edge's middle
		L.push_back(length);
		L3.push_back(pow(length,3));
		
	}
	
	// Compute total edge length and initial target distance	
	hbars3_sum = sum(hbars3);
	L3_sum = sum(L3);
	
	// Compute total force components
	for (int i = 0; i < mEdges.size(); i++) {
		
		std::vector<double> Fvec(3);
		std::vector<double> barvec(3);
		
		// Compute edge force
		double L0 = pow(hbars3[i],(double)1/3) * FSCALE * pow(L3_sum/hbars3_sum,(double)1/3);
		double F = std::max(L0-L[i],0.);
		
		// Get edge unit vector
		barvec[0] = ( mPoints[mEdges[i][0]].x() - mPoints[mEdges[i][1]].x() ) / L[i];
		barvec[1] = ( mPoints[mEdges[i][0]].y() - mPoints[mEdges[i][1]].y() ) / L[i];	
		barvec[2] = ( mPoints[mEdges[i][0]].z() - mPoints[mEdges[i][1]].z() ) / L[i];
		
		// Project force on edge unit vector
		Fvec[0] = F * barvec[0];
		Fvec[1] = F * barvec[1];
		Fvec[2] = F * barvec[2];
		
		// Assign force due to current edge to point indices in Ftot		
		
		Ftot[mEdges[i][0]][0] += Fvec[0];
		Ftot[mEdges[i][0]][1] += Fvec[1];
		Ftot[mEdges[i][0]][2] += Fvec[2];
		Ftot[mEdges[i][1]][0] += -Fvec[0];
		Ftot[mEdges[i][1]][1] += -Fvec[1];
		Ftot[mEdges[i][1]][2] += -Fvec[2];
		
	}
	
	// Move points & bring outside points back to boundary
	std::vector<double> displ_temp(3);
	
	for (int i=0; i<mNpoints; i++) {
		
		displ_temp[0] = DELTAT * Ftot[i][0];
		displ_temp[1] = DELTAT * Ftot[i][1];
		displ_temp[2] = DELTAT * Ftot[i][2];
		
		mPoints[i].setPosition(mPoints[i].x() + displ_temp[0], mPoints[i].y() + displ_temp[1], mPoints[i].z() + displ_temp[2]);
		
		// Check if point has moved outside
		double d = mShape_3D->fd_compute(mPoints[i].x(), mPoints[i].y(), mPoints[i].z());
		
		if (d > 0) {
			
			//Bring it back to boundary
			double dgradx = ( mShape_3D->fd_compute(mPoints[i].x() + deps, mPoints[i].y(), mPoints[i].z()) - d ) / deps;
			double dgrady = ( mShape_3D->fd_compute(mPoints[i].x(), mPoints[i].y() + deps, mPoints[i].z()) - d ) / deps;
			double dgradz = ( mShape_3D->fd_compute(mPoints[i].x(), mPoints[i].y(), mPoints[i].z() + deps) - d ) / deps;
			
			mPoints[i].setPosition(mPoints[i].x() - d * dgradx, mPoints[i].y() - d * dgrady, mPoints[i].z() - d * dgradz);
					
		}
		
		// Compute displacement relative to previous triangulation (_old) and previous step (_prev)
		double displ_old = euclDispl_3D(i);
		double displ_prev = sqrt(pow(displ_temp[0],2) + pow(displ_temp[1],2) + pow(displ_temp[2],2));
		
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

void UniSpring::resetPhysicalModel_3D(){
	
	max_displ_old = 0;
	max_displ_prev = 0;
	std::vector<double> F_temp(3,0); 
	for (int i=0; i<mNpoints; i++) {
		
		Ftot[i]=F_temp;
		
	}
	hbars3.clear();
	L.clear();
	L3.clear();
	
}
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
	
	double mindist = fd_disk(px, py, 1, 1, 1);
	return mindist;
	
}

/** 
 Definition of the signed distance function (disk)
 */
double Disk::fd_disk(double px, double py, double r, double cx, double cy){
	
	double mindist = sqrt(pow(px-cx,2)+pow(py-cy,2)) - r;
	return mindist;
	
}

/** 
 Pre-uniformize x-coordinates between mShape.ratio - RECT_SCALE*mShape.ratio/2 and mShape.ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void Disk::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints){
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		
		// Scale coordinates
		
		double scaledX = ratio + RECT_SCALE * ratio * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		
	}	
	
}

/** 
 Scale x-coordinates between ratio - RECT_SCALE*mShape.ratio/2 and ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2. Alternative to preUniformise().
 */
void Disk::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio + RECT_SCALE * ratio * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
	}	
	
}

Square::Square (float s, float llx, float lly) {
	
	type = shape_square;
	ratio = 1;
	scale_factor = s/2;
	shift_scaled_x = llx;
	shift_scaled_y = lly;
	
}

double Square::fd_compute (double px, double py) {
	
	double mindist = Rectangle::fd_rect(px, py, 0, 0, 2, 2);
	return mindist;
	
}

/** 
 Pre-uniformize x-coordinates between mShape.ratio - RECT_SCALE*mShape.ratio/2 and mShape.ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void Square::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints){
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		
		// Scale coordinates
		
		double scaledX = ratio + RECT_SCALE * ratio * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		
	}	
	
}

/** 
 Scale x-coordinates between ratio - RECT_SCALE*mShape.ratio/2 and ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2. Alternative to preUniformise().
 */
void Square::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio + RECT_SCALE * ratio * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
	}	
	
}

Rectangle::Rectangle (float llx, float lly, float urx, float ury) {
	
	type = shape_rect; 
	ratio = (urx-llx)/(ury-lly);
	scale_factor = (ury-lly)/2;
	shift_scaled_x = llx;
	shift_scaled_y = lly;
	
}

double Rectangle::fd_compute (double px, double py) {
	
	double mindist = fd_rect(px, py, 0, 0, 2*ratio, 2);
	return mindist;
	
}

/** 
 Definition of the signed distance function (rectangle)
 */
double Rectangle::fd_rect(double px, double py, double llx, double lly, double urx, double ury){
	
	double mindist = -std::min(std::min(std::min(-lly+py,ury-py),-llx+px),urx-px);
	return mindist;
	
}

/** 
 Pre-uniformize x-coordinates between mShape.ratio - RECT_SCALE*mShape.ratio/2 and mShape.ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 */
void Rectangle::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints){
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		
		// Scale coordinates
		
		double scaledX = ratio + RECT_SCALE * ratio * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * index_y / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		
	}	
	
}


/** 
 Scale x-coordinates between ratio - RECT_SCALE*mShape.ratio/2 and ratio + RECT_SCALE*mShape.ratio/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2. Alternative to preUniformise().
 */
void Rectangle::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio + RECT_SCALE * ratio * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio;
		double scaledY = 1 + RECT_SCALE * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
	}	
	
}


static bool PointEqual (std::vector<double>& a, std::vector<double>& b)
{
#   define EPSILON 1e-6
    return abs(a[0] - b[0]) < EPSILON  &&  abs(a[1] - b[1]) < EPSILON;
}


// check validity of point list for polygon constructor: number and unicity of points
bool UniSpringSpace::Polygon::CheckPolygon(std::vector< std::vector<double> > vertices)
{
    int num = vertices.size();
    // check number of points for at least a triangle
    if (num >= 3)
    {
	// check unicity of vertices
	for (int i = 0; i < num; i++)
	    for (int j = i + 1; j < num; j++)
		if (PointEqual(vertices[i], vertices[j]))
		    return false;

	return true;
    }
    else
	return false;
}

UniSpringSpace::Polygon::Polygon () {
	
	vx.clear();
	vy.clear();
	
	std::vector<double> P1(2);
	std::vector<double> P2(2);
	std::vector<double> P3(2);
	std::vector<double> P4(2);
	P1[0] = 0; P1[1] = 0;
	P2[0] = 0; P2[1] = 100;
	P3[0] = 100; P3[1] = 100;
	P4[0] = 100; P4[1] = 0;
	
	std::vector< std::vector<double> > vertices;
	vertices.push_back(P1);
	vertices.push_back(P2);
	vertices.push_back(P3);
	vertices.push_back(P4);
	
	type = shape_poly; 
	ratio = 1;
	scale_factor = 50;
	shift_scaled_x = 0;
	shift_scaled_y = 0;
	minVX = 0;
	maxVX = 100;
	minVY = 0;
	maxVY = 100;		
	
	// Scale polygon to fit a bounding box such as llx = 0, lly = 0, urx = 2*ratio, ury = 2)
	vx.push_back(0); vy.push_back(0);
	vx.push_back(0); vy.push_back(2);
	vx.push_back(2); vy.push_back(2);
	vx.push_back(2); vy.push_back(0);
	
	nVertices = 4;
	
	UniSpringSpace::Polygon::close();
	
}

UniSpringSpace::Polygon::Polygon (std::vector< std::vector<double> > vertices) {
	
	vx.clear();
	vy.clear();
	
        //todo: check if (vertices.size() > 0) --> exception
        
	minVX = vertices[0][0];
	maxVX = vertices[0][0];
	minVY = vertices[0][1];
	maxVY = vertices[0][1];
	
	nVertices = vertices.size();
	
	// Find min & max
	for (int i=0; i<nVertices; i++) {
		if (vertices[i][0] < minVX) minVX = vertices[i][0];	
		if (vertices[i][0] > maxVX) maxVX = vertices[i][0];
		if (vertices[i][1] < minVY) minVY = vertices[i][1];
		if (vertices[i][1] > maxVY) maxVY = vertices[i][1];
	}
	
	type = shape_poly; 
	ratio = (maxVX-minVX)/(maxVY-minVY);
	scale_factor = (maxVY-minVY)/2;
	shift_scaled_x = minVX;
	shift_scaled_y = minVY;
	
	// Scale polygon to fit a bounding box such as llx = 0, lly = 0, urx = 2*ratio, ury = 2
	for (int i=0; i<nVertices; i++) {
		vx.push_back( (vertices[i][0] - minVX) / (maxVX - minVX) * 2 * ratio );
		vy.push_back( (vertices[i][1] - minVY) / (maxVY - minVY) * 2 );
	}
	
	UniSpringSpace::Polygon::close();
	
}

/** 
 Definition of the signed distance function for a polygon of vertex coordinates vx and vy
 */
double UniSpringSpace::Polygon::fd_compute (double px, double py) {
	
	double mindist = fd_poly(px, py, vx, vy, nVertices);
	return mindist;
	
}


double UniSpringSpace::Polygon::fd_poly(double px, double py, std::vector<double> &vx, std::vector<double> &vy, int nVertices) {
	
	// Linear parameters of segments that connect the vertices (Ax + By + C = 0)
	std::vector<double> A;
	std::vector<double> B;
	std::vector<double> C;
	
	for (int i=1; i<nVertices; i++) {
		A.push_back( -(vy[i] - vy[i-1]) );
		B.push_back( (vx[i] - vx[i-1]) );
		C.push_back( ( vy[i]*vx[i-1] - vx[i]*vy[i-1] ) );
	}
	
	// Find the projection of point (x,y) on each rib. A, B are of size nVertices-1
	std::vector<double> projection_x;
	std::vector<double> projection_y;
	
	for (int i=0; i<nVertices-1; i++) {
		double AB = 1l / ( pow(A[i],2) + pow(B[i],2) );
		double vv =  A[i]*px + B[i]*py + C[i];		
		projection_x.push_back( px - A[i]*AB*vv );
		projection_y.push_back( py - B[i]*AB*vv );
	}
	
	// Test for the case where a polygon rib is either horizontal or vertical.
	for (int i=0; i<nVertices-1; i++) {
		if (B[i] == 0l) 
			projection_x[i] = vx[i];
		if (A[i] == 0l) 
			projection_y[i] = vy[i];		
	}
	
	// Find all cases where projected point is between the two rib's vertices (inside the rib's segment)
	std::vector<int> idx;
	std::vector<double> dp;
	
	for (int i=0; i<nVertices-1; i++) {
		bool x_condition = ( ( ( projection_x[i] >= vx[i]) && (projection_x[i] <= vx[i+1]) ) || ( ( projection_x[i] >= vx[i+1] ) && ( projection_x[i] <= vx[i] ) ) );
		bool y_condition = ( ( ( projection_y[i] >= vy[i]) && (projection_y[i] <= vy[i+1]) ) || ( ( projection_y[i] >= vy[i+1] ) && ( projection_y[i] <= vy[i] ) ) );
		if (x_condition && y_condition) {
			idx.push_back(i);
			dp.push_back(sqrt( pow(projection_x[i] - px, 2) + pow(projection_y[i] - py, 2) ));
		}
	}
	
	// Distance from point (x,y) to the vertices
	std::vector<double> dv;
	
	for (int i=0; i<nVertices-1; i++)
		dv.push_back( sqrt(pow(vx[i]-px,2) + pow(vy[i]-py,2)) );
	
	// Compute min of dp and dv
	double minDv = dv[0];
	
	for (int i=0; i<nVertices-1; i++) {
		if (dv[i]<minDv)
			minDv = dv[i];
	}
	
	double minDp;
	if (dp.size() > 0) {
		minDp = dp[0];
		
		for (int i=0; i<dp.size(); i++) {
			if (dp[i]<minDp)
				minDp = dp[i];
		}
	}
	
	// Compute distance from (px,py) to polygon
	double distance;
	
	if (idx.size() == 0) // If all projections are outside rib segments
		distance = minDv; // The distance is then the distance to the closest vertex
	else 
		distance = std::min(minDv,minDp);
	
	if (UniSpringSpace::Polygon::isInPoly(px, py, vx, vy))
		distance = - distance;
	
	return distance;
	
}

/** 
 Pre-uniformize x-coordinates within the square inscribed in the polygon's incenter.
 */
void UniSpringSpace::Polygon::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints){
	
	double inscribedCircleCenterX;
	double inscribedCircleCenterY;
	double inscribedCircleRadius;
	
	findInscribedCircle(&inscribedCircleCenterX, &inscribedCircleCenterY, &inscribedCircleRadius);
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// sort coordinates
	sort (mPointsX.begin(), mPointsX.end());
	sort (mPointsY.begin(), mPointsY.end());			
	
	// Pre-uniformize and scale coordinates within the square inscribed in the polygon's inscribed circle
	double squareSide = sqrt(2) * inscribedCircleRadius;
	std::vector<double>::iterator it_preuni_x;
	std::vector<double>::iterator it_preuni_y;
	
	for (int i=0; i<mNpoints; i++) {
		
		int sizeX_debug = mPointsX.size();
		int sizeY_debug = mPointsY.size();
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		int index_x = it_preuni_x - mPointsX.begin(); // To get the true index, add i which is the number of previously found and erased elements
		int index_y = it_preuni_y - mPointsY.begin();
		
		double scaledX = inscribedCircleCenterX + squareSide * index_x / (mNpoints - 1) - squareSide/2;
		double scaledY = inscribedCircleCenterY + squareSide * index_y / (mNpoints - 1) - squareSide/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
						
	}	
	
}


/** 
 Scale x-coordinates within the square inscribed in the polygon's incenter. Alternative to preUniformise().
 */
void UniSpringSpace::Polygon::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double inscribedCircleCenterX;
	double inscribedCircleCenterY;
	double inscribedCircleRadius;
	
	findInscribedCircle(&inscribedCircleCenterX, &inscribedCircleCenterY, &inscribedCircleRadius);
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		
	}
	
	// Scale coordinates within the square inscribed in the polygon's inscribed circle
	double squareSide = sqrt(2) * inscribedCircleRadius;
	
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = inscribedCircleCenterX + squareSide * ((*mPoints)[i].x()-minX)/(maxX-minX) - squareSide/2;
		double scaledY = inscribedCircleCenterY + squareSide * ((*mPoints)[i].y()-minY)/(maxY-minY) - squareSide/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY);
		
	}
	
}

/** 
 Approximate search of the interior point with the largest distance to the polygon (defined as the minimal distance to its ribs)
 */
void UniSpringSpace::Polygon::findInscribedCircle(double *inscribedCircleCenterX, double *inscribedCircleCenterY, double *inscribedCircleRadius) {
	
	double x = 0; // The polygon has been scaled in the constructor within a bounding box such as llx = 0, lly = 0, urx = 2*ratio, ury = 2
	double y = 0;
	*inscribedCircleRadius = 0;
	
	//Create grid points inside the polygon, with a resolution of POLYGON_GRID_RES on each side of the bounding box
	while (x < 2*ratio) {
		
		while (y < 2) {
			
			if ( fd_compute(x, y) < 0 && -fd_compute(x, y) > *inscribedCircleRadius ) {
				*inscribedCircleRadius = -fd_compute(x, y);
				*inscribedCircleCenterX = x;
				*inscribedCircleCenterY = y;
			}
			
			y += POLYGON_GRID_RES * 2; // 2 = (2 - 0) = ury - lly of bounding box
			
		}
		
		x += POLYGON_GRID_RES * (2*ratio); // 2*ratio = (2*ratio - 0) = urx - llx of bounding box
		y = 0;
		
	}
	
}


void UniSpringSpace::Polygon::close() {
	
	if (vx[nVertices-1] != vx[0] || vx[nVertices-1] != vx[0]) {
		
		vx.push_back(vx[0]);
		vy.push_back(vy[0]);
		
	}
	
	nVertices++;
	
}

/*
 * Winding number test for a point in a polygon
 * Input:   P = a point,
 *          V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
 * Return:  wn = the winding number (=0 only if P is outside V[])
 */ 
bool UniSpringSpace::Polygon::isInPoly(double px, double py, std::vector<double> vx, std::vector<double> vy) {
	
	int wn = 0;    // the winding number counter
	
	// loop through all edges of the polygon
	for (int i=0; i<vx.size()-1; i++) {   // edge from V[i] to V[i+1] // Last vertex is "virtual" (closes the polygon) : ignore it
		if (vy[i] <= py) {         // start y <= P.y
			if (vy[i+1] > py)      // an upward crossing
				if (isLeft( vx[i], vy[i], vx[i+1], vy[i+1], px, py) > 0)  // P left of edge
					++wn;            // have a valid up intersect
		}
		else {                       // start y > P.y (no test needed)
			if (vy[i+1] <= py)     // a downward crossing
				if (isLeft( vx[i], vy[i], vx[i+1], vy[i+1], px, py) < 0)  // P right of edge
					--wn;            // have a valid down intersect
		}
	}
	return wn != 0;
	
}


/* 
 * Tests if a point is Left|On|Right of an infinite line.
 * Input:  three points P0, P1, and P2
 * Return: >0 for P2 left of the line through P0 and P1
 *         =0 for P2 on the line
 *         <0 for P2 right of the line
 */
double UniSpringSpace::Polygon::isLeft( double P0x, double P0y, double P1x, double P1y, double P2x, double P2y ) {
	
	return ( (P1x - P0x) * (P2y - P0y) - (P2x - P0x) * (P1y - P0y) );
	
}



Sphere::Sphere (float r, float cx, float cy, float cz) {
	
	type = shape_3D_sphere; 
	ratio_1 = 1;
	ratio_2 = 1;
	scale_factor = r;
	shift_scaled_x = cx - r;
	shift_scaled_y = cy - r;
	shift_scaled_y = cz - r;
	
}

double Sphere::fd_compute (double px, double py, double pz) {
	
	double mindist = Sphere::fd_sphere(px, py, pz, 1, 1, 1, 1); 
	return mindist;
	
}

/** 
 Definition of the signed distance function (sphere)
 */
double Sphere::fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz){
	
	double mindist = sqrt(pow(px-cx,2)+pow(py-cy,2)+pow(pz-cz,2)) - r;
	return mindist;
	
}


/** 
 Pre-uniformize x-coordinates between mShape_3D->ratio1 - RECT_SCALE*mShape_3D.ratio1/2 and 1 + RECT_SCALE*mShape_3D.ratio1/2, 
 y-coordinates between mShape_3D->ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 z-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2,
 */
void Sphere::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints) {
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	std::vector<double> mPointsZ;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	double maxZ = (*mPoints)[0].z();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		
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
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		it_preuni_z = find (mPointsZ.begin(), mPointsZ.end(), (*mPoints)[i].z());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		int index_z = it_preuni_z - mPointsZ.begin();
		
		// Scale coordinates
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * index_y / (mNpoints - 1) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * index_z / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		*it_preuni_z = maxZ + 1;
		
	}	
	
}

/** 
 Scale x-coordinates between ratio1 - RECT_SCALE*ratio1/2 and ratio + RECT_SCALE*mShape_3D.ratio1/2,
 y-coordinates between ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 Alternative to preUniformise().
 */
void Sphere::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	double minZ = (*mPoints)[0].z();
	double maxZ = (*mPoints)[0].z();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() < minZ) minZ = (*mPoints)[i].z();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * ((*mPoints)[i].z()-minZ)/(maxZ-minZ) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
	}	
	
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
	
	double mindist = RParallel::fd_rparallel(px, py, pz, 0, 0, 0, 2, 2, 2);
	return mindist;
	
}

/** 
 Pre-uniformize x-coordinates between ratio1 - RECT_SCALE*mShape_3D.ratio1/2 and 1 + RECT_SCALE*mShape_3D.ratio1/2, 
 y-coordinates between ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 z-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2,
 */
void Cube::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints) {
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	std::vector<double> mPointsZ;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	double maxZ = (*mPoints)[0].z();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		
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
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		it_preuni_z = find (mPointsZ.begin(), mPointsZ.end(), (*mPoints)[i].z());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		int index_z = it_preuni_z - mPointsZ.begin();
		
		// Scale coordinates
		
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * index_y / (mNpoints - 1) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * index_z / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		*it_preuni_z = maxZ + 1;
		
	}	
	
}

/** 
 Scale x-coordinates between ratio1 - RECT_SCALE*ratio1/2 and ratio + RECT_SCALE*mShape_3D.ratio1/2,
 y-coordinates between ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 Alternative to preUniformise().
 */
void Cube::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	double minZ = (*mPoints)[0].z();
	double maxZ = (*mPoints)[0].z();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() < minZ) minZ = (*mPoints)[i].z();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * ((*mPoints)[i].z()-minZ)/(maxZ-minZ) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
	}	
	
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
	
	double mindist = RParallel::fd_rparallel(px, py, pz, 0, 0, 0, 2*ratio_1, 2*ratio_2, 2);
	return mindist;
	
}

/** 
 Definition of the signed distance function (right parallelepiped)
 */
double RParallel::fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz){
	
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
 Pre-uniformize x-coordinates between ratio1 - RECT_SCALE*mShape_3D.ratio1/2 and 1 + RECT_SCALE*mShape_3D.ratio1/2, 
 y-coordinates between ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 z-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2,
 */
void RParallel::preUniformize(std::vector<hed::Node> *mPoints, int mNpoints){
	
	std::vector<double> mPointsX;
	std::vector<double> mPointsY;
	std::vector<double> mPointsZ;
	double maxX = (*mPoints)[0].x();
	double maxY = (*mPoints)[0].y();
	double maxZ = (*mPoints)[0].z();
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	for (int i=0; i<mNpoints; i++) {
		
		mPointsX.push_back((*mPoints)[i].x());
		mPointsY.push_back((*mPoints)[i].y());
		mPointsZ.push_back((*mPoints)[i].z());
		
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
		
		it_preuni_x = find (mPointsX.begin(), mPointsX.end(), (*mPoints)[i].x());
		it_preuni_y = find (mPointsY.begin(), mPointsY.end(), (*mPoints)[i].y());
		it_preuni_z = find (mPointsZ.begin(), mPointsZ.end(), (*mPoints)[i].z());
		
		int index_x = it_preuni_x - mPointsX.begin();
		int index_y = it_preuni_y - mPointsY.begin();
		int index_z = it_preuni_z - mPointsZ.begin();
		
		// Scale coordinates
		
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * index_x / (mNpoints - 1) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * index_y / (mNpoints - 1) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * index_z / (mNpoints - 1) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
		// Assign values such as these coordinates won't be found again
		*it_preuni_x = maxX + 1;
		*it_preuni_y = maxY + 1;
		*it_preuni_z = maxZ + 1;
		
	}	
	
}

/** 
 Scale x-coordinates between ratio1 - RECT_SCALE*ratio1/2 and ratio + RECT_SCALE*mShape_3D.ratio1/2,
 y-coordinates between ratio2 - RECT_SCALE*mShape_3D.ratio2/2 and 1 + RECT_SCALE*mShape_3D.ratio2/2,
 y-coordinates between 1 - RECT_SCALE/2 and 1 + RECT_SCALE/2.
 Alternative to preUniformise().
 */
void RParallel::scale (std::vector<hed::Node> *mPoints, int mNpoints) {
	
	double minX = (*mPoints)[0].x();
	double maxX = (*mPoints)[0].x();
	double minY = (*mPoints)[0].y();
	double maxY = (*mPoints)[0].y();
	double minZ = (*mPoints)[0].z();
	double maxZ = (*mPoints)[0].z();
	
	// Find min & max
	for (int i=0; i<mNpoints; i++) {
		
		if ((*mPoints)[i].x() < minX) minX = (*mPoints)[i].x();	
		if ((*mPoints)[i].x() > maxX) maxX = (*mPoints)[i].x();
		if ((*mPoints)[i].y() < minY) minY = (*mPoints)[i].y();
		if ((*mPoints)[i].y() > maxY) maxY = (*mPoints)[i].y();
		if ((*mPoints)[i].z() < minZ) minZ = (*mPoints)[i].z();
		if ((*mPoints)[i].z() > maxZ) maxZ = (*mPoints)[i].z();
		
	}
	
	// Scale coordinates
	for (int i=0; i<mNpoints; i++) {
		
		double scaledX = ratio_1 + RECT_SCALE * ratio_1 * ((*mPoints)[i].x()-minX)/(maxX-minX) - RECT_SCALE/2 * ratio_1;
		double scaledY = ratio_2 + RECT_SCALE * ratio_2 * ((*mPoints)[i].y()-minY)/(maxY-minY) - RECT_SCALE/2 * ratio_2;
		double scaledZ = 1 + RECT_SCALE * ((*mPoints)[i].z()-minZ)/(maxZ-minZ) - RECT_SCALE/2;
		
		(*mPoints)[i].setPosition(scaledX, scaledY, scaledZ);
		
	}	
	
}


UniSpring::UniSpring() {
	
	max_displ_old = 0;
	max_displ_prev = 0;
	dptol = 0.0016;
	stop = 0;
	
};

// TODO : check if mNpoints is 0. Also check later
void UniSpring::set_points (int n, int cols, float *points, Shape *shape, bool preUni) {
	
	// Get pointer to Shape. Use of pointer allows to use preserve the dynamic type of shape, as defined elsewhere (herited classes Disk, Rectangle, Square)
	mShape = shape;	
	
	triang.cleanAll();
	
	// Copy point data
	mNpoints = n;
	//mPoints = new coordT[mNpoints*(DIM+1)];	//+1 for convex hull
	
	// Allocate space
	mPoints.resize(mNpoints);
	mPointsOld.resize(mNpoints);
	nodes.resize(mNpoints);
	
	for (int i = 0; i < mNpoints; i++) {
		mPoints[i].init(i, points[i * cols], points[i * cols + 1]);
		mPointsOld[i].init(i, points[i * cols], points[i * cols + 1]);
		nodes[i] = &mPoints[i];	
	}
	
	// Init total force vectors
	std::vector<double> F_temp(2,0); 
	Ftot.resize(mNpoints, F_temp);
	
	// Scaling, optional pre-uniformisation
	if (preUni == true) mShape->preUniformize(&mPoints, mNpoints);
	else mShape->scale(&mPoints, mNpoints);	
	
	// Triangulate
	triang.createDelaunay(nodes.begin(), nodes.end()); //new_end if using unique
	getEdgeVector();
	
}

void UniSpring::set_points_3D(int n, float *points, Shape_3D *shape, bool preUni) {
		
	// Get pointer to Shape_3D. Use of pointer allows to use preserve the dynamic type of shape, as defined elsewhere (herited classes Cube, Sphere, RParrallel)
	mShape_3D = shape;	
	
	triang.cleanAll();
	
	// Copy point data
	mNpoints = n;
	mPoints.resize(mNpoints);
	mPointsOld.resize(mNpoints);
	
	for (int i=0; i<mNpoints; i++) {
		mPoints[i].init(i, points[i*DIM],points[i*DIM+1],points[i*DIM+2]);
		mPointsOld[i].init(i, points[i*DIM],points[i*DIM+1],points[i*DIM+2]);
		nodes[i] = &mPoints[i];	
	}
	
	// Init total force vectors
	std::vector<double> F_temp(3,0); 
	Ftot.resize(mNpoints, F_temp);
	
	// Scaling, optional pre-uniformisation
	if (preUni == true) mShape_3D->preUniformize(&mPoints, mNpoints);
	else mShape_3D->scale(&mPoints, mNpoints);	

	// Triangulate
	triang.createDelaunay(nodes.begin(), nodes.end()); //new_end if using unique
	getEdgeVector_3D();
	
}

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
		mPointsOld = mPoints; // Copy old points positions
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
		mPointsOld = mPoints; // Copy old points positions
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

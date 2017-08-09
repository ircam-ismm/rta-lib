/**
 * @file rta_unispring.h
 * @author Riccardo Borghesi
 *
 * @copyright
 * Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * All rights reserved.
 * 
 * License (BSD 3-clause)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _RTA_UNISPRING
#define _RTA_UNISPRING

#include <limits>
#include <numeric>
#include <vector>
#include "halfedge/HeTriang.h"
#include "halfedge/HeDart.h"
#include "halfedge/HeTraits.h"

#define RTA_UNISPRING_NDIM 2

// Physical model parameters / TODO change name to match RTA names. TODO : include elsewhere ? (auto selection 2D/3D)
#define H0 0.5 // Mean initial distance between points (for scaling)
#define TTOL 0.1
#define EPS 2.2204e-16
#define GEPS 0.001*H0
#define MAX_EDGE_LENGTH 0.2 // Max length edge allowed in triangulated structure
#define FSCALE 1.2 // Must be >1 to help points spread accross the whole target region. 1.2 is ok for 2D
//#define FSCALE 1.1 // given by distmesh_3D empirical formula for 3D
#define DELTAT 0.2 // ok for 2D
//#define DELTAT 0.1 // 0.1 better in 3D

// Scale factor
#define POLYGON_GRID_RES 0.01
#define RECT_SCALE sqrt(2) //ok for 2D
//#define RECT_SCALE sqrt(2)/2 // Must be < 1. sqrt(2) is ok for 2, too much for 3D (points are pre-uniformized outside the target region). sqrt(2)/2 seems ok

namespace UniSpringSpace 
{
	
	
typedef enum { shape_disk, shape_square, shape_rect, shape_poly } shape_enum_t; // 2D
typedef enum { shape_3D_sphere, shape_3D_cube, shape_3D_rparallel } shape_3D_enum_t; // 3D

class Shape {
public:
	virtual double fd_compute(double px, double py) { return 0.0; };
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints) { };
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints) { };
	shape_enum_t type;
	float ratio; // width/height (of bounding box)
	float scale_factor;
	float shift_scaled_x;
	float shift_scaled_y;
//private:
//	virtual bool isInPoly(double px, double py) { };
};

class Disk : public Shape 
{
public:
    Disk (float r = 0.5, float cx = 0.5, float cy = 0.5);
	virtual double fd_compute (double px, double py);
	static double fd_disk(double px, double py, double r, double cx, double cy);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};

class Square : public Shape 
{
public:
    Square (float s = 1, float llx = 0, float lly = 0);
	virtual double fd_compute (double px, double py);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};

class Rectangle : public Shape
{
public:
    Rectangle (float llx = 0, float lly = 0, float urx = 1, float ury = 1);
	virtual double fd_compute (double px, double py);
	static double fd_rect(double px, double py, double llx, double lly, double urx, double ury);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};
	
class Polygon : public Shape
{
public:
    // check validity of point list for polygon constructor: number and unicity of points
    static bool CheckPolygon(std::vector< std::vector<double> > vertices);

	Polygon ();
	Polygon (std::vector< std::vector<double> > vertices);
	virtual double fd_compute (double px, double py);
	static double fd_poly(double px, double py, std::vector<double> &vx, std::vector<double> &vy, int nVertices);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
	static bool isInPoly(double px, double py, std::vector<double> vx, std::vector<double> vy);
	static double isLeft( double P0x, double P0y, double P1x, double P1y, double P2x, double P2y );
private:
	void findInscribedCircle (double *inscribedCircleCenterX, double *inscribedCircleCenterY, double *inscribedCircleRadius);
	void close(); // Check if polygon is closed, close it if necessary
	std::vector<double> vx; // vertex x-coordinates
	std::vector<double> vy; // vertex x-coordinates
	int nVertices;
	double minVX; // Min/max coordinates of vertices (i.e. bounding box, i.e. interior points), before scaling
	double maxVX;
	double minVY;
	double maxVY;
		
};
	
	
class Shape_3D {
public:
	virtual double fd_compute(double px, double py, double pz) { return 0.0; };
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints) { };
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints) { };
	shape_3D_enum_t type;
	float ratio_1; // width/depth (of bounding box)
	float ratio_2; // height/depth (of bounding box)
	float scale_factor;
	float shift_scaled_x;
	float shift_scaled_y;
	float shift_scaled_z;
		
};
		
class Sphere : public Shape_3D 
{
public:
	Sphere (float r = 1, float cx = 1, float cy = 1, float cz = 1);
	virtual double fd_compute (double px, double py, double pz);
	static double fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};
	
class Cube : public Shape_3D
{
public:
	Cube (float s = 1, float llbx = 0, float llby = 0, float llbz = 0);
	virtual double fd_compute (double px, double py, double pz);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};	
	
class RParallel : public Shape_3D 
{
public:
	
	RParallel (float llbx = 0, float llby = 0, float llbz = 0, float urtx = 1, float urty = 1, float urtz = 1);
	virtual double fd_compute (double px, double py, double pz);
	static double fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz);
	virtual void preUniformize(std::vector<hed::Node> *mPoints, int mNpoints);
    virtual void scale (std::vector<hed::Node> *mPoints, int mNpoints);
};


class UniSpring
{
public:
    UniSpring();	// constructor
    ~UniSpring()	// destructor
    {
	triang.cleanAll();
	//TODO: clean up all members?
    }
					  
    /** set points and initialise unispring algorithm:
	copy points array, pre-uniformise (if preUni = true), do first triangulation
	@param n	number of points
	@param points	pointer to points x, y data
	@param shape	defines shape
     */
	void set_points (int n, int cols, float *points, Shape *shape, bool preUni = true);
	
	/** set points and initialise unispring algorithm:
	 copy points array, pre-uniformise, do first triangulation
	 @param n	number of points
	 @param points	pointer to points x, y, z data
	 @param shape_3D	defines 3D shape
     */
	void set_points_3D (int n, float *points, Shape_3D *shape, bool preUni = true);
	
	//void set_points (int n, double *points, Shape shape);

    /** copy points to given pointer, scaled to dimension given by shape definition
     */
    void get_points_scaled (float *points);
	void get_points_scaled_3D (float *points);
	//void get_points_scaled (double *points);
	
	/** get triangulation edges as vector of pairs of points ids.
	 */
	std::vector< std::vector<int> > get_edges();
	int get_num_edges() { return mEdges.size(); };

    /** run one update step
	@return stop flag, true if movement under tolerance
     */
    int  update ();
	int  update_3D ();

    void set_tolerance (float tol);
	
	//static double fd_disk(double px, double py, double r, double cx, double cy); // TODO: redefine as Shape methods
	//static double fd_rect(double px, double py, double llx, double lly, double urx, double ury);
	//static double fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz);
	static double fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz);
	//static double fd_poly(double px, double py, std::vector<double> vx, std::vector<double> vx);

private:
	double euclDistance(int i1, int i2);
	double euclDistance_3D(int i1, int i2);
	double euclDispl(int i1);
	double euclDispl_3D(int i1);
	double fh(double px, double py);
	double fh_3D(double px, double py, double pz);
	double sum(std::vector<double> v);
	void removeDuplicateEdges();
	void loadData();
	void updatePositions();
	void updatePositions_3D();
	void getEdgeVector();
	void getEdgeVector_3D();
	void resetPhysicalModel();
	void resetPhysicalModel_3D();
	void print_summary();
	
	// TTL
	hed::Triangulation      triang;
	std::vector<hed::Node*> nodes; // vector of pointers to point coordinates data
	std::vector<hed::Node>  mPoints;
	std::vector<hed::Node>  mPointsOld;
	
	// Physical model
	Shape *mShape;
	Shape_3D *mShape_3D;
	std::vector<double> mPointsX; // For preuniformisation step
	std::vector<double> mPointsY;
	std::vector<double> mPointsZ;
	std::vector< std::vector<int> > mEdges; // Edge point indexes
	int mNpoints; // Number of points
	std::vector<double> hbars2; // Desired lengths (squared) (2D)
	double hbars2_sum;
	std::vector<double> hbars3; // Desired lengths (cubed) (3D)
	double hbars3_sum;
	std::vector<double> L2; // Current lengths (squared) (2D)
	double L2_sum;
	std::vector<double> L3; // Current lengths (cubed) (3D)
	double L3_sum;
	std::vector<double> L; // Current lengths
	std::vector< std::vector<double> > Ftot; // Current total force components
	std::vector<double> displacements; // Distance traveled during previous iteration
	double max_displ_old, max_displ_prev; // Maximum displacement (interior points). Old / prev: since last triangulation / previous iteration
	double dptol; // Stop criterion
	int stop;    

};

} // end namespace UniSpring

bool operator<(std::vector<int> const& v1, std::vector<int> const& v2);
bool operator==(std::vector<int> const& v1, std::vector<int> const& v2);
#endif

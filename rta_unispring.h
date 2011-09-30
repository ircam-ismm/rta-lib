#ifndef _RTA_UNISPRING
#define _RTA_UNISPRING

extern "C" {
#define qh_QHimport
#include "qhull_a.h"
}

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

// Scale factor for display
#define RECT_SCALE sqrt(2) //ok for 2D
//#define RECT_SCALE sqrt(2)/2 // Must be < 1. sqrt(2) is ok for 2, too much for 3D (points are pre-uniformized outside the target region). sqrt(2)/2 seems ok

namespace UniSpringSpace 
{
	
	
typedef enum { shape_disk, shape_square, shape_rect } shape_enum_t; // 2D
typedef enum { shape_3D_sphere, shape_3D_cube, shape_3D_rparallel } shape_3D_enum_t; // 3D

class Shape {
public:
	virtual double fd_compute (double px, double py) {printf("base");}; // TODO remove printf "base" (for debug: base should'nt be print)
    shape_enum_t type;
	float ratio; // width/height (of bounding box)
	float scale_factor;
	float shift_scaled_x;
	float shift_scaled_y;

};

class Disk : public Shape 
{
public:
		
    Disk (float r = 1, float cx = 1, float cy = 1);
	virtual double fd_compute (double px, double py);
};

class Square : public Shape 
{
public:
    Square (float s = 1, float llx = 0, float lly = 0);
	virtual double fd_compute (double px, double py);
};

class Rectangle : public Shape
{
public:
    Rectangle (float llx = 0, float lly = 0, float urx = 1, float ury = 1);
	virtual double fd_compute (double px, double py);
};
	
class Shape_3D {
public:
	virtual double fd_compute (double px, double py, double pz) {printf("base");}; // TODO remove printf "base" (for debug)
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
};
	
class Cube : public Shape_3D
{
public:
	Cube (float s = 1, float llbx = 0, float llby = 0, float llbz = 0);
	virtual double fd_compute (double px, double py, double pz);
};	
	
class RParallel : public Shape_3D 
{
public:
	
	RParallel (float llbx = 0, float llby = 0, float llbz = 0, float urtx = 1, float urty = 1, float urtz = 1);
	virtual double fd_compute (double px, double py, double pz);
};


class UniSpring
{
public:
	UniSpring();
					  
    /** set points and initialise unispring algorithm:
	copy points array, pre-uniformise (if preUni = true), do first triangulation
	@param n	number of points
	@param points	pointer to points x, y data
	@param shape	defines shape
     */
    void set_points (int n, float *points, Shape *shape, bool preUni = true);
	
	/** set points and initialise unispring algorithm:
	 copy points array, pre-uniformise, do first triangulation
	 @param n	number of points
	 @param points	pointer to points x, y, z data
	 @param shape_3D	defines 3D shape
     */
	void set_points_3D (int n, float *points, Shape_3D *shape);
	
	//void set_points (int n, double *points, Shape shape);

    /** copy points to given pointer, scaled to dimension given by shape definition
     */
    void get_points_scaled (float *points);
	void get_points_scaled_3D (float *points);
	//void get_points_scaled (double *points);
	
	/** get triangulation edges as vector of pairs of points ids.
	 */
	std::vector< std::vector<int> > get_edges();

    /** run one update step
	@return stop flag, true if movement under tolerance
     */
    int  update ();
	int  update_3D ();

    void set_tolerance (float tol);
	
	static double fd_disk(double px, double py, double r, double cx, double cy); // TODO: redefine as Shape methods
	static double fd_rect(double px, double py, double llx, double urx, double lly, double ury);
	static double fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz);
	static double fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz);


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
	void preUniformize();
	void preUniformize_3D();
	void scale();
	void updatePositions();
	void updatePositions_3D();
	void getEdgeVector();
	void getEdgeVector_3D();
	void resetPhysicalModel();
	void resetPhysicalModel_3D();
	void print_summary();

	
	// TTL
	hed::Triangulation triang;
	std::vector<hed::Node*> nodes; // vector of pointers to point coordinates data
	hed::Node *mPoints;
	hed::Node *mPointsOld;
	
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

#ifndef _RTA_UNISPRING
#define _RTA_UNISPRING

extern "C" {
#define qh_QHimport
#include "qhull_a.h"
}

#include <limits>
#include <numeric>
#include <vector>

#define RTA_UNISPRING_NDIM 3

// Physical model parameters / TODO change name to match RTA names
#define H0 0.5 // Mean initial distance between points (for scaling)
#define TTOL 0.1
#define EPS 2.2204e-16
#define GEPS 0.001*H0
#define FSCALE 1.2 // Must be >1 to help points spread accross the whole target region
#define DELTAT 0.2

// Scale factor for display
#define RECT_SCALE sqrt(2) // Must be < 1

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
	copy points array, pre-uniformise, do first triangulation
	@param n	number of points
	@param points	pointer to points x, y data
	@param shape	defines shape
     */
    void set_points (int n, float *points, Shape *shape);
	
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
	//void get_points_scaled (double *points);

    /** run one update step
	@return stop flag, true if movement under tolerance
     */
    int  update ();

    void set_tolerance (float tol);
	
	static double fd_disk(double px, double py, double r, double cx, double cy);
	static double fd_rect(double px, double py, double llx, double urx, double lly, double ury);
	static double fd_sphere(double px, double py, double pz, double r, double cx, double cy, double cz);
	static double fd_rparallel(double px, double py, double pz, double llbx, double llby, double llbz, double urtx, double urty, double urtz);


private:
	double euclDistance(int i1, int i2);
	double euclDispl(int i1);
	double fh(double px, double py);
	double sum(std::vector<double> v);
	void loadData();
	void preUniformize();
	void preUniformize_3D();
	void setupQhull();
	void triangulate();
	void updatePositions();
	void getEdgeVector();
	void retriangulate();
	void freeQhullMemory();
	void resetPhysicalModel();
	void print_summary();
	int qh_new_qhull2(int dim, int numpoints, coordT *points, boolT ismalloc,
					  char *qhull_cmd, FILE *outfile, FILE *errfile);

	// Qhull
	coordT *mPoints; // Current point coordinates
	coordT *mPointsOld; // Point coordinates used in last triangulation
	coordT **rows;
	boolT ismalloc;
	char *flags;
	FILE *outfile;
	FILE *errfile;
	int exitcode;
	facetT *facet;
	qhT *oldqh;
	int curlong, totlong;
	facetT *neighbor;
	vertexT *vertex, **vertexp;
	ridgeT *ridge, **ridgep;
	mergeT *merge, **mergep;
	
	// Physical model
	Shape *mShape;
	Shape_3D *mShape_3D;
	std::vector<double> mPointsX; // For preuniformisation step
	std::vector<double> mPointsY;
	std::vector<double> mPointsZ;
	std::vector< std::vector<int> > mEdges; // Edge point indexes
	int mNpoints; // Number of points
	std::vector<double> hbars2; // Desired lengths (squared)
	double hbars2_sum;
	std::vector<double> L2; // Current lengths (squared)
	double L2_sum;
	std::vector<double> L; // Current lengths (squared)
	std::vector< std::vector<double> > Ftot; // Current total force components
	std::vector<double> displacements; // Distance traveled during previous iteration
	double max_displ_old, max_displ_prev; // Maximum displacement (interior points). Old / prev: since last triangulation / previous iteration
	double dptol; // Stop criterion
	int stop;    

};

} // end namespace UniSpring
#endif

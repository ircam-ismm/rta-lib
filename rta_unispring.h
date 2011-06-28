
namespace UniSpring 
{

#define RTA_UNISPRING_NDIM 2

typedef enum { shape_disk, shape_square, shape_rect } shape_enum_t;

class Shape {
public:
    shape_enum_t type;
};

class Disk : public Shape 
{
public:
    Disk (float r = 1) { type = shape_disk;  radius_ = r; };
    float radius_;
};

class Square : public Shape 
{
public:
    Square (float s = 1) { type = shape_square;  size_ = s; };
    float size_;
};

class Rect : public Shape
{
public:
    Rect (float llx = 0, float lly = 0, float urx = 1, float ury = 1) 
    {
	type = shape_rect; 
	llx_ = llx;
	lly_ = lly;
	urx_ = urx;
	ury_ = ury;
    };
    float llx_, lly_, urx_, ury_;
};


class UniSpring
{
public:
    /** set points and initialise unispring algorithm:
	copy points array, pre-uniformise, do first triangulation
	@param n	number of points
	@param points	pointer to points x, y data
	@param shape	defines shape
     */
    void set_points (int n, float *points, Shape shape);

    /** copy points to given pointer, scaled to dimension given by shape definition
     */
    void get_points_scaled (float *points);

    /** run one update step
	@return stop flag, true if movement under tolerance
     */
    int  update ();

    void set_tolerance (float tol);

private:
    
};

} // end namespace UniSpring



#define RTA_UNISPRING_NDIM 2

typedef enum { shape_disk, shape_square, shape_rect } shape_enum_t;

class Shape {
    shape_enum_t type;
};

class Disk : Shape {
    Disk () { type = shape_disk };
    float radius;
};

class Square : Shape {
    Square () { type = shape_square };
    float size;
};

class Rect : Shape{
    Rect () { type = shape_rect };
    float llx, lly, urx, ury;
};


class UniSpring
{
public:
    void set_points (int n, float *points, Shape *shape);
    int  unpdate ();

private:
    
}

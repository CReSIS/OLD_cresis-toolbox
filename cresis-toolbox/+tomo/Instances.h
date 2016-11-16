// Instances.h: Some very helpful variables and functions.
// By Mingze Xu, July 2016
//
#ifndef __INSTANCES_H__
#define __INSTANCES_H__

#include <sstream>
#include <tuple>

typedef vector< pair<size_t, size_t> > CoordType;
typedef vector< tuple<size_t, size_t, size_t> > PointType;
typedef vector< vector<size_t> > LayerType;

// Middle coordinate
#define MID 33

// Directions
#define dir_up 0
#define dir_down 1
#define dir_left 2
#define dir_right 3
#define dir_all 4

// TRWS gamma
#define gamma 0.5

// TRWS large cost
#define LARGE 1000000000

// Compute square value
template <class T>
T sqr(T x) { return x*x; }

// Convert integer to string
string itos(int i) {
    stringstream s;
    s << i;
    return s.str();
}

// Dynamic smoothness
double norm_pdf(double x, double scale=5.0, double mu=MID, double s=12.0) {
    return scale * (1.0/(s*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/s));
}

// Distance transform
// THE CODE BELOW THIS POINT WAS TAKEN FROM DAVID CRANDALL
void dt(const double *src, double *dst, double *dst_ind, int s1, int s2, int d1, int d2, double scale, int off=0) {
    int d = (d1+d2) >> 1;
    int s = s1;
    for (int p = s1; p <= s2; p++)
        if (src[s] + sqr(s-d-off) * scale > src[p] + sqr(p-d-off) * scale)
            s = p;
    dst[d] = src[s] + sqr(s-d-off) * scale;
    dst_ind[d] = s;

    if(d-1 >= d1)
        dt(src, dst, dst_ind, s1, s, d1, d-1, scale, off);
    if(d2 >= d+1)
        dt(src, dst, dst_ind, s, s2, d+1, d2, scale, off);
}

void dt_1d(const double *f, double scale, double *result, double *dst_ind, int beg, int end, int off=0) {
    dt(f, result, dst_ind, beg, end-1, beg, end-1, scale, off);
}
// END CODE FROM DAVID CRANDALL

#endif

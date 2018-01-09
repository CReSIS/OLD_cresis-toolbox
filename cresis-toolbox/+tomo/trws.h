// instances.h: Some very helpful variables and functions.
// By Mingze Xu, November 2016
//
#ifndef __INSTANCES_H__
#define __INSTANCES_H__

#define _USE_MATH_DEFINES
#include <cmath>
#include <stddef.h>

#include <sstream>
#include <tuple>

// Directions
#define dir_up 0
#define dir_down 1
#define dir_left 2
#define dir_right 3
#define dir_all 4
// Default mu (mean) of smooth
#define MID 32
// Default sigma (variance) of smooth (input argument: smooth_var)
#define SIGMA 24
// Default scale of smooth (input argument: smooth_weight)
#define SCALE 5
// TRWS: large cost
#define LARGE 1000000000
// TRWS: gamma
#define gamma 0.5
// TRWS: max_loops
#define MAX_LOOP 50

typedef vector< pair<size_t, size_t> > CoordType;
typedef vector< tuple<size_t, size_t, size_t> > PointType;
typedef vector< vector<size_t> > LayerType;

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
double norm_pdf(double x, double mu=MID, double sigma=SIGMA, double scale=SCALE) {
    return scale * (1.0/(sigma*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/sigma));
}

// THE CODE BELOW THIS POINT WAS TAKEN FROM DAVID CRANDALL
// Distance transform
// -- Every index from d1 to d2 will be set in dst and dst_ind
// -- dst will contain the minimum value for that destination
// -- dst_ind will contain the minimum source index for that destination
void dt(const double *src, double *dst, double *dst_ind, int s1, int s2, int d1, int d2, double scale, int off=0) {
    int d = (d1+d2) >> 1; // Find the midpoint of the destination
    int s = s1;
    for (int p = s1; p <= s2; p++) // Search through all the sources and find the minimum
        if (src[s] + sqr(s-d-off) * scale > src[p] + sqr(p-d-off) * scale)
            s = p;
    dst[d] = src[s] + sqr(s-d-off) * scale; // Minimum value to the midpoint
    dst_ind[d] = s; // Minimum source index for the midpoint

    if(d-1 >= d1) {
      // Tree search: bottom half of destinations, also constrain sources to less than the minimum src index
        dt(src, dst, dst_ind, s1, s, d1, d-1, scale, off);
    }
    if(d2 >= d+1) {
      // Tree search: top half of destinations, also constrain sources to more than the minimum src index
        dt(src, dst, dst_ind, s, s2, d+1, d2, scale, off);
    }
}

void dt_1d(const double *f, double scale, double *result, double *dst_ind, int beg, int end, int off=0) {
    dt(f, result, dst_ind, beg, end-1, beg, end-1, scale, off);
}
// END CODE FROM DAVID CRANDALL

#endif

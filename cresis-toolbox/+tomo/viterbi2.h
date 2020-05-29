// viterbi2.h
//
// Layer-tracking program based on the Viterbi algorithm
//  for MUSIC-processed 2D and 3D data
//
// Adapted from original code by Mingze Xu, David Crandall, and John Paden
//
// Author: Victor Berger
// See also: viterbi2.cpp

#ifndef _VITERBI_H_
#define _VITERBI_H_

#include "mex.h"
#include <cmath>
#include <limits>


double INF = std::numeric_limits<double>::infinity();

class viterbi2
{
public:
  viterbi2(const int d_row,
          const int d_col,
          const double *d_image,
          const double *d_along_track_slope,
          const double d_along_track_weight,
          const double *d_upper_bounds,
          const double *d_lower_bounds,
          double *d_result) : f_row(d_row),
                              f_col(d_col),
                              f_image(d_image),
                              f_along_track_slope(d_along_track_slope),
                              f_along_track_weight(d_along_track_weight),
                              f_upper_bounds(d_upper_bounds),
                              f_lower_bounds(d_lower_bounds),
                              f_result(d_result)
  {
    find_path();
  }

  // VARIABLES
  const int f_row, f_col;
  const double *f_image, *f_along_track_slope, f_along_track_weight, *f_upper_bounds, *f_lower_bounds;
  double *f_result;

  // METHODS
  int calculate_best(double *path_prob);
  double *find_path(void);
  void viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index);

  // Compute square value
  template <class T>
  inline T sqr(T x) { return x * x; }

  int encode(int x, int y) { return x * f_row + y; }
  int vic_encode(int row, int col) { return f_row * col + row; }

  // THE CODE BELOW THIS POINT WAS TAKEN FROM DAVID CRANDALL
  // With minor adjustments by Reece Mathews (transition_weight)
  // Distance transform
  // -- Every index from d1 to d2 will be set in dst and dst_ind
  // -- dst will contain the minimum value for that destination
  // -- dst_ind will contain the minimum source index for that destination
  void dt(const double *src, double *dst, double *dst_ind, int s1, int s2,
          int d1, int d2, double transition_weight, int off = 0)
  {

    int d = (d1 + d2) >> 1, s = ((s1 + s2) >> 1); // Find the midpoint of the destination
    for (int p = s1; p <= s2; p++)
    { // Search through all the sources and find the minimum
      if (src[p] + sqr(p - d - off) * transition_weight < src[s] + sqr(s - d - off) * transition_weight)
      {
        s = p;
      }
    }
    dst[d] = src[s] + sqr(s - d - off) * transition_weight; // Minimum value to the midpoint
    dst_ind[d] = s;                                         // Minimum source index for the midpoint
    if (d2 >= d + 1)
    { // Recursive call, binary search (top half of destinations)
      dt(src, dst, dst_ind, s, s2, d + 1, d2, transition_weight, off);
    }
    if (d - 1 >= d1)
    { // Recursive call, binary search (bottom half of destinations)
      dt(src, dst, dst_ind, s1, s, d1, d - 1, transition_weight, off);
    }
  }
  // END CODE FROM DAVID CRANDALL
};
#endif

// TODO[reece]: dt must update bounds of next column to refer to rows within bounds of previous column
//              Does not appear to be happening. Due to dt's search-space optimization? Does it not update entire search space?

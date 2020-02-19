// viterbi.h
//
// Layer-tracking program based on the Viterbi algorithm
//  for MUSIC-processed 2D and 3D data
//
// Adapted from original code by Mingze Xu, David Crandall, and John Paden
//
// Author: Victor Berger
// See also: viterbi.cpp

#ifndef _VITERBI_H_
#define _VITERBI_H_

#include "mex.h"
#include <cmath>
#include <limits>

// Default mu (mean) of smoothness
const int    MID   	     = 33;
// Default sigma (variance) of smoothness
const double SIGMA 		 = 24; 
// Scaling for smoothness constraint
const double SCALE 		 = 12;
// Relative weight of extra GT
const double EGT_WEIGHT  = 10;
// Relative weight of intersection with expected surface multiple bin
// TODO[reece]: Allow multiple bin params to be passed into viterbi
const double MULTIPLE_BIN_WEIGHT = 100;
// Maximum multiple number which has any effect on cost
const double MULTIPLE_MAX_NUM = 5;
// Maximum travel time distance from closest multiple which has any effect on cost
const double MULTIPLE_MAX_DIST = 10;
// Mass conservation weight
const double MC_WEIGHT   = 10; 
// Large cost
const double LARGE       = 1000000000; 

class viterbi {
public:
	viterbi( const int d_row,          
             const int d_col,
             const double *d_image,         
             const int *d_sgt,
             const double *d_mask,
		     const double d_img_mag_weight,            
		     const double *d_smooth_slope, 
             const ptrdiff_t *d_bounds,
             const int d_num_extra_tr,
		     const double *d_egt_x,         
             const double *d_egt_y,
             const double *d_gt_weights,
             const double *d_mask_dist,
             const double *d_costmatrix,
             const int d_costmatrix_X,
             const int d_costmatrix_Y,
             const double d_transition_weight,
             const int d_zero_bin,
             double *d_result
	) : 
	f_row(d_row),                   
	f_col(d_col),
	f_image(d_image),
	f_sgt(d_sgt),
	f_mask(d_mask),
	f_img_mag_weight(d_img_mag_weight),
	f_smooth_slope(d_smooth_slope), 
	f_bounds(d_bounds),
	f_num_extra_tr(d_num_extra_tr),
	f_egt_x(d_egt_x),        
	f_egt_y(d_egt_y),
	f_gt_weights(d_gt_weights),
    f_mask_dist(d_mask_dist),
    f_costmatrix(d_costmatrix),
    f_costmatrix_X(d_costmatrix_X),
    f_costmatrix_Y(d_costmatrix_Y),
    f_transition_weight(d_transition_weight),
    f_zero_bin(d_zero_bin),
    f_result(d_result)
	{  
		find_path();
	}

	// VARIABLES
	const int f_row, f_col, *f_sgt, f_num_extra_tr, f_costmatrix_X, f_costmatrix_Y, f_zero_bin;
	const double *f_image, *f_mask, f_img_mag_weight, *f_smooth_slope,
                 *f_egt_x, *f_egt_y, *f_gt_weights, *f_mask_dist,  
                 *f_costmatrix, f_transition_weight;
	const ptrdiff_t *f_bounds;
	double *f_result, *f_cost;

	int num_col_vis, start_col, end_col;
  double multiple_cost_base;

	// METHODS
	int calculate_best(double *path_prob);
	double unary_cost(int x, int y), *find_path(void);
	void viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index);

    // Dynamic smoothness
	double norm_pdf(int x, double mu = MID, double sigma = SIGMA, double scale = SCALE) 
    {
        if(sigma != std::numeric_limits<double>::infinity()) {
			return scale * (1.0/(sigma*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/sigma));
        }
		return scale;
    }
    
    // Compute square value
    template <class T> inline T sqr(T x) { return x*x; }
    
	int encode(int x, int y) { return x * f_row + y; }
	 
    int vic_encode(int row, int col) { return f_row * col + row; }

    // THE CODE BELOW THIS POINT WAS TAKEN FROM DAVID CRANDALL
    // Distance transform
    // -- Every index from d1 to d2 will be set in dst and dst_ind
    // -- dst will contain the minimum value for that destination
    // -- dst_ind will contain the minimum source index for that destination
    void dt(const double *src, double *dst, double *dst_ind, int s1, int s2,
            int d1, int d2, double transition_weight, int off = 0) 
    {
        
        int d = (d1 + d2) >> 1, s = s1; // Find the midpoint of the destination
        for (int p = s1; p <= s2; p++) { // Search through all the sources and find the minimum
            if (src[s] + sqr(s-d-off) * transition_weight > src[p] + sqr(p-d-off) * transition_weight) {
                s = p;
            }
        }
        dst[d] = src[s] + sqr(s-d-off) * transition_weight; // Minimum value to the midpoint
        dst_ind[d] = s; // Minimum source index for the midpoint
        if(d2 >= d + 1) { // Recursive call, binary search (top half of destinations)
            dt(src, dst, dst_ind, s, s2, d+1, d2, transition_weight, off);
        }
        if(d-1 >= d1) { // Recursive call, binary search (bottom half of destinations)
            dt(src, dst, dst_ind, s1, s, d1, d-1, transition_weight, off);
        }
    }
    void dt_1d(const double *f, double transition_weight, double *result, double *dst_ind, 
               int beg, int end, int off = 0) {
        dt(f, result, dst_ind, beg, end-1, beg, end-1, transition_weight, off);
    }
    // END CODE FROM DAVID CRANDALL	
};
#endif

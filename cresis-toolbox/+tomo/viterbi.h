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
const double EGT_WEIGHT  = 10;  // TODO[reece]: is this supposed to be positive?
                                //              appears to divide added cost for distance from egt
                                //              Repulsion is positive and increases cost
// Relative weight of intersection with expected surface multiple bin
// TODO[reece]: Allow multiple bin params to be passed into viterbi
const double MULTIPLE_BIN_WEIGHT = 100;
// Maximum multiple number which has any effect on cost
const double MULTIPLE_MAX_NUM = 5;
// Maximum travel time distance from closest multiple which has any effect on cost
const double MULTIPLE_MAX_DIST = 10;
// Icemask proximity scan threshold
const double ICE_BIN_THR = 3;
// Repulsion from surface
const double REPULSION   = 150000; 
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
             const int d_bgt,               
             const double *d_mask,
		     const double *d_mu,            
             const double *d_sigma,
		     const int d_mid,              
             const double d_egt_weight,
		     const double d_smooth_weight,  
             const double d_smooth_var,
		     const double *d_smooth_slope, 
             const ptrdiff_t *d_bounds,
		     const size_t d_ms,             
             const int d_num_extra_tr,
		     const double *d_egt_x,         
             const double *d_egt_y,
             const double *d_weight_points,
		     const double d_repulsion,      
             const double d_ice_bin_thr,
             const double *d_mask_dist,
             const double *d_costmatrix,
             const int d_costmatrix_X,
             const int d_costmatrix_Y,
             const double d_scale,
             const int d_plane_bin,
             double *d_result
	) : 
	f_row(d_row),                   
	f_col(d_col),
	f_image(d_image),
	f_sgt(d_sgt),
	f_bgt(d_bgt),
	f_mask(d_mask),
	f_mu(d_mu),
	f_sigma(d_sigma),
	f_mid(d_mid),    
	f_egt_weight(d_egt_weight),
	f_smooth_weight(d_smooth_weight),
	f_smooth_var(d_smooth_var),
	f_smooth_slope(d_smooth_slope), 
	f_bounds(d_bounds),
	f_ms(d_ms),         
	f_num_extra_tr(d_num_extra_tr),
	f_egt_x(d_egt_x),        
	f_egt_y(d_egt_y),
	f_weight_points(d_weight_points),
	f_repulsion(d_repulsion),    
	f_ice_bin_thr(d_ice_bin_thr), 
    f_mask_dist(d_mask_dist),
    f_costmatrix(d_costmatrix),
    f_costmatrix_X(d_costmatrix_X),
    f_costmatrix_Y(d_costmatrix_Y),
    f_scale(d_scale),
    f_plane_bin(d_plane_bin),
    f_result(d_result)
	{  
		find_path();
	}

	// VARIABLES
	const int f_row, f_col, f_mid, f_bgt, *f_sgt, f_num_extra_tr, f_costmatrix_X, f_costmatrix_Y, f_plane_bin;
	const double *f_image, *f_mask, *f_mu, *f_sigma, f_egt_weight, *f_smooth_slope,
                 *f_egt_x, *f_egt_y, *f_weight_points, f_smooth_weight,
                 f_smooth_var, f_repulsion, f_ice_bin_thr, *f_mask_dist,  
                 *f_costmatrix, f_scale;
	const ptrdiff_t *f_bounds;
	const size_t f_ms;
	double *f_result, *f_cost;

	int depth, num_col_vis, t, start_col, end_col;
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
	 
    int vic_encode(int row, int col) { return depth * col + row; }

    // THE CODE BELOW THIS POINT WAS TAKEN FROM DAVID CRANDALL
    // Distance transform
    // -- Every index from d1 to d2 will be set in dst and dst_ind
    // -- dst will contain the minimum value for that destination
    // -- dst_ind will contain the minimum source index for that destination
    void dt(const double *src, double *dst, double *dst_ind, int s1, int s2,
            int d1, int d2, double scale, int off = 0) 
    {
        
        int d = (d1 + d2) >> 1, s = s1; // Find the midpoint of the destination
        for (int p = s1; p <= s2; p++) { // Search through all the sources and find the minimum
            if (src[s] + sqr(s-d-off) * scale > src[p] + sqr(p-d-off) * scale) {
                s = p;
            }
        }
        dst[d] = src[s] + sqr(s-d-off) * scale; // Minimum value to the midpoint
        dst_ind[d] = s; // Minimum source index for the midpoint
        if(d2 >= d + 1) { // Recursive call, binary search (top half of destinations)
            dt(src, dst, dst_ind, s, s2, d+1, d2, scale, off);
        }
        if(d-1 >= d1) { // Recursive call, binary search (bottom half of destinations)
            dt(src, dst, dst_ind, s1, s, d1, d-1, scale, off);
        }
    }
    void dt_1d(const double *f, double scale, double *result, double *dst_ind, 
               int beg, int end, int off = 0) {
        dt(f, result, dst_ind, beg, end-1, beg, end-1, scale, off);
    }
    // END CODE FROM DAVID CRANDALL	
};
#endif

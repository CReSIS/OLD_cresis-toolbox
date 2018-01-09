// viterbi_lib.h
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
const int    DEF_MID   	     = 33;
const int    DEF_SIGMA 		 = 24;
const int    DEF_SCALE 		 = 12; // scaling for smoothness constraint
const int    DEF_EGT_WEIGHT  = 10; // relative weight of extra GT
const int 	 DEF_ICE_BIN_THR = 3;  // icemask proximity scan threshold
const int    DEF_REPULSION   = 150000; // repulsion from surface
const double LARGE           = 1000000000;

class detect
{
public:
	detect( const int d_row,               const int d_col,
		const double *d_in_data,       const int *d_slayer,
		const int d_blayer,            const double *d_mask,
		const double *d_mu,            const double *d_sigma,
		const int d_mid,               const double d_egt_weight,
		const double d_smooth_weight,  const double d_smooth_var,
		const double *d_smooth_slope,  const ptrdiff_t *d_bounds,
		const size_t d_ms,             const int d_num_extra_tr,
		const double *d_extra_truth_x, const double *d_extra_truth_y,
		double *d_result,			   const double *d_weight_points,
		const double d_repulsion,      const double d_ice_bin_thr 
	)
	: f_row(d_row),                     f_col(d_col),
	f_in_data(d_in_data),             f_slayer(d_slayer),
	f_blayer(d_blayer),               f_mask(d_mask),
	f_mu(d_mu),                       f_sigma(d_sigma),
	f_mid(d_mid),                     f_egt_weight(d_egt_weight),
	f_smooth_weight(d_smooth_weight), f_smooth_var(d_smooth_var),
	f_smooth_slope(d_smooth_slope),   f_bounds(d_bounds),
	f_ms(d_ms),                       f_num_extra_tr(d_num_extra_tr),
	f_extra_truth_x(d_extra_truth_x), f_extra_truth_y(d_extra_truth_y),
	f_result(d_result),               f_weight_points(d_weight_points),
	f_repulsion(d_repulsion),         f_ice_bin_thr(d_ice_bin_thr) 
	{
		find_path();
	}
	// VARIABLES
	const int f_row, f_col, f_mid, f_blayer,
		*f_slayer, f_num_extra_tr;
	const double *f_in_data, *f_mask, *f_mu, *f_sigma, f_egt_weight, *f_smooth_slope,
                 *f_extra_truth_x, *f_extra_truth_y, *f_weight_points, f_smooth_weight,
                 f_smooth_var, f_repulsion, f_ice_bin_thr; 
	const ptrdiff_t *f_bounds;
	const size_t f_ms;
	double *f_result;
	int depth, num_col_vis, t, start_col, end_col;

	// METHODS
	int calculate_best(double *path_prob);
	double unary_cost(int x, int y), *find_path(void);
	void viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index, int path_index),
	reset_arrays(double *path_prob, double *path_prob_next, double *index);

	template <class T>
	inline T sqr(T x) { return x*x; }
	
	template <class T>
	inline T sqr2(T x) 
    {
        x = abs(x);
        if (x>3) 
          return pow(x*x,2)*0.11111;
        return x*x;
        
    }

	double norm_pdf(int x, double mu = DEF_MID,
		double sigma = DEF_SIGMA, double scale = DEF_SCALE)
		{
			if(sigma != std::numeric_limits<double>::infinity())
				return scale * (1.0/(sigma*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/sigma));
			return scale;
		}

		int encode(int x, int y)
		{
			return x * f_row + y;
		}

		int vic_encode(int row, int col)
		{
			return depth * col + row;
		}

		//  Distance transform calculation, adapted from David Crandall
		//  Every index from d1 to d2 will be set in dst and dst_ind
		//  dst will contain the minimum value for that destination
		//  dst_ind will contain the minimum source index for that destination
		void dt(const double *src, double *dst, double *dst_ind, int s1, int s2,
			int d1, int d2, double scale, int off = 0)
			{

				int d = (d1 + d2) >> 1, s = s1; // Find the midpoint of the destination
				for (int p = s1; p <= s2; p++) // Search through all the sources and find the minimum
// 					if (src[s] + sqr2(s-d-off) * scale > src[p] + sqr2(p-d-off) * scale)
                    if (src[s] + sqr(s-d-off) * scale > src[p] + sqr(p-d-off) * scale)
						s = p;

// 				dst[d] = src[s] + sqr2(s-d-off) * scale; 
                dst[d] = src[s] + sqr(s-d-off) * scale; // Minimum value to the midpoint

				dst_ind[d] = s; // Minimum source index for the midpoint

				if(d2 >= d + 1) // Recursive call, binary search (top half of destinations)
					dt(src, dst, dst_ind, s, s2, d+1, d2, scale, off);

				if(d-1 >= d1) // Recursive call, binary search (bottom half of destinations)
					dt(src, dst, dst_ind, s1, s, d1, d-1, scale, off);
			}

			void dt_1d(const double *f, double scale, double *result,
				double *dst_ind, int beg, int end, int off = 0)
				{
					dt(f, result, dst_ind, beg, end-1, beg, end-1, scale, off);
				}
			};
			#endif

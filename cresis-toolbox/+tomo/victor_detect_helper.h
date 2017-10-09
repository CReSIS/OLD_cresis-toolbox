#ifndef _VICTOR_DETECT_HELPER_H_
#define _VICTOR_DETECT_HELPER_H_

#include "mex.h"
#include <cmath>
#include <string.h>

const int    def_mid   = 33;
const int    def_sigma = 24;
const int    def_scale = 5;
const double large     = 1000000000;

class detect
{
	public:
		detect(const int d_row,              const int d_col,
				const double *d_in_data,      const int *d_slayer,
				const int d_blayer,           const double *d_mask,
				const double *d_mu,           const double *d_sigma,
				const int d_mid,              const int d_egt_weight,
				const int d_smooth_weight,    const int d_smooth_var,
				const double *d_smooth_slope, const ptrdiff_t *d_bounds,
				const size_t d_ms,            const int d_num_extra_tr,
				const double *d_extra_truth_x,const double *d_extra_truth_y,
				double* d_result)
			: f_row(d_row),                    f_col(d_col),
			  f_in_data(d_in_data),            f_slayer(d_slayer),
			  f_blayer(d_blayer),              f_mask(d_mask),
			  f_mu(d_mu),                      f_sigma(d_sigma),
			  f_mid(d_mid),                    f_egt_weight(d_egt_weight),
			  f_smooth_weight(d_smooth_weight),f_smooth_var(d_smooth_var),
			  f_smooth_slope(d_smooth_slope),  f_bounds(d_bounds),
			  f_ms(d_ms),                      f_num_extra_tr(d_num_extra_tr),
			  f_extra_truth_x(d_extra_truth_x),f_extra_truth_y(d_extra_truth_y),
			  f_result(d_result)
	{
		find_path();
	}
		// VARIABLES
		const int f_row, f_col, f_mid, f_smooth_weight, f_egt_weight,
		      		f_smooth_var, f_blayer, *f_slayer, f_num_extra_tr;
		const double *f_in_data, *f_mask, *f_mu, *f_sigma,
		      			 *f_smooth_slope,*f_extra_truth_x, *f_extra_truth_y;
		const ptrdiff_t *f_bounds;
		const size_t f_ms;
		double *f_result;

		// METHODS
		int num_nonzero(double *arr, int length);
		double unary_cost(int x, int y), *find_path(void);
		void ac(double* s, double* d, int n);

		template <class T>
			T sqr(T x) { return x*x; }

		double norm_pdf(int x, double mu = def_mid,
				double sigma = def_sigma, double scale = def_scale)
		{
			return scale * (1.0/(sigma*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/sigma));
		}

		int encode(int x, int y)
		{
			return ((x * f_row) + y);
		}

		// Distance transform calculation, adapted from David Crandall
		//  Every index from d1 to d2 will be set in dst and dst_ind
		//  dst will contain the minimum value for that destination
		//  dst_ind will contain the minimum source index for that destination
		void dt(const double *src, double *dst, double *dst_ind, int s1, int s2,
						int d1, int d2, double scale, int off = 0)
		{
			int d = (d1 + d2) / 2, s = s1; // Find the midpoint of the destination
			for (int p = s1; p <= s2; p++) // Search through all the sources and find the minimum
				if (src[s] + sqr(s-d-off) * scale > src[p] + sqr(p-d-off) * scale)
					s = p;

			dst[d] = src[s] + sqr(s-d-off) * scale; // Minimum value to the midpoint
			dst_ind[d] = s; // Minimum source index for the midpoint

			if(d2 >= d+1) // Recursive call, binary search (top half of destinations)
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

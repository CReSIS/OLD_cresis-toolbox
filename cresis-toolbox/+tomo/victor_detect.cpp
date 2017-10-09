#include "victor_detect_helper.h"
#include <algorithm>

const bool debug = 0;

// Entry point from Matlab
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
								 const mxArray *prhs[])
{
	if(nrhs != 13)
		mexErrMsgTxt("Usage: [labels] = detect(input_img, surface_gt, bottom_gt, extra_gt, ice_mask, mean, var, mid, egt_weight, smooth_weight, smooth_var, smooth_slope, bounds)\n");

	if(debug)
		mexPrintf("\n$$$ DEBUG MODE ON $$$\nRunning new Viterbi detection routine...");

	const int    _row             = mxGetM(prhs[0]);
	const int    _col             = mxGetN(prhs[0]);
	const double *_in_data        = mxGetPr(prhs[0]);
	const double *_bins           = mxGetPr(prhs[1]);
	const double *t_truth         = mxGetPr(prhs[2]);
	const double *t_extra_truth   = mxGetPr(prhs[3]);
	const double *_mask           = mxGetPr(prhs[4]);
	const double *_mu             = mxGetPr(prhs[5]); // mean
	const double *_sigma          = mxGetPr(prhs[6]); // variance
	const double *t_mid           = mxGetPr(prhs[7]);
	const double *t_egt_weight    = mxGetPr(prhs[8]);
	const double *t_smooth_weight = mxGetPr(prhs[9]);
	const double *t_smooth_var    = mxGetPr(prhs[10]);
	const double *_smooth_slope   = mxGetPr(prhs[11]);
	const size_t _ms              = mxGetNumberOfElements(prhs[5]);
	const int    _num_extra_tr    = mxGetNumberOfElements(prhs[3]);

	if(_num_extra_tr % 2 != 0)
		mexErrMsgTxt("Extra truth vector not properly defined.");

	// Bounds
	ptrdiff_t _bounds[2];
	if (nrhs >= 13 && mxGetNumberOfElements(prhs[12]))
	{
		if (!mxIsInt64(prhs[12]))
			mexErrMsgTxt("Usage: bounds must be type int64");
		if (mxGetNumberOfElements(prhs[12]) != 2)
			mexErrMsgTxt("Usage: bounds must be a 2 element vector");

		ptrdiff_t *tmp = (ptrdiff_t*)mxGetPr(prhs[12]);
		_bounds[0] = tmp[0];
		_bounds[1] = tmp[1];
		if (_bounds[0] < 0)
			_bounds[0] = 0;
		if (_bounds[1] < 0)
			_bounds[1] = _col - 1;
		if (_bounds[0] >= _col)
			mexErrMsgTxt("Usage: bounds[0] < size(input,2)");
		if (_bounds[1] >= _col)
			mexErrMsgTxt("Usage: bounds[1] < size(input,2)");
		if (_bounds[1] < _bounds[0])
			mexErrMsgTxt("Usage: bounds[1] must be greater than bounds[0]");
	}
	else
	{
		// Default setting is to process all columns
		_bounds[0] = 0;
		_bounds[1] = _col - 1;
	}

	// Initialize surface layer array
	int _slayer[_col];
	for(int k = 0; k < _col; ++k)
		_slayer[k] = (int)_bins[k];

	// Initialize variables to default values if temporary values not set
	const int _smooth_weight = (t_smooth_weight[0] < 0) ? def_scale : (int)t_smooth_weight[0];
	const int _smooth_var    = (t_smooth_var[0] < 0) ? def_sigma : (int)t_smooth_var[0];
	const int _mid           = (t_mid[0] < 0) ? (_col / 2) : (int)t_mid[0];
	const int _blayer        = (t_truth[0] > 0) ? (t_truth[0]) : _slayer[(int)t_mid[0]] + 50;
	const int _egt_weight    = (int)t_egt_weight[0];

	double _extra_truth_x[(_num_extra_tr / 2)],
	       _extra_truth_y[(_num_extra_tr / 2)];
	for (int p = 0; p < (_num_extra_tr / 2); ++p)
	{
		_extra_truth_x[p] = t_extra_truth[p * 2];
		_extra_truth_y[p] = t_extra_truth[(p * 2) + 1];
	}

	plhs[0] = mxCreateDoubleMatrix(1, _col, mxREAL);
	double *_result = mxGetPr(plhs[0]);

	// Call detect class constructor
	detect obj(_row, _col, _in_data, _slayer, _blayer, _mask,
			_mu, _sigma, _mid, _egt_weight, _smooth_weight,
			_smooth_var, _smooth_slope, _bounds, _ms,
			_num_extra_tr, _extra_truth_x, _extra_truth_y,
			_result);

}

//  Used to define unary cost of target of position x, y
double detect::unary_cost(int x, int y)
{
	if (debug)
		mexPrintf("\nRunning unary cost function...");

	int t = (f_ms - 1) / 2;

	// Set cost to large if bottom is above surface
	if((f_slayer[x] > t) && (f_slayer[x] < f_row - t) && ((y + t + 1) < f_slayer[x]))
		return large;

	// Set cost of center point to large if far from center ground truth
	if((x == f_mid) && (y + t < f_blayer || y + t > f_blayer + 500))
		return large;

	double cost = 0;

	// Increase cost if point is far from extra ground truth
	for(int f = 0; f < (f_num_extra_tr / 2); ++f)
		if(f_extra_truth_x[f] == x)
		{
			cost += 2 * sqr(abs((int)f_extra_truth_y[f] - (int)(t + y)) / f_egt_weight);
			break;
		}

	// Increased cost if bottom is close to surface
	if(abs(y + t) - f_slayer[x] < 10)
		cost += 100 - 10 * abs((int)(y + t) - (int)f_slayer[x]);

	// Add template quadratic distance cost
	for(int k = 0; k < f_ms; ++k)
		cost += sqr(f_in_data[encode(x, y + k)] - f_mu[k]) / f_sigma[k];

	if(debug)
		mexPrintf("\nUnary cost of target (x = %d, y = %d) -> %f", x, y, cost);

	return cost;
}

// Returns Viterbi solution of optimal path
double* detect::find_path(void)
{
	if (debug)
		mexPrintf("\nRunning path finding function...");

	int       t           = (f_ms - 1) / 2;
	int       depth       = f_row - f_ms;
	const int start_col   = f_bounds[0];
	const int end_col     = f_bounds[1];
	const int num_col_vis = end_col - start_col + 1;

	int loop = 0, next = 1, num;

	double pathA[depth][num_col_vis], pathB[depth][num_col_vis];
	double path_prob[2][depth], index[depth];

	// Initialize path and path_prob matrices
	for(int k = 0; k < depth; ++k)
	{
		index[k]        = 0;
		path_prob[0][k] = 0;
		path_prob[1][k] = 0;
		for(int l = 0; l < num_col_vis; ++l)
		{
			pathA[k][l] = -1;
			pathB[k][l] = -1;
		}
	}

	// Initial column (leftmost)
	if(f_mask[start_col] == 0 && f_slayer[start_col] > t)
		for(int row = 0; row < depth; ++row)
		{
			pathA[row][0] = (double)row;
			if ((row + t) != f_slayer[start_col])
				path_prob[0][row] = large;
		}
	else // Set unary costs for all rows of initial column
		for(int row = 0; row < depth; ++row)
		{
			pathA[row][0] = (double)row;
			path_prob[0][row] = unary_cost(start_col, row);
		}

	// Other columns
	double norm = 0;
	for(int col = start_col +1; col <=end_col; ++col)
	{
		norm = norm_pdf(col, (double)f_mid, f_smooth_var, f_smooth_weight);

		// Calculate cost of each path with distance transform function
		dt_1d(path_prob[loop % 2], norm, path_prob[next % 2], index, 0, depth, f_smooth_slope[col - 1]);

		if(f_mask[col] == 0 && f_slayer[col] > t)
			for(int row = 0; row < depth; ++row)
			{
				if(next % 2 != 0)
				{
					num = num_nonzero(pathA[(int)index[row]], num_col_vis);
					ac(pathA[(int)index[row]], pathB[row], num);
					pathB[row][num] = (double)row;
				}
				else
				{
					num = num_nonzero(pathB[(int)index[row]], num_col_vis);
					ac(pathB[(int)index[row]], pathA[row], num);
					pathA[row][num] = (double)row;
				}

				if(row + t != f_slayer[col])
					path_prob[next % 2][row] = large;
			}
		else
			for(int row = 0; row < depth; ++row)
			{
				if(next % 2 != 0)
				{
					num = num_nonzero(pathA[(int)index[row]], num_col_vis);
					ac(pathA[(int)index[row]], pathB[row], num);
					pathB[row][num] = (double)row;
				}
				else
				{
					num = num_nonzero(pathB[(int)index[row]], num_col_vis);
					ac(pathB[(int)index[row]], pathA[row], num);
					pathA[row][num] = (double)row;
				}
				path_prob[next % 2][row] += unary_cost(col, row);
			}
		++loop; ++next;
	}

	// Find minimum cost final result
	int pos, min = large;
	for(int k = 0; k < depth; ++k)
		if(path_prob[loop % 2][k] < min)
		{
			min = path_prob[loop % 2][k];
			pos = k;
		}

	// Set 'result' array
	for(int v = 0; v < f_col; ++v)
		f_result[v] = (v < f_bounds[0] || v > end_col) ? 0 : ((loop % 2 == 0) ? (pathA[pos][v - f_bounds[0]] + t) : (pathB[pos][v - f_bounds[0]] + t));

	return f_result;
}

// Count number of non-zero values starting from the last
int detect::num_nonzero(double *ar, int l)
{
	for(int k = l - 1; k >= 0; --k)
		if(ar[k] != -1)
			return (k + 1);
	return 0;
}

// Array copy
void detect::ac(double* s, double* d, int n)
{
	std::copy(s, s + n, d);
}

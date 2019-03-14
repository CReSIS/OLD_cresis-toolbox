// viterbi.cpp
//
// Layer-tracking program based on the Viterbi algorithm
// 
// Adapted from original code by Mingze Xu, David Crandall, and John Paden
//
// Author: Victor Berger
// See also: viterbi_lib.h

#include "viterbi_lib.h"
const bool debug = 0;

// Entry point from Matlab
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if(nrhs != 14 && nrhs != 15)
	mexErrMsgTxt("Usage: [labels] = viterbi(input_img, surface_gt, bottom_gt, extra_gt, ice_mask, mean, var, mid, egt_weight, smooth_weight, smooth_var, smooth_slope, bounds, viterbi_weight, repulsion, ice_bin_thr)\n");

  if(debug)
	mexPrintf("\n$$$ DEBUG MODE ON $$$\n");

  const int    _row             = mxGetM(prhs[0]);
  const int    _col             = mxGetN(prhs[0]);
  const double *_in_data        = mxGetPr(prhs[0]);
  const double *_surf_tr        = mxGetPr(prhs[1]);
  const double *t_bott_tr       = mxGetPr(prhs[2]);
  const double *t_extra_truth   = mxGetPr(prhs[3]);
  const double *_mask           = mxGetPr(prhs[4]);
  const double *_mu             = mxGetPr(prhs[5]); // mean
  const double *_sigma          = mxGetPr(prhs[6]); // variance
  const double *t_egt_weight    = mxGetPr(prhs[7]);
  const double *t_smooth_weight = mxGetPr(prhs[8]);
  const double *t_smooth_var    = mxGetPr(prhs[9]);
  const double *_smooth_slope   = mxGetPr(prhs[10]);
  const double *_weight_points  = mxGetPr(prhs[12]);
  const double *t_repulsion     = mxGetPr(prhs[13]);
  const double *t_ice_bin_thr   = mxGetPr(prhs[14]);
  const int    _num_extra_tr    = mxGetNumberOfElements(prhs[3]);
  const size_t _ms              = mxGetNumberOfElements(prhs[5]);

  if(_num_extra_tr % 2 != 0)
	mexErrMsgTxt("\nERROR:\nExtra truth vector not properly defined.\n");

  // Bounds
  ptrdiff_t _bounds[2];
  if(nrhs >= 13 && mxGetNumberOfElements(prhs[11]))
  {
    if(!mxIsInt64(prhs[11]))
		mexErrMsgTxt("Usage: bounds must be type int64");
    if(mxGetNumberOfElements(prhs[11]) != 2)
		mexErrMsgTxt("Usage: bounds must be a 2 element vector");

    ptrdiff_t *tmp = (ptrdiff_t*)mxGetPr(prhs[11]);
    _bounds[0] = tmp[0];
    _bounds[1] = tmp[1];
    if(_bounds[0] < 0)
		_bounds[0] = 0;
    if(_bounds[1] < 0)
		_bounds[1] = _col;
    if(_bounds[0] > _col)
		mexErrMsgTxt("Usage: bounds[0] <= size(input,2)");
    if(_bounds[1] > _col)
		mexErrMsgTxt("Usage: bounds[1] <= size(input,2)");
    if(_bounds[1] < _bounds[0])
		mexErrMsgTxt("Usage: bounds[1] must be greater than bounds[0]");

  }
  else
  {
    // Default setting is to process all columns
    _bounds[0] = 0;
    _bounds[1] = _col;
  }

  // Initialize surface layer array
  int _slayer[_col];
  for(int k = 0; k < _col; ++k)
	_slayer[k] = (int)_surf_tr[k];

  // Initialize variables to default values if temporary values not set
  const int _mid              = floor(_col / 2);
  const int _blayer           = ((t_bott_tr) ? (t_bott_tr[0] > 0 ? round(t_bott_tr[0]) : -1) : -1);
  const double _egt_weight    = t_egt_weight && t_egt_weight < 0 ? DEF_EGT_WEIGHT : t_egt_weight[0];
  const double _smooth_weight = t_smooth_weight[0] < 0 ? DEF_SIGMA : t_smooth_weight[0];
  const double _smooth_var    = t_smooth_var[0] < 0 ? DEF_SIGMA : t_smooth_var[0];
  const double _repulsion     = t_repulsion[0] < 0 ? DEF_REPULSION : t_repulsion[0];
  const double _ice_bin_thr   = t_ice_bin_thr[0] < 0 ? DEF_ICE_BIN_THR : t_ice_bin_thr[0];

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
    _result,  _weight_points, _repulsion, _ice_bin_thr);
  }

  //  Used to define unary cost of target of position x, y
  double detect::unary_cost(int x, int y)
  {
    if (debug)
		mexPrintf("\nRunning unary cost function...");

    // Set cost to large if bottom is above surface
    if((f_slayer[x] > t) && (f_slayer[x] < f_row - t) && ((y + t + 1) < f_slayer[x]))
		return LARGE;

    // Set cost of center point to large if far from center ground truth
    if((f_blayer != -1) && (x == f_mid) && (y + t < f_blayer || y + t > f_blayer + 500))
		return LARGE;

    double cost = 0;

    // Increase cost if point is far from extra ground truth
    for(int f = 0; f < (f_num_extra_tr / 2); ++f)
		if(f_extra_truth_x[f] == x && x < f_bounds[1] && x > f_bounds[0])
		{
		  cost += f_weight_points[x] * sqr(abs((int)f_extra_truth_y[f] - (int)(t + y)) / f_egt_weight);
		  break;
		}

    // Reduce cost for bins near a region with negative ice mask
    for(int dist = 0, range = x - f_ice_bin_thr; range <= x + f_ice_bin_thr; ++range)
    {
      if(range == x || abs(y - f_slayer[x]) > 40 || range < f_bounds[0] || range > f_bounds[1])
		continue;
      if(f_mask[range] == 0)
      {
        dist = abs(x - range);
        cost -= 1000 * dist;
      }
    }

    // Increased cost if bottom is close to surface
    if(abs(y + t - f_slayer[x]) < 10)
		cost += f_repulsion * (100 - 10 * abs((int)(y + t) - (int)f_slayer[x]));

    // Image magnitude correlation
    double tmp_cost = 0;
    for (size_t i = 0; i < f_ms; i++)
		tmp_cost  -= 2.5 * (f_in_data[encode(x, y + i)]) * f_mu[i] / f_sigma[i];
    tmp_cost = tmp_cost*12/10 - 60;

    cost += tmp_cost;

    if(debug)
		mexPrintf("\nUnary cost of target (x = %d, y = %d) -> %f", x, y, cost);

    return cost;
  }

  // Returns Viterbi solution of optimal path
  double* detect::find_path(void)
  {
    if (debug)
    mexPrintf("\nRunning path finding function...");

    start_col   = f_bounds[0];
    end_col     = f_bounds[1];
    t           = (f_ms - 1) / 2;
    depth       = f_row - f_ms;
    num_col_vis = end_col - start_col;

    int *path = new int[depth * (num_col_vis + 2)];
    double path_prob[depth], path_prob_next[depth], index[depth];

    for(int k = 0; k < f_col; ++k)
		f_result[k] = 0;

    for(int k = 0; k < depth * (num_col_vis + 2); ++k)
		path[k] = 0;

    reset_arrays(path_prob, path_prob_next, index);
    viterbi_right(path, path_prob, path_prob_next, index, 0);

    int encode;
    int viterbi_index     = calculate_best(path_prob);
    int idx               = end_col;
    f_result[end_col - 1] = f_mask[end_col - 1] == 1 ? viterbi_index + t : f_slayer[end_col - 1];

    // Set result vector
    for(int k = start_col + 1; k <= end_col; ++k)
    {
      encode = vic_encode(viterbi_index, num_col_vis + start_col - k);
      viterbi_index = path[encode];
      f_result[idx - 2] = f_mask[idx - 2] == 1 ? viterbi_index + t : f_slayer[idx - 2];
      --idx;
      if(encode < 0 || idx < 2)
		break;
    }
    delete[] path;
    return f_result;
  }

  // Select path with lowest overall cost
  int detect::calculate_best(double *path_prob)
  {
    double min = LARGE;
    int viterbi_index = 0;
    for(int k = 0; k < depth; ++k)
		if(path_prob[k] < min)
		{
		  min = path_prob[k];
		  viterbi_index = k;
		}
    return viterbi_index;
  }

  // Reset values of given arrays to 0
  void detect::reset_arrays(double *path_prob, double *path_prob_next, double *index)
  {
    for(int k = 0; k < depth; ++k)
    {
      path_prob[k] = 0;
      path_prob_next[k] = 0;
      index[k] = 0;
    }
  }

  // Perform detection to the right
  void detect::viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index, int path_index)
  {
    bool next = 0;
    double norm = 0;
    for(int col = start_col; col <= end_col + 1; ++col)
    {
      if(path_index >= depth * (num_col_vis + 2) ||col >= f_col || col < 0)
		continue;

      norm = norm_pdf(col, (double)f_mid, f_smooth_var, f_smooth_weight);

      if(!next)
		dt_1d(path_prob, norm, path_prob_next, index, 0, depth, f_smooth_slope[col-1]);
      else
		dt_1d(path_prob_next, norm, path_prob, index, 0, depth, f_smooth_slope[col-1]);

      if(f_mask[col] == 0 && f_slayer[col] > t)
		  for(int row = 0; row < depth; ++row)
		  {
			path[path_index] = index[row];
			if(row + t != f_slayer[col])
			{
			  if(!next) path_prob_next[row] += LARGE;
			  else      path_prob[row] += LARGE;
			}
			++path_index;
		  }
      else for(int row = 0; row < depth; ++row)
      {
        path[path_index] = index[row];
        if(!next) path_prob_next[row] += unary_cost(col, row);
        else      path_prob[row] += unary_cost(col, row);
        ++path_index;
      }
      next = !next;
    }
  }

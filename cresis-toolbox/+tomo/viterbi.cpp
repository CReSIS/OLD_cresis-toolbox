// viterbi.cpp
//
// Layer-tracking program based on the Viterbi algorithm
//  for MUSIC-processed 2D and 3D data
//
// Authors: Victor Berger and John Paden
//           Center for Remote Sensing of Ice Sheets
//           2017-2018
//          Adapted from original code by Mingze Xu and David Crandall
//
// Changes to cost function (shifted exponential decay): Victor Berger and John Paden 2018
// Changes to cost function (geostatistical analysis): Victor Berger and John Paden 2019
// Various changes throughout (purge dead code, fix bugs, add more weights): Reece Mathews 2020 
//
// See also: viterbi.h
//
// mex -v -largeArrayDims viterbi.cpp


#include "viterbi.h"

//  Used to define unary cost of target at position x, y
double viterbi::unary_cost(int x, int y)
{

  // Merge layers when no ice exists
  if (f_mask[x] == 0 && y != f_sgt[x])
  {
    return LARGE;
  }

  // Set cost to large if bottom is above surface
  if (y < f_sgt[x])
  {
    return LARGE;
  }

  double cost = 0;

  // Increase cost if far from extra ground truth
  for (int f = 0; f < (f_num_extra_tr / 2); ++f)
  {
    if (f_egt_x[f] == x)
    {
      mexPrintf("x:%d y:%d cutoff:%d dist:%d weight:%d \n", x, y, f_gt_cutoffs[x], abs((int)f_egt_y[f] - (int)y), f_gt_weights[x]);
      if (f_gt_cutoffs[x] > -1 && abs((int)f_egt_y[f] - (int)y) > f_gt_cutoffs[x])
      {
        // Must be within cutoff range
        return LARGE;
      }
      else
      {
        cost += f_gt_weights[x] * 10 * sqr(((int)f_egt_y[f] - (int)y));
        break;
      }
    }
  }

  // Increase cost if near surface or surface multiple bin
  const int travel_time = f_sgt[x] - f_zero_bin; // Between multiples
  int dist_to_surf = abs(y - f_sgt[x]);
  int multiple_bin = travel_time == 0 ? -1 : (y - f_sgt[x]) / travel_time - 1;
  int dist_to_bin = travel_time == 0 ? dist_to_surf : dist_to_surf % travel_time;
  // If closer to next multiple, use that distance instead
  if (dist_to_bin > travel_time / 2 && multiple_bin >= 0)
  {
    dist_to_bin = travel_time - dist_to_bin;
    multiple_bin++;
  }
  double current_mult_weight = f_mult_weight;
  if (multiple_bin < 0)
  {
    // multiple bin -1 is surface, others are multiples
    multiple_bin = 0;
    current_mult_weight = f_surf_weight;
  }

  double multiple_cost = current_mult_weight;
  multiple_cost *= pow(f_mult_weight_decay, multiple_bin);
  multiple_cost *= pow(f_mult_weight_local_decay, dist_to_bin);

  multiple_cost = multiple_cost < 0 ? 0 : multiple_cost;
  cost += multiple_cost;

  // Distance from nearest ice-mask
  // Probabilistic measurement
  int DIM = (std::isinf(f_mask_dist[x]) || f_mask_dist[x] >= f_costmatrix_Y) ? f_costmatrix_Y - 1 : f_mask_dist[x];
  int matrix_idx = f_costmatrix_X * DIM + y - f_sgt[x];
  int matrix_final_idx = f_costmatrix_X * f_costmatrix_Y - 1;
  if (matrix_idx >= 0 && matrix_idx <= matrix_final_idx)
  {
    // Index within bounds of costmatrix
    cost += f_costmatrix[matrix_idx];
  }
  else if (f_costmatrix_X > 0)
  {
    // costmatrix provided but index outside bounds
    if (matrix_idx < 0)
    {
      // Use first cost in matrix when index is before matrix start
      cost += f_costmatrix[0];
    }
    else
    {
      // Use last cost in matrix when index is after matrix end
      cost += f_costmatrix[matrix_final_idx];
    }
  }
  // Image magnitude
  cost -= f_image[encode(x, y)] * f_img_mag_weight;

  return cost;
}

// Returns Viterbi solution of optimal path
double *viterbi::find_path(void)
{
  start_col = f_bounds[0];
  end_col = f_bounds[1];
  num_col_vis = end_col - start_col + 1;

  int *path = new int[f_row * num_col_vis];
  double path_prob[f_row], path_prob_next[f_row], index[f_row];

  for (int k = 0; k < f_col; ++k)
    f_result[k] = 0;

  for (int k = 0; k < f_row * num_col_vis; ++k)
  {
    path[k] = 0;
  }

  for (int k = 0; k < f_row; ++k)
  {
    path_prob[k] = 0;
    path_prob_next[k] = 0;
    index[k] = 0;
  }

  viterbi_right(path, path_prob, path_prob_next, index);

  int encode;
  int viterbi_index = calculate_best(path_prob);
  int idx = end_col;
  f_result[end_col] = (f_mask[end_col] == 1 || std::isinf(f_mask[end_col])) ? viterbi_index + 1 : f_sgt[end_col] + 1;

  // Set result vector
  for (int k = start_col + 1; k <= end_col; ++k)
  {
    encode = vic_encode(viterbi_index, num_col_vis + start_col - k);
    viterbi_index = path[encode];
    f_result[idx - 1] = viterbi_index + 1;  // Account for matlab 1-indexing
    --idx;
    if (encode < 0 || idx < 1)
    {
      break;
    }
  }

  delete[] path;
  return f_result;
}

// Select path with lowest overall cost
int viterbi::calculate_best(double *path_prob)
{
  double min = LARGE;
  int viterbi_index = 0;
  for (int k = 0; k < f_row; ++k)
  {
    if (path_prob[k] < min)
    {
      min = path_prob[k];
      viterbi_index = k;
    }
  }
  return viterbi_index;
}

// Perform viterbi to the right
void viterbi::viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index)
{
  int idx = 0;
  bool next = 0;

  for (int col = start_col; col <= end_col; ++col)
  {
    if (idx >= f_row * num_col_vis || col >= f_col || col < 0)
    {
      continue;
    }
    // Have to add unary cost to first column before calculating best prev index for next column
    for (int row = 0; row < f_row; ++row)
    {
      if (col > start_col)
      {
        path[idx] = index[row];
      }
      if (next)
      {
        path_prob_next[row] += unary_cost(col, row);
      }
      else
      {
        path_prob[row] += unary_cost(col, row);
      }
      ++idx;
    }

    if (col >= end_col)
    {
      // Allow addition of unary cost to final column but do not
      // add binary cost for next column since there are no more columns

      // Overwrite path_prob with path_prob_next if that is more recent
      // path_prob is always searched for the best index after the viterbi_right call
      if (next) 
      {
        for (int i = 0; i < f_row; i++) 
        {
          path_prob[i] = path_prob_next[i];
        }
      }

      continue;
    }

    if (!next)
    {
      dt_1d(path_prob, f_transition_weights[col], path_prob_next, index, 0, f_row-1, f_smooth_slope[col], f_max_slope);
    }
    else
    {
      dt_1d(path_prob_next, f_transition_weights[col], path_prob, index, 0, f_row-1, f_smooth_slope[col], f_max_slope);
    }
    next = !next;
  }
}

// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int max_args = 18;           // All args including optional
  int min_args = max_args - 1; // Number of required arguments
  int arg = 0;

  if (nrhs < min_args || nrhs > max_args || nlhs != 1)
  {
    mexErrMsgTxt("Usage: [labels] = viterbi(input_img, surface_gt, extra_gt, ice_mask, img_mag_weight, smooth_slope, max_slope, bounds, gt_weights, gt_cutoffs, mask_dist, costmatrix, transition_weights, surf_weight, mult_weight, mult_weight_decay, mult_weight_local_decay, [zero_bin])\n");
  }

  // Input checking
  // input image ========================================================
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: image must be type double");
  }
  if (mxGetNumberOfDimensions(prhs[arg]) != 2)
  {
    mexErrMsgTxt("usage: image must be a 2D matrix");
  }
  const int _row = mxGetM(prhs[arg]);
  const int _col = mxGetN(prhs[arg]);
  const double *_image = mxGetPr(prhs[arg]);

  // surface ground truth ===============================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: sgt must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col)
  {
    mexErrMsgTxt("usage: sgt must have numel(sgt)=size(image,2)");
  }
  const double *_surf_tr = mxGetPr(prhs[arg]);

  // extra ground truth =================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: egt must be type double");
  }
  const int _num_extra_tr = mxGetNumberOfElements(prhs[arg]);
  if (_num_extra_tr % 2 != 0)
  {
    mexErrMsgTxt("usage: egt size must be a multiple of 2");
  }
  if (mxGetNumberOfElements(prhs[arg]) > 0)
  {
    if (mxGetNumberOfDimensions(prhs[arg]) != 2)
    {
      mexErrMsgTxt("usage: egt must be a 2xNgt array");
    }
  }
  const double *t_egt = mxGetPr(prhs[arg]);

  // mask ===============================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: mask must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col)
  {
    mexErrMsgTxt("usage: mask must have numel(mask)=size(image,2)");
  }
  const double *_mask = mxGetPr(prhs[arg]);

  // img_mag_weight ==========================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: img_mag_weight must be scalar double");
  }
  const double _img_mag_weight = *(double *)mxGetPr(prhs[arg]);

  // smooth_slope =======================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: smooth_slope must be type double");
  }
  if (_col - 1 != mxGetNumberOfElements(prhs[arg]))
  {
    mexErrMsgTxt("usage: smooth_slope must have numel(smooth_slope)=size(image,2)-1");
  }
  const double *_smooth_slope = mxGetPr(prhs[arg]);

  // max_slope =======================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: max_slope must be type scalar double");
  }
  const double _max_slope = *(double *)mxGetPr(prhs[arg]);

  // bounds =============================================================
  arg++;
  ptrdiff_t _bounds[2];
  if (mxGetNumberOfElements(prhs[arg]) != 0)
  {
    if (!mxIsInt64(prhs[arg]))
    {
      mexErrMsgTxt("Usage: bounds must be type int64");
    }
    if (mxGetNumberOfElements(prhs[arg]) != 2)
    {
      mexErrMsgTxt("Usage: bounds must be a 2 element vector");
    }
    ptrdiff_t *tmp = (ptrdiff_t *)mxGetPr(prhs[arg]);
    _bounds[0] = tmp[0] - 1;  // Account for matlab 1-indexing
    _bounds[1] = tmp[1] - 1;
    if (_bounds[0] < 0)
    {
      _bounds[0] = 0;
    }
    if (_bounds[1] < 0)
    {
      _bounds[1] = _col - 1;
    }
    if (_bounds[0] > _col)
    {
      mexErrMsgTxt("Usage: bounds[0] <= size(input,2)");
    }
    if (_bounds[1] > _col)
    {
      mexErrMsgTxt("Usage: bounds[1] <= size(input,2)");
    }
    if (_bounds[1] < _bounds[0])
    {
      mexErrMsgTxt("Usage: bounds[1] must be greater than bounds[0]");
    }
  }
  else
  {
    // Default setting is to process all columns
    _bounds[0] = 0;
    _bounds[1] = _col - 1;
  }

  // gt_weights ======================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: gt_weights must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col)
  {
    mexErrMsgTxt("usage: gt_weights must have numel(gt_weights)=size(image,2)");
  }
  const double *_gt_weights = mxGetPr(prhs[arg]);

  // gt_cutoffs ======================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: gt_cutoffs must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col)
  {
    mexErrMsgTxt("usage: gt_cutoffs must have numel(gt_weights)=size(image,2)");
  }
  const double *_gt_cutoffs = mxGetPr(prhs[arg]);

  // mask_dist ==========================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: mask_dist must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col)
  {
    mexErrMsgTxt("usage: mask_dist must have numel(mask_dist)=size(image,2)");
  }
  const double *_mask_dist = mxGetPr(prhs[arg]);

  // costmatrix =========================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: costmatrix must be type double");
  }
  const double *_costmatrix = mxGetPr(prhs[arg]);
  const int _costmatrix_X = mxGetM(prhs[arg]);
  const int _costmatrix_Y = mxGetN(prhs[arg]);

  // transition_weights ===================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgTxt("usage: transition_weights must be type double");
  }
  if (mxGetNumberOfElements(prhs[arg]) != _col - 1)
  {
    mexErrMsgTxt("usage: transition_weights have numel(transition_weights)=size(image,2)-1");
  }
  const double *_transition_weights = mxGetPr(prhs[arg]);

  // surf_weight ===================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: surf_weight must be type scalar double");
  }
  double _surf_weight = *(double *)mxGetPr(prhs[arg]);
  if (_surf_weight == -1)
  {
    _surf_weight = SURF_WEIGHT;
  }

  // mult_weight ===================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: mult_weight must be type scalar double");
  }
  double _mult_weight = *(double *)mxGetPr(prhs[arg]);
  if (_mult_weight == -1)
  {
    _mult_weight = MULT_WEIGHT;
  }

  // mult_weight_decay ===================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: mult_weight_decay must be type scalar double");
  }
  double _mult_weight_decay = *(double *)mxGetPr(prhs[arg]);
  if (_mult_weight_decay == -1)
  {
    _mult_weight_decay = MULT_WEIGHT_DECAY;
  }

  // mult_weight_local_decay ===================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgTxt("usage: mult_weight_local_decay must be type scalar double");
  }
  double _mult_weight_local_decay = *(double *)mxGetPr(prhs[arg]);
  if (_mult_weight_local_decay == -1)
  {
    _mult_weight_local_decay = MULT_WEIGHT_LOCAL_DECAY;
  }

  // zero bin ===================================================
  arg++;
  int _zero_bin = 0;
  if (nrhs >= min_args + 1)
  {
    if (!mxIsInt64(prhs[arg]))
    {
      mexErrMsgTxt("usage: zero bin must be type int64");
    }
    _zero_bin = *(int *)mxGetPr(prhs[arg]) - 1;  // Account for matlab 1-indexing
  }

  // ====================================================================

  // Assert arg == max_args
  if (arg != max_args - 1)
  {
    mexErrMsgTxt("BUG: Viterbi.cpp mex function args incorrectly ordered.");
  }

  // Initialize surface layer array
  int _sgt[_col];
  for (int k = 0; k < _col; ++k)
  {
    _sgt[k] = (int)_surf_tr[k] - 1;  // Account for matlab 1-indexing
  }

  double _egt_x[(_num_extra_tr / 2)], _egt_y[(_num_extra_tr / 2)];
  for (int p = 0; p < (_num_extra_tr / 2); ++p)
  {
    _egt_x[p] = round(t_egt[(p * 2)]) - 1;  // Account for matlab 1-indexing
    _egt_y[p] = round(t_egt[(p * 2) + 1]) - 1;
  }

  // Allocate output
  plhs[0] = mxCreateDoubleMatrix(1, _col, mxREAL);
  double *_result = mxGetPr(plhs[0]);
  viterbi obj(_row, _col, _image, _sgt, _mask, _img_mag_weight,
              _smooth_slope, _max_slope, _bounds, _num_extra_tr, _egt_x, _egt_y, 
              _gt_weights, _gt_cutoffs, _mask_dist, _costmatrix, _costmatrix_X, _costmatrix_Y, 
              _transition_weights, _surf_weight, _mult_weight, _mult_weight_decay, 
              _mult_weight_local_decay, _zero_bin, _result);
}

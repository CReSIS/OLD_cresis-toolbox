/* 
viterbi2.cpp

Layer-tracking program based on the Viterbi algorithm
 for MUSIC-processed 2D and 3D data

Authors: Victor Berger and John Paden
          Center for Remote Sensing of Ice Sheets
          2017-2018
         Adapted from original code by Mingze Xu and David Crandall

Changes to cost function (shifted exponential decay): Victor Berger and John Paden 2018
Changes to cost function (geostatistical analysis): Victor Berger and John Paden 2019
Various changes throughout (purge dead code, fix bugs, reduce inputs, remove unary cost): Reece Mathews 2020

See also: viterbi2.h

mex -v -largeArrayDims viterbi2.cpp 
*/

#include "viterbi2.h"

// Returns Viterbi solution of optimal path
double *viterbi2::find_path(void)
{

  int *path = new int[f_row * (f_col - 1)];
  double *path_prob = new double[f_row];
  double *path_prob_next = new double[f_row];
  double *index = new double[f_row];

  for (int k = 0; k < f_col; ++k)
    f_result[k] = NAN;

  for (int k = 0; k < f_row * (f_col - 1); ++k)
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
  f_result[f_col - 1] = viterbi_index + 1; // Account for matlab 1-indexing

  // Set result vector
  for (int k = f_col - 2; k >= 0; k--)
  {
    encode = vic_encode(viterbi_index, k);
    viterbi_index = path[encode];    // Follow path backwards
    f_result[k] = viterbi_index + 1; // Account for matlab 1-indexing
  }

  delete[] path;
  delete[] path_prob;
  delete[] path_prob_next;
  delete[] index;
  
  return f_result;
}

// Select path-terminal with lowest overall cost
int viterbi2::calculate_best(double *path_prob)
{
  double min = INF;
  int viterbi_index = f_upper_bounds[f_col - 1];
  // Only search within bounds of last column
  for (int k = viterbi_index; k <= f_lower_bounds[f_col - 1]; k++)
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
void viterbi2::viterbi_right(int *path, double *path_prob, double *path_prob_next, double *index)
{
  bool next = false;

  for (int col = 0; col < f_col; col++)
  {
    for (int row = f_upper_bounds[col]; row <= f_lower_bounds[col]; row++)
    {
      // Path points to best index in previous column -- not defined for first column

      // Add unary cost
      if (next)
      {
        path_prob_next[row] -= f_image[vic_encode(row, col)];
      }
      else
      {
        path_prob[row] -= f_image[vic_encode(row, col)];
      }
    }
    // No binary cost for last col -- added between columns
    if (col == f_col - 1)
    {
      break;
    }

    // Add binary cost from current to next
    if (next)
    {
      // Upper bounds are smaller than lower bounds (smallest -> largest, top -> bottom)
      dt(path_prob_next, path_prob, index, f_upper_bounds[col], f_lower_bounds[col],
         f_upper_bounds[col + 1], f_lower_bounds[col + 1], f_along_track_weight,
         f_along_track_slope[col]);
    }
    else
    {
      dt(path_prob, path_prob_next, index, f_upper_bounds[col], f_lower_bounds[col],
         f_upper_bounds[col + 1], f_lower_bounds[col + 1], f_along_track_weight,
         f_along_track_slope[col]);
    }
    next = !next;

    // Update path with updated index post-binary calculation.
    // Must use bounds from this column, cannot condense path update into next outer loop iteration
    for (int row = f_upper_bounds[col + 1]; row <= f_lower_bounds[col + 1]; row++)
    {
      path[vic_encode(row, col)] = index[row];
    }
  }

  // path_prob is checked for best index but path_prob_next is more recent,
  // update path_prob to reflect path_prob_next
  // TODO[reece]: Change pointer path instead of copying
  if (next)
  {
    for (int i = 0; i < f_row; i++)
    {
      path_prob[i] = path_prob_next[i];
    }
  }
}

/* Check bounds for invalid indices. */
void verify_bounds(double *bounds_dest, const mxArray *mx_bounds, char *bounds_name, int default_bound, int total_rows, int total_cols)
{
  char err_str[100];
  if (mxGetNumberOfElements(mx_bounds) != 0)
  {
    // Must be double for NaNs
    if (!mxIsDouble(mx_bounds))
    {
      sprintf(err_str, "Usage: %s must be type double", bounds_name);
      mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
    }
    if (mxGetNumberOfElements(mx_bounds) != total_cols)
    {
      sprintf(err_str, "Usage: numel(%s)=size(image,2)", bounds_name);
      mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
    }

    for (int i = 0; i < total_cols; i++)
    {
      bounds_dest[i] = ((double *)mxGetPr(mx_bounds))[i] - 1; // Account for matlab 1-indexing
      if (bounds_dest[i] < 0)
      {
        sprintf(err_str, "Usage: %s must be >= 1 (displaying with matlab 1-indexing: index %d is %.0f)", bounds_name, i + 1, bounds_dest[i] + 1);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (bounds_dest[i] >= total_rows)
      {
        sprintf(err_str, "Usage: %s must be within image (displaying with matlab 1-indexing: index %d is %.0f and size(image, 1) is %d)", bounds_name, i + 1, bounds_dest[i] + 1, total_rows);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (mxIsNaN(bounds_dest[i]))
      {
        bounds_dest[i] = default_bound;
      }
    }
  }
  else
  {
    // Default setting if empty array passed is to process all rows
    for (int i = 0; i < total_rows; i++)
    {
      bounds_dest[i] = default_bound;
    }
  }
}

// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int arg = 0;

  if (nrhs != 5 || nlhs != 1)
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "Usage: [labels] = viterbi2(image, along_track_slope, along_track_weight, upper_bounds, lower_bounds)\n");
  }

  // Input checking
  // input image ========================================================
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "usage: image must be type double");
  }
  if (mxGetNumberOfDimensions(prhs[arg]) != 2)
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "usage: image must be a 2D matrix");
  }
  const int _row = mxGetM(prhs[arg]);
  const int _col = mxGetN(prhs[arg]);
  const double *_image = mxGetPr(prhs[arg]);

  // along_track_slope =======================================================
  arg++;
  if (!mxIsDouble(prhs[arg]))
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "usage: along_track_slope must be type double");
  }
  if (_col - 1 != mxGetNumberOfElements(prhs[arg]))
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "usage: along_track_slope must have numel(along_track_slope)=size(image,2)-1");
  }
  const double *_along_track_slope = mxGetPr(prhs[arg]);

  // along_track_weight ==========================================================
  arg++;
  if (!mxIsDouble(prhs[arg]) || mxGetNumberOfElements(prhs[arg]) != 1)
  {
    mexErrMsgIdAndTxt("MATLAB:inputError", "usage: along_track_weight must be scalar double");
  }
  const double _along_track_weight = *(double *)mxGetPr(prhs[arg]);

  // Upper bounds =============================================================
  arg++;
  double _upper_bounds[_col];
  verify_bounds(_upper_bounds, prhs[arg], "upper_bounds", 0, _row, _col);

  // Lower bounds =============================================================
  arg++;
  double _lower_bounds[_col];
  verify_bounds(_lower_bounds, prhs[arg], "lower_bounds", _row - 1, _row, _col);

  // Verify lower_bounds below upper_bounds
  char err_str[100];
  for (int i = 0; i < _col; i++)
  {
    if (_lower_bounds[i] < _upper_bounds[i])
    {
      sprintf(err_str, "Usage: lower_bounds must be lower (i.e. >=) than upper_bounds (displaying with matlab 1-indexing: index %d of lower_bounds is %.0f and corresponding upper_bound is %.0f)", i + 1, _lower_bounds[i] + 1, _upper_bounds[i] + 1);
      mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
    }
  }

  // Allocate output
  plhs[0] = mxCreateDoubleMatrix(1, _col, mxREAL);
  double *_result = mxGetPr(plhs[0]);
  viterbi2 obj(_row, _col, _image, _along_track_slope, _along_track_weight, _upper_bounds, _lower_bounds, _result);
}
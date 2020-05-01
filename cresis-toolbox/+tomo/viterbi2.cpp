// viterbi2.cpp
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
// Various changes throughout (purge dead code, fix bugs, reduce inputs, remove unary cost): Reece Mathews 2020
//
// See also: viterbi2.h
//
// mex -v -largeArrayDims viterbi2.cpp

#include "viterbi2.h"


// Returns Viterbi solution of optimal path
double *viterbi2::find_path(void)
{

  int path[f_row * (f_col - 1)];
  double path_prob[f_row], path_prob_next[f_row], index[f_row];

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
    viterbi_index = path[encode];  // Follow path backwards
    f_result[k] = viterbi_index + 1; // Account for matlab 1-indexing
  }

  return f_result;
}

// Select path with lowest overall cost
int viterbi2::calculate_best(double *path_prob)
{
  double min = INF;
  int viterbi_index = 0;
  for (int k = 0; k < f_row; k++)
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
  bool next = 0;

  for (int col = 0; col <= f_col; ++col)
  {
    for (int row = 0; row < f_row; ++row)
    {
      // Path points to best index in previous column -- not defined for first column
      if (col > 0)
      {
        path[vic_encode(row, col-1)] = index[row];
      }

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
    if (col == f_col)
    {
      break;
    }

    // Add binary cost from current to next
    if (next)
    {
      dt_1d(path_prob_next, f_along_track_weight, path_prob, index, 0, f_row - 1, f_along_track_slope[col]);
    }
    else
    {
      dt_1d(path_prob, f_along_track_weight, path_prob_next, index, 0, f_row - 1, f_along_track_slope[col]);
    }
    next = !next;
  }

  // path_prob is checked for best index but path_prob_next is more recent,
  // update path_prob to reflect path_prob_next
  if (next)
  {
    for (int i = 0; i < f_row; i++)
    {
      path_prob[i] = path_prob_next[i];
    }
  }
}

// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int arg = 0;
  char err_str[100];

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
  if (mxGetNumberOfElements(prhs[arg]) != 0)
  {
    // Must be double for NaNs
    if (!mxIsDouble(prhs[arg]))
    {
      mexErrMsgIdAndTxt("MATLAB:inputError", "Usage: upper_bounds must be type double");
    }
    if (mxGetNumberOfElements(prhs[arg]) != _col)
    {
      mexErrMsgIdAndTxt("MATLAB:inputError", "Usage: numel(upper_bounds)=size(image,2)");
    }

    for (int i = 0; i < _col; i++)
    {
      _upper_bounds[i] = ((double *)mxGetPr(prhs[arg]))[i] - 1; // Account for matlab 1-indexing
      if (_upper_bounds[i] < 0)
      {
        sprintf(err_str, "Usage: upper_bounds must be >= 1 (displaying with matlab 1-indexing: index %d is %f)", i + 1, _upper_bounds[i] + 1);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (_upper_bounds[i] >= _row)
      {
        sprintf(err_str, "Usage: upper_bounds must be within image (displaying with matlab 1-indexing: index %d is %f and size(image, 1) is %d)", i + 1, _upper_bounds[i] + 1, _row);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (mxIsNaN(_upper_bounds[i]))
      {
        _upper_bounds[i] = 0;
      }
    }
  }
  else
  {
    // Default setting if empty array passed is to process all rows
    for (int i = 0; i < _row; i++)
    {
      _upper_bounds[i] = 0;
    }
  }

  // Lower bounds =============================================================
  arg++;
  double _lower_bounds[_col];
  if (mxGetNumberOfElements(prhs[arg]) != 0)
  {
    // Must be double for NaNs
    if (!mxIsDouble(prhs[arg]))
    {
      mexErrMsgIdAndTxt("MATLAB:inputError", "Usage: lower_bounds must be type double");
    }
    if (mxGetNumberOfElements(prhs[arg]) != _col)
    {
      mexErrMsgIdAndTxt("MATLAB:inputError", "Usage: numel(lower_bounds)=size(image,2)");
    }

    for (int i = 0; i < _col; i++)
    {
      _lower_bounds[i] = ((double *)mxGetPr(prhs[arg]))[i] - 1; // Account for matlab 1-indexing
      if (_lower_bounds[i] < 0)
      {
        sprintf(err_str, "Usage: lower_bounds must be >= 1 (displaying with matlab 1-indexing: index %d is %f)", i + 1, _lower_bounds[i] + 1);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (_lower_bounds[i] >= _row)
      {
        sprintf(err_str, "Usage: lower_bounds must be within image (displaying with matlab 1-indexing: index %d is %f and size(image, 1) is %d)", i + 1, _lower_bounds[i] + 1, _row);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
      if (mxIsNaN(_lower_bounds[i]))
      {
        _lower_bounds[i] = _row - 1;
      }
      if (_lower_bounds[i] < _upper_bounds[i])
      {
        sprintf(err_str, "Usage: lower_bounds must be lower (i.e. >=) than upper_bounds (displaying with matlab 1-indexing: index %d of lower_bounds is %f and corresponding upper_bound is %f)", i + 1, _lower_bounds[i] + 1, _upper_bounds[i] + 1);
        mexErrMsgIdAndTxt("MATLAB:inputError", err_str);
      }
    }
  }
  else
  {
    // Default setting if empty array passed is to process all rows
    for (int i = 0; i < _row; i++)
    {
      _lower_bounds[i] = _row - 1;
    }
  }

  // Allocate output
  plhs[0] = mxCreateDoubleMatrix(1, _col, mxREAL);
  double *_result = mxGetPr(plhs[0]);
  viterbi2 obj(_row, _col, _image, _along_track_slope, _along_track_weight, _upper_bounds, _lower_bounds, _result);
}

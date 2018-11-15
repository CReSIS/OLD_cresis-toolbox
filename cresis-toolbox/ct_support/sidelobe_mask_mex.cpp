// sidelobe_mask_mex.cpp
//
// C implementation of sidelobe_mask
//
// mex -v -largeArrayDims sidelobe_mask_mex.cpp
//
// Authors: John Paden

#include "matrix.h"
#include "mex.h"

// function mask = sidelobe_mask(data, sl_rows, sl)
// % mask = sidelobe_mask(data, sl_rows, sl)
// %
// % Author: John Paden
//
// mask = false(size(data));
// num_rows = size(data,1);
// num_cols = size(data,2);
//
// for col = 1:num_cols
//   for row = 1:num_rows
//     for sl_rows_idx = 1:length(sl_rows)
//       cur_row = row+sl_rows(sl_rows_idx);
//       if cur_row >= 1 && cur_row < num_rows && data(cur_row,col)+sl(sl_rows_idx) < data(row,col)
//         mask(cur_row,col) = true;
//       end
//     end
//   end
// end


// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  if (nrhs != 3 || nlhs != 1) {
    mexErrMsgTxt("Usage: mask = sidelobe_mask_mex(data, sl_rows, sl)");
  }
  
  // data ===============================================================
  if (!mxIsSingle(prhs[0])) {
    mexErrMsgTxt("usage: data must be type single");
  }
  if (mxGetNumberOfDimensions(prhs[0]) != 2) {
    mexErrMsgTxt("usage: data must be a 2D matrix [rows, columns]");
  }
  const size_t *dim_data = mxGetDimensions(prhs[0]);
  float *data = (float*)mxGetData(prhs[0]);
  // dim_data[0]: rows of data
  // dim_data[1]: cols of data
  
  // mask ===============================================================
  plhs[0] = mxCreateLogicalMatrix(dim_data[0], dim_data[1]);
  mxLogical *mask = (mxLogical*)mxGetData(plhs[0]);
  
  // sl_rows ============================================================
  if (!mxIsInt32(prhs[1])) {
    mexErrMsgTxt("usage: sl_rows must be type int32");
  }
  size_t numel_sl_rows = mxGetNumberOfElements(prhs[1]);
  int *sl_rows = (int*)mxGetData(prhs[1]);
  
  // sl =================================================================
  if (!mxIsSingle(prhs[2])) {
    mexErrMsgTxt("usage: sl must be type single");
  }
  if (mxGetNumberOfElements(prhs[2]) != sl_rows[numel_sl_rows-1]-sl_rows[0]+1) {
    mexErrMsgTxt("usage: sl must have numel equal to sl_rows(end)-sl_rows(1)+1.");
  }
  float *sl = (float*)mxGetData(prhs[2]);
  
  // Check simple case ==================================================
  if (numel_sl_rows == 0) {
    // No sidelobes so just return
    return;
  }
  
  // Sidelobe check =====================================================
  ptrdiff_t idx = 0;
  for (ptrdiff_t col = 0; col < dim_data[1]; col++) {
    for (ptrdiff_t row = 0; row < dim_data[0]; row++) {
      ptrdiff_t sl_idx;
      ptrdiff_t cur_row = row + sl_rows[0];
      // Make sure sidelobe rows do not extend before the first row
      if (cur_row < 0) {
        cur_row = 0;
        sl_idx = -cur_row;
      } else {
        sl_idx = 0;
      }
      if (cur_row < dim_data[0]) {
        // There is at least one valid bin to check
        ptrdiff_t cur_idx = idx + cur_row - row;
        ptrdiff_t stop_row = row + sl_rows[numel_sl_rows-1];
        // Make sure sidelobe rows do not extend beyond the last row
        if (stop_row >= dim_data[0]) {
          stop_row = dim_data[0]-1;
        }
        
        for (; cur_row <= stop_row; cur_row++) {
          if (data[cur_idx] + sl[sl_idx] < data[idx]) {
            // data[cur_idx] is low enough that it may be a sidelobe of data[idx]
            mask[cur_idx] = true;
          }
          cur_idx++;
          sl_idx++;
        }
        
      }
      idx++;
    }
  }
}

// function gap_idxs = data_gaps_check_cpp(master_dist, slave_dist, thresh_dist, trim_dist)
// % gap_idxs = data_gaps_check_cpp(master_dist, slave_dist, thresh_dist, trim_dist)
// %
// % Returns indices corresponding to gaps in a dataset that is being
// % interpolated. This is useful when you interpolate one dataset (e.g.
// % ATM data) onto another master time reference (e.g. radar data)
// % and there might be large gaps in the ATM dataset that should be
// % treated differently. Process:
// %   1. Interpolate data onto master axis using interp1 with extrapolation
// %   2. Set all the "gap_idxs" from this function to NaN
// %
// % master_dist = this is the master distance axis (e.g. radar along track)
// % slave_dist = this is the slave distance axis (e.g. ATM along track)
// % thresh_dist = this is the threshold distance that will be considered
// %   a gap in the data (e.g. 50 m)
// % trim_dist = if a gap is found, setting this trim input argument will
// %   cause all indices where extrapolation occurs to be considered part
// %   of the gap past this value in distance (e.g. trim == 0 means no
// %   extrapolation into a gap will be allowed, 20 m would allow extrapolation
// %   up to 20 m from the closest slave point). trim is assumed to be less
// %   than threshold.
// %
// % See example at bottom of data_gaps_check.m
// %
// % Author: John Paden
//
// mex -v -largeArrayDims data_gaps_check_mex.cpp

#include <math.h>
#include <mex.h>

#ifndef isnan
inline bool isnan(double x) {
    return x != x;
}
#endif

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  
  double *master, *slave, thresh_dist, trim_dist;
  master = mxGetPr(prhs[0]);
  slave = mxGetPr(prhs[1]);
  thresh_dist = *mxGetPr(prhs[2]);
  trim_dist = *mxGetPr(prhs[3]);
  
  ptrdiff_t num_master = (ptrdiff_t)mxGetNumberOfElements(prhs[0]);
  ptrdiff_t num_slave   = (ptrdiff_t)mxGetNumberOfElements(prhs[1]);

  // Allocate the output
  ptrdiff_t num_row = (ptrdiff_t)mxGetM(prhs[0]);
  ptrdiff_t num_col = (ptrdiff_t)mxGetN(prhs[0]);
  plhs[0] = mxCreateLogicalMatrix(num_row,num_col);
  mxLogical *gap_idxs = (mxLogical *)mxGetPr(plhs[0]);
  
  ptrdiff_t low_idx;
  ptrdiff_t high_idx;
  ptrdiff_t first_idx = 0;
  ptrdiff_t last_idx = num_slave-1;
  for (; first_idx < num_slave && isnan(slave[first_idx]); first_idx++);
  for (; last_idx >= 0 && isnan(slave[last_idx]); last_idx--);
  low_idx = first_idx;
  high_idx = first_idx;
  
  for (ptrdiff_t idx = 0; idx < num_master; idx++)
  {
    // Move forward to the first slave index beyond the current master index
    ptrdiff_t orig_high_idx = high_idx;
    for (; high_idx <= last_idx && (isnan(slave[high_idx]) || slave[high_idx] < master[idx]); high_idx++);
    if (high_idx != orig_high_idx)
    {
      low_idx = orig_high_idx;
    }
    
    // Check to see if this master index is in a gap or not
    if (high_idx > first_idx)
    {
      if (high_idx < last_idx)
      {
        // Handles all master points in between two slave indexes
        if (slave[high_idx] - slave[low_idx] > thresh_dist
                && abs(master[idx] - slave[low_idx]) > trim_dist
                && abs(master[idx] - slave[high_idx]) > trim_dist)
        {
          // The master point is between a gap (i.e. slave points that
          // contain this master point are separated by more than
          // thresh_dist). Also, the master is trim_dist away from either
          // slave point.
          gap_idxs[idx] = 1;
        }
      }
      else
      {
        // Handles master points after last slave index
        if (master[idx] - slave[last_idx] > trim_dist)
        {
          gap_idxs[idx] = 1;
        }
      }
    }
    else
    {
      // Handles master points before first slave index, also handles
      // empty slave vector
      if (num_slave == 0 || slave[first_idx] - master[idx] > trim_dist)
      {
        gap_idxs[idx] = 1;
      }
    }
  }
}

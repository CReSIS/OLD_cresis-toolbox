// Compile with mex -v -largeArrayDims detect.cpp
//   -v: verbose
//   -largeArrayDims: required for 64 bit
//
// detect.cpp: Detect the ice-bed layer in each MUSIC slice.
//
// By Mingze Xu, July 2016
// Added bounds: John Paden
//
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

using namespace std;

#include "mex.h"
#include "Instances.h"

class HMM {
public:
  // Data
  size_t width;
  size_t height;
  ptrdiff_t mid;
  vector<double> matrix;
  // Ground truth
  vector<size_t> sgt;
  size_t bgt;
  CoordType egt;
  vector<int> ice_mask;
  // Model
  size_t ms;
  vector<double> mu;
  vector<double> sigma;
  // Shape
  double egt_weight;
  double smooth_weight;
  double smooth_var;
  // 2-element vector containing start/stop index limits of the 
  // columns to run the algorithm on
  const ptrdiff_t *bounds;
  vector<double> smooth_slope;
  
  HMM(const double *input, const vector<size_t> &slayer, size_t blayer, const CoordType &elayer,
          const double *_ice_mask, const size_t _ms, const double *mean, const double *var, const size_t _width, const size_t _height,
          const ptrdiff_t _mid, const double _egt_weight, const double _smooth_weight, const double _smooth_var,
          const double *_smooth_slope, const ptrdiff_t *_bounds)
          : bgt(blayer), width(_width), height(_height), mid(_mid), egt_weight(_egt_weight),
                  smooth_weight(_smooth_weight), smooth_var(_smooth_var), ms(_ms), bounds(_bounds) {
            
            if (mid < 0)
              mid = width/2;

            // Init data
            matrix = vector<double>(width*height, 0.0);
            for (size_t i = 0; i < width*height; i++) {
              matrix[i] = input[i];
            }
            
            // Init surface ground truth
            sgt.assign(slayer.begin(), slayer.end());
            
            // Init extra ground truth
            egt.assign(elayer.begin(), elayer.end());
            
            // Init ice mask
            ice_mask = vector<int>(width, 0);
            for (size_t i = 0; i < width; i++) {
              ice_mask[i] = (int)_ice_mask[i];
            }
            
            // Init mu and sigma
            mu = vector<double>(ms, 0.0);
            sigma = vector<double>(ms, 0.0);
            for (size_t i = 0; i < ms; i++) {
              mu[i] = mean[i];
              sigma[i] = var[i];
            }
            
            for (size_t i = 0; i < width-1; i++) {
              smooth_slope.push_back(_smooth_slope[i]);
            }
          }
          
          // Index of data
          size_t encode(size_t x, size_t y) { return x*height + y; }
          // Unary cost
          double unary_cost(size_t x, size_t y);
          // Layer labeling
          vector<double> layer_labeling();
};

double HMM::unary_cost(size_t x, size_t y) {
  double cost = 0.0;
  size_t t = (ms-1)/2;
  
  // Bottom layer should be below the surface layer
  if (sgt[x] > t && sgt[x] < height-t && y+t+1 < sgt[x]) {
    return LARGE;
  }
  
  // Using bottom ground truth
  if (x == mid && (y+t < bgt || y+t > bgt+500)) {
    return LARGE;
  }
  
  // Using extra ground truth (uncomment these lines if using extra ground truth)
  for (size_t i = 0; i < egt.size(); i++) {
    if (x == egt[i].first) {
      cost += 2*sqr(abs((int)egt[i].second - (int)(y+t))/egt_weight);
      /*
       * if (y+t == egt[i].second) {
       * return cost;
       * } else {
       * return LARGE;
       * }
       */
    }
  }
  
  // Penalty if too close to surface layer
  if (abs((int)(y+t) - (int)sgt[x]) < 10) {
    cost += 100 - 10*abs((int)(y+t) - (int)sgt[x]);
  }
  
  // Template quadratic distance
  for (size_t i = 0; i < ms; i++) {
    cost += sqr(matrix[encode(x, y+i)] - mu[i]) / sigma[i];
  }
  
  return cost;
}

vector<double> HMM::layer_labeling() {
  size_t loop = 0;
  size_t next = loop+1;
  size_t depth = height-ms;
  size_t t = (ms-1)/2;
  vector<string> path[2];
  double path_prob[2][depth];
  double index[depth];
  
  path[0] = vector<string>(depth, "");
  path[1] = vector<string>(depth, "");
  for (size_t i = 0; i < depth; i++) {
    path_prob[0][i] = 0.0;
    path_prob[1][i] = 0.0;
    index[i] = 0.0;
  }
  
  ptrdiff_t start_col = bounds[0];
  ptrdiff_t stop_col = bounds[1];
  
  // First column is the left column ======================================
  {
    int col = start_col;
    if (ice_mask[col] == 0 && sgt[col] > t) {
      for (size_t row = 0; row < depth; row++) {
        path[loop%2][row] = itos((int)row);
        if (row+t != sgt[col]) {
          path_prob[loop%2][row] = LARGE;
        }
      }
    } else {
      for (size_t row = 0; row < depth; row++) {
        path[loop%2][row] = itos((int)row);
        path_prob[loop%2][row] = unary_cost(col, row);
      }
    }
  }
  
  // Continued columns ====================================================
  for (size_t col = start_col+1; col <= stop_col; col++) {
    double beta = norm_pdf((double)col, (double)mid, smooth_var, smooth_weight);
    
    // Distance transform
    dt_1d(path_prob[loop%2], beta, path_prob[next%2], index, 0, depth, smooth_slope[col-1]);
    
    if (ice_mask[col] == 0 && sgt[col] > t) {
      for (size_t row = 0; row < depth; row++) {
        path[next%2][row] = path[loop%2][(size_t)index[row]] + " " + itos((int)row);
        if (row+t != sgt[col]) {
          path_prob[next%2][row] += LARGE;
        }
      }
    } else {
      for (size_t row = 0; row < depth; row++) {
        path[next%2][row] = path[loop%2][(size_t)index[row]] + " " + itos((int)row);
        path_prob[next%2][row] += unary_cost(col, row);
      }
    }
    
    loop++;
    next++;
  }
  
  // Set solution
  double min_val = INFINITY;
  int flag = -1;
  for (size_t i = 0; i < depth; i++) {
    if (path_prob[loop%2][i] < min_val || flag == -1) {
      min_val = path_prob[loop%2][i];
      flag = (int)i;
    }
  }
  
  vector<double> result;
  size_t f = (size_t)flag;
  string labels = path[loop%2][f];
  if (!labels.empty()) {
    stringstream stream(labels);
    for (size_t col = 0; col < width; col++) {
      if (col >= start_col && col <= stop_col) {
        double n;
        stream >> n;
        if(!stream)
          break;
        result.push_back(n+t);
      } else {
        result.push_back(0);
      }
    }
  }
  
  return result;
}

// MATLAB FUNCTION START
// Convert vector to mex array
mxArray * getMexArray(const vector<double> &v) {
  mxArray * mx = mxCreateDoubleMatrix(1, v.size(), mxREAL);
  copy(v.begin(), v.end(), mxGetPr(mx));
  return mx;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  if (nrhs != 13) {
    mexErrMsgTxt("Usage: [labels] = detect(input_img, surface_gt, bottom_gt, extra_gt, ice_mask, mean, var, mid, egt_weight, smooth_weight, smooth_var, smooth_slope, bounds)\n");
  }
  
  size_t rows = mxGetM(prhs[0]);
  size_t cols = mxGetN(prhs[0]);
  size_t m = mxGetN(prhs[3]);
  double *input = mxGetPr(prhs[0]);
  double *surface = mxGetPr(prhs[1]);
  double *bottom = mxGetPr(prhs[2]);
  double *extra = mxGetPr(prhs[3]);
  double *mask = mxGetPr(prhs[4]);
  size_t ms = mxGetNumberOfElements(prhs[5]);
  double *mean = mxGetPr(prhs[5]);
  double *var = mxGetPr(prhs[6]);
  double *mid = mxGetPr(prhs[7]);
  double *egt_weight = mxGetPr(prhs[8]);
  double *smooth_weight = mxGetPr(prhs[9]);
  double *smooth_var = mxGetPr(prhs[10]);
  double *smooth_slope = mxGetPr(prhs[11]);
  
  // bounds ===============================================================
  ptrdiff_t bounds[2];
  if (nrhs >= 13 && mxGetNumberOfElements(prhs[12])) {
    if (!mxIsInt64(prhs[12])) {
      mexErrMsgTxt("usage: bounds must be type int64");
    }
    if (mxGetNumberOfElements(prhs[12]) != 2) {
      mexErrMsgTxt("usage: bounds must be a 2 element vector");
    }
    ptrdiff_t *tmp = (ptrdiff_t*)mxGetPr(prhs[12]);
    bounds[0] = tmp[0];
    bounds[1] = tmp[1];
    if (bounds[0] < 0)
      bounds[0] = 0;
    if (bounds[1] < 0)
      bounds[1] = cols-1;
    if (bounds[0] >= cols)
      mexErrMsgTxt("usage: bounds[0] < size(input,2)");
    if (bounds[1] >= cols)
      mexErrMsgTxt("usage: bounds[1] < size(input,2)");
    if (bounds[1] < bounds[0])
      mexErrMsgTxt("usage: bounds[1] must be greater than bounds[0]");
  } else {
    // Default setting is to process all columns
    bounds[0] = 0;
    bounds[1] = cols-1;
  }
  // mexPrintf("%lld %lld\n", bounds[0], bounds[1]);
  
  if (smooth_weight[0] < 0)
    smooth_weight[0] = SCALE;
  if (smooth_var[0] < 0)
    smooth_var[0] = SIGMA;
  
  // Convert surface coordinate to integer
  vector<size_t> slayer;
  for (size_t i = 0; i < cols; i++) {
    slayer.push_back(floor(surface[i]));
  }
  
  // Convert bottom coordinate to integer
  size_t blayer;
  if (bottom[0] > 0) {
    blayer = bottom[0];
  } else {
    blayer = slayer[(size_t)mid[0]]+50;
  }
  
  // Convert extra coordinate to integer
  CoordType elayer;
  for (size_t i = 0; i < m; i++) {
    elayer.push_back(pair<size_t, size_t>(floor(extra[i*2]), floor(extra[i*2+1])));
  }
  
  // Doing labeling ...
  HMM viterbi(input, slayer, blayer, elayer, mask, ms, mean, var, cols, rows, (ptrdiff_t)mid[0], egt_weight[0], smooth_weight[0], smooth_var[0], smooth_slope, bounds);
  vector<double> labels = viterbi.layer_labeling();
  
  plhs[0] = getMexArray(labels);
}
// END

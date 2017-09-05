// Compile with mex -v -largeArrayDims detect.cpp
//   -v: verbose
//   -largeArrayDims: required for 64 bit
//
// detect.cpp: Detect the ice-bed layer in each MUSIC slice.
//
// By Mingze Xu, July 2016
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
        size_t mid;
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
        vector<double> smooth_slope;

        HMM(const double *input, const vector<size_t> &slayer, size_t blayer, const CoordType &elayer, 
                const double *_ice_mask, const double *mean, const double *var, size_t _width, size_t _height,
                size_t _mid, double _egt_weight, double _smooth_weight, double _smooth_var,
                double *_smooth_slope, size_t _ms=11)
                : bgt(blayer), width(_width), height(_height), mid(_mid), egt_weight(_egt_weight),
                        smooth_weight(_smooth_weight), smooth_var(_smooth_var), ms(_ms) {
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
               if (y+t == egt[i].second) {
               return cost;
               } else {
               return LARGE;
               }
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

    // First column
    if (ice_mask[0] == 0 && sgt[0] > t) {
        for (size_t i = 0; i < depth; i++) {
            path[loop%2][i] = itos((int)i);
            if (i+t != sgt[0]) {
                path_prob[loop%2][i] = LARGE;
            }
        }
    } else {
        for (size_t i = 0; i < depth; i++) {
            path[loop%2][i] = itos((int)i);
            path_prob[loop%2][i] = unary_cost(0, i);
        }
    }

    // Continued columns
    for (size_t i = 1; i < width; i++) {
        double beta = norm_pdf((double)i, (double)mid, smooth_var, smooth_weight);

        // Distance transform
        dt_1d(path_prob[loop%2], beta, path_prob[next%2], index, 0, depth, smooth_slope[i-1]);

        if (ice_mask[i] == 0 && sgt[i] > t) {
            for (size_t j = 0; j < depth; j++) {
                path[next%2][j] = path[loop%2][(size_t)index[j]] + " " + itos((int)j);
                if (j+t != sgt[i]) {
                    path_prob[next%2][j] += LARGE;
                }
            }
        } else {
            for (size_t j = 0; j < depth; j++) {
                path[next%2][j] = path[loop%2][(size_t)index[j]] + " " + itos((int)j);
                path_prob[next%2][j] += unary_cost(i, j);
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
        for (size_t i = 0; i < width; i++) {
            double n;
            stream >> n;
            if(!stream)
                break;
            result.push_back(n+t);
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
    if (nrhs != 12) {
        mexPrintf("[labels] = detect(input_img, surface_gt, bottom_gt, extra_gt, ice_mask\n");
        mexPrintf("  mean, var, mid, egt_weight, smooth_weight, smooth_var, smooth_slope)\n");
        mexPrintf("\n");
        mexPrintf("HMM Viterbi detection function for ice-bed layer in MUSIC slices\n");
        mexPrintf("Compile with mex -v -largeArrayDims detect.cpp\n");
        mexPrintf("  -v: verbose   -largeArrayDims: required for 64 bit\n");
        mexPrintf("\n");
        mexPrintf("Inputs\n");
        mexPrintf(" input_img: radar echogram to be analyzed (double matrix Nt by N)\n");
        mexPrintf(" surface_gt: surface bins matrix (double matrix 1 by N)\n");
        mexPrintf(" bottom_gt: bottom bin (scalar for bottom at column 32)\n");
        mexPrintf(" extra_gt: manual bottom points (double matrix 2 by Ngt)\n");
        mexPrintf("   extra_gt(1,:): indicates the column in input_img of the ground truth point\n");
        mexPrintf("   extra_gt(2,:): indicates the row in input_img  of the ground truth point\n");
        mexPrintf(" ice_mask: ice(true) or rock(false) mask (double matrix 1 by N)\n");
        mexPrintf(" mean: mean of image peak template (double matrix 1 by Ntemplate)\n");
        mexPrintf(" var: variance of image peak template (double matrix 1 by Ntemplate)\n");
        mexPrintf(" mid: -1\n");
        mexPrintf(" egt_weight: weight attributed to manual ground truth points (e.g. 10).\n");
        mexPrintf(" smooth_weight: -1\n");
        mexPrintf(" smooth_var: -1\n");
        mexPrintf(" smooth_slope: to assume topography is flat (no slope): zeros(1, N - 1)\n");
        mexPrintf("\n");
        mexPrintf("Outputs\n");
        mexPrintf(" labels: location of bottom layer for each column (double 1 by N)\n");
        mexPrintf("\n");
        mexErrMsgTxt("ERROR WITH NUMBER OF INPUT ARGUMENTS - EXITING.");
    }
    
    double *input = mxGetPr(prhs[0]);
    double *surface = mxGetPr(prhs[1]);
    double *bottom = mxGetPr(prhs[2]);
    double *extra = mxGetPr(prhs[3]);
    double *mask = mxGetPr(prhs[4]);
    double *mean = mxGetPr(prhs[5]);
    double *var = mxGetPr(prhs[6]);
    double *mid = mxGetPr(prhs[7]);
    double *egt_weight = mxGetPr(prhs[8]);
    double *smooth_weight = mxGetPr(prhs[9]);
    double *smooth_var = mxGetPr(prhs[10]);
    double *smooth_slope = mxGetPr(prhs[11]);
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);
    size_t m = mxGetN(prhs[3]);

    if (mid[0] < 0)
        mid[0] = MID;
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
    HMM viterbi(input, slayer, blayer, elayer, mask, mean, var, cols, rows, (size_t)mid[0], egt_weight[0], smooth_weight[0], smooth_var[0], smooth_slope);
    vector<double> labels = viterbi.layer_labeling();

    plhs[0] = getMexArray(labels);
}
// END

// train_params.cpp: Train parameters using ground truth of surface layers.
// By Mingze Xu, July 2016
//
#include <iostream>
#include <cmath>
#include <vector>
#include "mex.h"

using namespace std;

// Return x*x
template <class T>
T sqr(T x) { return x*x; }

// Index of data
size_t encode(size_t x, size_t y, size_t height) {
    return x*height + y;
}

vector< vector<double> > train(double *input, size_t width, size_t height, double *sgt, size_t ms=11) {
    vector< vector<double> > result;
    vector<double> mu = vector<double>(ms, 0.0);
    vector<double> sigma = vector<double>(ms, 0.0);
    result.push_back(mu);
    result.push_back(sigma);

    size_t count = 0;
    size_t t = (ms-1)/2;

    // Mu
    for (size_t i = 0; i < width; i++) {
        if (sgt[i] > 0 && sgt[i]+t < height) {
            count++;
            for (size_t j = 0; j < ms; j++) {
                result[0][j] += input[encode(i, (int)(sgt[i]+j-t), height)];
            }
        }
    }

    for (size_t i = 0; i < ms; i++) {
        result[0][i] /= count;
    }

    // Sigma
    for (size_t i = 0; i < width; i++) {
        if (sgt[i] > 0 && sgt[i]+t < height) {
            for (size_t j = 0; j < ms; j++) {
                result[1][j] += sqr(input[encode(i, (int)(sgt[i]+j-t), height)] - result[0][j]);
            }
        }
    }

    for (size_t i = 0; i < ms; i++) {
        result[1][i] /= count;
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
    if (nrhs != 2) {
        cerr << "nrhs: " << nrhs << endl;
        mexErrMsgTxt("usage: train_params(input_image, surface_layer)");
    }

    double *input = mxGetPr(prhs[0]);
    double *surface = mxGetPr(prhs[1]);
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);

    vector< vector<double> > params;
    params = train(input, cols, rows, surface);

    plhs[0] = getMexArray(params[0]);
    plhs[1] = getMexArray(params[1]);
}
// END

// train_model.cpp: Train template model using ground truth of surface layers
// By Mingze Xu, July 2016
//
#include <iostream>
#include <cmath>
#include <vector>
#include "mex.h"

using namespace std;

// Compute square value
template <class T>
T sqr(T x) { return x*x; }

size_t encode(size_t w, size_t h, size_t height) {
    return w*height + h;
}

vector< vector<double> > train(double *input, size_t width, size_t height, const vector<size_t> &sgt, size_t ms=11) {
    vector< vector<double> > results;
    vector<double> mu = vector<double>(ms, 0.0);
    vector<double> sigma = vector<double>(ms, 0.0);
    results.push_back(mu);
    results.push_back(sigma);
    size_t t = (ms-1)/2;
    size_t count = 0;

    // Mu
    for (size_t i = 0; i < width; i++) {
        if (sgt[i] > t && sgt[i]+t < height) {
            count++;
            for (size_t j = 0; j < ms; j++) {
                results[0][j] += input[encode(i, sgt[i]-t+j, height)];
            }
        }
    }

    for (size_t i = 0; i < ms; i++) {
        results[0][i] /= count;
    }

    // Sigma
    for (size_t i = 0; i < width; i++) {
        if (sgt[i] > t && sgt[i]+t < height) {
            for (size_t j = 0; j < ms; j++) {
                results[1][j] += sqr(input[encode(i, sgt[i]-t+j, height)] - results[0][j]);
            }
        }
    }

    for (size_t i = 0; i < ms; i++) {
        results[1][i] /= count;
    }

    return results;
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
        mexErrMsgTxt("usage: train_model(input_image, surface_layer)");
    }

    double *input = mxGetPr(prhs[0]);
    double *surface = mxGetPr(prhs[1]);
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);

    // Change surface coordinate to integer
    vector<size_t> slayer;
    for (size_t i = 0; i < cols; i++) {
        slayer.push_back(floor(surface[i]));
    }

    // Get results
    vector< vector<double> > model;
    model = train(input, cols, rows, slayer);

    plhs[0] = getMexArray(model[0]);
    plhs[1] = getMexArray(model[1]);
}
// END

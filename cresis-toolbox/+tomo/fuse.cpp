// fuse.cpp: Fuse MUSIC slices from different beams into one.
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

// Compute pdf of normal distribution
double norm_pdf(double x, double mu=0.0, double s=10.0) {
    return (1.0/(s*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/s));
}

vector<double> fuse(double *right, double *center, double *left, size_t width, size_t height) {
    vector<double> result(width*height, 0.0);

    int l_mu = (int) (width/6);
    int c_mu = 3*l_mu;
    int r_mu = 5*l_mu;

    size_t cp = 0;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++, cp++) {
            double l_pdf = norm_pdf(i-l_mu);
            double c_pdf = norm_pdf(i-c_mu);
            double r_pdf = norm_pdf(i-r_mu);
            result[cp] = (l_pdf*left[cp] + c_pdf*center[cp] + r_pdf*right[cp]) / (l_pdf+c_pdf+r_pdf);
        }
    }

    return result;
}

// MATLAB FUNCTION START
mxArray * getMexArray(const vector<double> &v) {
    mxArray * mx = mxCreateDoubleMatrix(1, v.size(), mxREAL);
    copy(v.begin(), v.end(), mxGetPr(mx));
    return mx;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nrhs != 3) {
        cerr << "nrhs: " << nrhs << endl;
        mexErrMsgTxt("usage: fuse(right_image, center_image, left_image)");
    }

    double *input1 = mxGetPr(prhs[0]);
    double *input2 = mxGetPr(prhs[1]);
    double *input3 = mxGetPr(prhs[2]);
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);

    vector<double> fusion = fuse(input1, input2, input3, cols, rows);
    plhs[0] = getMexArray(fusion);
}
// END

// detect.cpp: Detect the ice-bed layer in the MUSIC slice.
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
        size_t width;
        size_t height;
        // Size of template
        size_t ms;
        // Data
        vector<double> matrix;
        // Surface ground truth
        vector<size_t> sgt;
        // Bottom ground truth
        size_t bgt;
        // Extra ground truth
        CoordType egt;
        // Ice mask
        vector<int> ice_mask;
        // Model
        vector<double> mu;
        vector<double> sigma;

        HMM(const double *input, const vector<size_t> &slayer, size_t blayer, const CoordType &elayer, const double *mask, const double *mean, 
                const double *var, size_t _width, size_t _height, size_t _ms=11) : bgt(blayer), width(_width), height(_height), ms(_ms) {
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
                ice_mask[i] = (int)mask[i];
            }

            // Init mu and sigma
            mu = vector<double>(ms, 0.0);
            sigma = vector<double>(ms, 0.0);
            for (size_t i = 0; i < ms; i++) {
                mu[i] = mean[i];
                sigma[i] = var[i];
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
    if (x == MID && (y+t < bgt || y+t > bgt+500)) {
        return LARGE;
    }

    // Using extra ground truth (uncomment these codes if use extra ground truth)
    for (size_t i = 0; i < egt.size(); i++) {
        if (x == egt[i].first) {
            if (y+t == egt[i].second) {
                return cost;
            } else {
                return LARGE;
            }
        }
    }

    // Penalty if too close to surface layer
    if (abs((int)y - (int)sgt[x]) < 20) {
        cost += 200;
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
        double beta = norm_pdf((double)i);

        // Distance transform
        dt_1d(path_prob[loop%2], beta, path_prob[next%2], index, 0, depth);

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
    if (nrhs != 7) {
        cerr << "nrhs: " << nrhs << endl;
        mexErrMsgTxt("usage: detect(input_image, surface_gt, bottom_gt, extra_gt, ice_mask, mean, var)");
    }

    double *input = mxGetPr(prhs[0]);
    double *surface = mxGetPr(prhs[1]);
    double *bottom = mxGetPr(prhs[2]);
    double *extra = mxGetPr(prhs[3]);
    double *mask = mxGetPr(prhs[4]);
    double *mean = mxGetPr(prhs[5]);
    double *var = mxGetPr(prhs[6]);
    size_t rows = mxGetM(prhs[0]);
    size_t cols = mxGetN(prhs[0]);
    size_t m = mxGetN(prhs[3]);

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
        blayer = slayer[MID]+50;
    }

    // Convert extra coordinate to integer
    CoordType elayer;
    for (size_t i = 0; i < m; i++) {
        elayer.push_back(pair<size_t, size_t>(floor(extra[i*2]), floor(extra[i*2+1])));
    }

    // Doing labeling ...
    HMM viterbi(input, slayer, blayer, elayer, mask, mean, var, cols, rows);
    vector<double> labels = viterbi.layer_labeling();

    plhs[0] = getMexArray(labels);
}
// END

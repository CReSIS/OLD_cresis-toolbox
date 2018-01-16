// extract.cpp: Extract 3D surface of ice-bed layers.
// By Mingze Xu, July 2016
//
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#include "mex.h"
#include "Instances.h"

class TRWSNode {
    public:
        double *prior;
        double *messages[4];

        TRWSNode(size_t max_disp) {
            prior = new double[max_disp];

            messages[dir_up] = new double[max_disp];
            messages[dir_down] = new double[max_disp];
            messages[dir_left] = new double[max_disp];
            messages[dir_right] = new double[max_disp];

            for (size_t d = 0; d < max_disp; d++) {
                messages[dir_up][d] = messages[dir_down][d] = messages[dir_left][d] = messages[dir_right][d] = 0.0;
            }
        }
        
        ~TRWSNode() {
            delete [] prior;
            delete [] messages[dir_up];
            delete [] messages[dir_down];
            delete [] messages[dir_left];
            delete [] messages[dir_right];
        }

        double get_msg(size_t dir, size_t disp) {
            return messages[dir][disp];
        }

        void set_msg(size_t dir, size_t disp, double val) {
            messages[dir][disp] = val;
        }
};

class TRWS {
    public:
        // Model size
        size_t ms;
        // Rows of one slice
        size_t depth;
        // Cols of one slice
        size_t height;
        // Number of slices
        size_t width;
        // Number of states
        size_t max_disp;

        // Dataset
        vector<double> dataset;
        // Matrix
        vector<TRWSNode *> matrix;
        // Surface ground truth
        LayerType sgt;
        // Bottom ground truth
        vector<size_t> bgt;
        // Extra ground truth
        PointType egt;
        // Ice mask
        vector<int> ice_mask;
        // Model
        vector<double> mu;
        vector<double> sigma;
        // Temporary messages
        vector<double> incomes[5];
        // Result
        vector<double> result;

        TRWS(const double *input, const size_t *dim, const LayerType &slayer, const vector<size_t> &blayer, 
                const PointType &elayer, const double *_ice_mask, const double *_mu, const double *_sigma, size_t _ms=11) : ms(_ms) {
            // Set dimensions
            depth = dim[0];
            height = dim[1];
            width = dim[2];
            max_disp = depth-ms;

            // Set dataset
            dataset = vector<double>(depth*height*width, 0.0);
            for (size_t i = 0; i < depth*height*width; i++) {
                dataset[i] = input[i];
            }

            // Set matrix and ice mask
            for (size_t i = 0; i < height*width; i++) {
                matrix.push_back(new TRWSNode(max_disp));
                ice_mask.push_back((int)_ice_mask[i]);
            }

            // Set ground truth
            sgt.assign(slayer.begin(), slayer.end());
            bgt.assign(blayer.begin(), blayer.end());
            egt.assign(elayer.begin(), elayer.end());

            // Set model
            for (size_t i = 0; i < ms; i++) {
                mu.push_back(_mu[i]);
                sigma.push_back(_sigma[i]);
            }

            // Set incomes
            incomes[dir_up] = vector<double>(max_disp, 0.0);
            incomes[dir_down] = vector<double>(max_disp, 0.0);
            incomes[dir_left] = vector<double>(max_disp, 0.0);
            incomes[dir_right] = vector<double>(max_disp, 0.0);
            incomes[dir_all] = vector<double>(max_disp, 0.0);

            // Set result
            result = vector<double>(height*width, 0.0);
        }

        ~TRWS() {
            for(size_t i = 0; i < matrix.size(); i++) {
                delete matrix[i];
            }
        }
        
        // The same to get element at 2D matrix (height, width)
        size_t encode(size_t h, size_t w);
        // The same to get element at 3D matrix in MATLAB (depth, height, width)
        size_t encode(size_t d, size_t h, size_t w);
        // Unary cost
        double unary_cost(size_t d, size_t h, size_t w);
        // Set prior
        void set_prior();
        // Set message
        double set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta);
        // Set result
        void set_result();
        // Extract surface
        void surface_extracting();
};

size_t TRWS::encode(size_t h, size_t w) {
    return h + w*height;
}

size_t TRWS::encode(size_t d, size_t h, size_t w) {
    return d + h*depth + w*depth*height;
}

double TRWS::unary_cost(size_t d, size_t h, size_t w) {
    size_t t = (ms-1)/2;

    // Ice mask
    if (ice_mask[encode(h, w)] == 0 && sgt[w][h] > t) {
        if (d+t == sgt[w][h]) {
            return 0.0;
        } else {
            return LARGE;
        }
    }

    // Bottom layer should be below surface layer
    if (sgt[w][h] > t && sgt[w][h]+t < depth && d+t+1 < sgt[w][h]) {
        return LARGE;
    }

    // Bottom ground truth
    if (h == MID && (d+t < bgt[w] || d+t > bgt[w]+500)) {
        return LARGE;
    }

    // Extra ground truth (uncomment these if use extra ground truth)
    /*
    for (size_t i = 0; i < egt.size(); i++) {
        if (w == get<0>(egt[i]) && h == get<1>(egt[i])) {
            if (d+t == get<2>(egt[i])) {
                return 0.0;
            } else {
                return LARGE;
            }
        }
    }
    */

    double cost = 0.0;
    if (abs((int)d - (int)sgt[w][h]) < 20) {
        cost += 200;
    }

    // Model quadratic distance
    for (size_t i = 0; i < ms; i++) {
        cost += sqr(dataset[encode(d+i,h,w)] - mu[i]) / sigma[i];
    }

    return cost;
}

void TRWS::set_prior() {
    for (size_t i = 0, cp = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++, cp++) {
            for (size_t d = 0; d < max_disp; d++) {
                matrix[cp]->prior[d] = unary_cost(d, j, i);
                for (size_t dir = 0; dir < 4; dir++) {
                    matrix[cp]->set_msg(dir, d, matrix[cp]->prior[d]);
                }
            }
        }
    }
}

double TRWS::set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta) {
    double message_in[max_disp];
    double message_out[max_disp];
    double path[max_disp];

    // First, delete message from dir_me
    for (size_t d = beg1; d < max_disp; d++) {
        message_in[d] = incomes[dir_all][d] - incomes[dir_me][d];
    }

    // Second, prepare message
    dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta);

    // Finally, normalize message
    double min_val = INFINITY;
    for (size_t d = beg2; d < max_disp; d++) {
        if (message_out[d] < min_val) {
            min_val = message_out[d];
        }
    }

    for (size_t d = beg2; d < max_disp; d++) {
        nd_me->set_msg(dir_me, d, message_out[d]-min_val);
    }

    return min_val;
}

void TRWS::set_result() {
    size_t t = (ms-1)/2;
    double temp = 0.0;
    double min_val = INFINITY;
    size_t flag = max_disp+1;

    for (size_t h = 0; h < height; h++) {
        for (size_t w = 0; w < width; w++) {
            size_t center = encode(h, w);
            min_val = INFINITY;
            flag = max_disp+1;

            // assert(sgt[w][MID] >= t);
            for (size_t d = sgt[w][MID]-t; d < max_disp; d++) {
                temp = matrix[center]->prior[d];

                if (h > 0) {
                    size_t up = encode(h-1, w);
                    temp += matrix[up]->get_msg(dir_down, d);
                    temp += abs(result[up] - (int)(d+t));
                }

                if (h+1 < height) {
                    size_t down = encode(h+1, w);
                    temp += matrix[down]->get_msg(dir_up, d);
                }

                if (w > 0) {
                    size_t left = encode(h, w-1);
                    temp += matrix[left]->get_msg(dir_right, d);
                    temp += abs(result[left] - (int)(d+t));
                }

                if (w+1 < width) {
                    size_t right = encode(h, w+1);
                    temp += matrix[right]->get_msg(dir_left, d);
                }

                if (temp < min_val || flag > max_disp) {
                    min_val = temp;
                    flag = d;
                }
            }

            result[center] = (int)(flag+t);
        }
    }
}

void TRWS::surface_extracting() {
    int loop = 0;
    int max_loop = 50;
    size_t t = (ms-1)/2;

    while (loop < max_loop) {
        mexPrintf("loop: %d\n", loop);
        mexEvalString("drawnow;");
        // Forward
        for (size_t h = 1; h < height-1; h++) {
            for (size_t w = 1; w < width-1; w++) {
                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                // assert(sgt[w][MID] >= t);
                for (size_t d = sgt[w][MID]-t; d < max_disp; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d]+incomes[dir_up][d]+incomes[dir_down][d]+incomes[dir_left][d]+incomes[dir_right][d];
                    incomes[dir_all][d] *= gamma;
                }

                double beta = 1.0;
                size_t beg1 = sgt[w][MID]-t;
                size_t beg2 = 0;
                // Right
                beta = norm_pdf(MID, 6.0);
                beg2 = sgt[w+1][MID]-t;
                set_message(matrix[center], dir_right, beg1, beg2, beta);
                // Down
                beta = norm_pdf(h, 6.0);
                beg2 = sgt[w][MID]-t;
                set_message(matrix[center], dir_down, beg1, beg2, beta);
            }
        }

        // Backward
        for (size_t h = height-2; h > 0; h--) {
            for (size_t w = width-2; w > 0; w--) {
                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                double min_val = INFINITY;
                // assert(sgt[w][MID] >= t);
                for (size_t d = sgt[w][MID]-t; d < max_disp; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d]+incomes[dir_up][d]+incomes[dir_down][d]+incomes[dir_left][d]+incomes[dir_right][d];
                    if (incomes[dir_all][d] < min_val) {
                        min_val = incomes[dir_all][d];
                    }
                }

                // Normalize message
                for (size_t d = sgt[w][MID]-t; d < max_disp; d++) {
                    incomes[dir_all][d] -= min_val;
                    incomes[dir_all][d] *= gamma;
                }

                double beta = 1.0;
                size_t beg1 = sgt[w][MID]-t;
                size_t beg2 = 0;
                // Left
                beta = norm_pdf(MID, 6.0);
                beg2 = sgt[w-1][MID]-t;
                set_message(matrix[center], dir_left, beg1, beg2, beta);
                // Up
                beta = norm_pdf(h, 6.0);
                beg2 = sgt[w][MID]-t;
                set_message(matrix[center], dir_up, beg1, beg2, beta);
            }
        }

        loop++;
    }

    set_result();
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
        mexErrMsgTxt("usage: extract(dataset, surface_gt, bottom_gt, extra_gt, ice_mask, mean, variance)");
    }

    double *input = mxGetPr(prhs[0]);
    // dim[0]: rows of one slice
    // dim[1]: cols of one slice
    // dim[2]: number of slices
    const size_t *dim = mxGetDimensions(prhs[0]);
    double *surface = mxGetPr(prhs[1]);
    double *bottom = mxGetPr(prhs[2]);
    double *extra = mxGetPr(prhs[3]);
    size_t ne = mxGetN(prhs[3]);
    double *mask = mxGetPr(prhs[4]);
    double *mean = mxGetPr(prhs[5]);
    double *var = mxGetPr(prhs[6]);

    // mexPrintf("rows of one slice (dim[0]): %d\n", dim[0]);
    // mexPrintf("cols of one slice (dim[1]): %d\n", dim[1]);
    // mexPrintf("number of slices (dim[2]): %d\n", dim[2]);

    // Convert surface coordinate to integer
    LayerType slayer;
    for (size_t i = 0, cp = 0; i < dim[2]; i++) {
        slayer.push_back(vector<size_t>());
        for (size_t j = 0; j < dim[1]; j++, cp++) {
            slayer[i].push_back(floor(surface[cp]));
        }
    }

    // Convert bottom coordinate to integer
    vector<size_t> blayer;
    for (size_t i = 0; i < dim[2]; i++) {
        if (bottom[i] > 0) {
            blayer.push_back(floor(bottom[i]));
        } else {
            blayer.push_back(floor(slayer[i][MID]+50));
        }
    }

    // Convert extra ground truth to integer
    PointType elayer;
    for (size_t i = 0; i < ne; i++) {
        elayer.push_back(tuple<size_t, size_t, size_t>(floor(extra[i*3]), floor(extra[i*3+1]), floor(extra[i*3+2])));
    }

    mexPrintf("Initing TRWS ...\n");
    mexEvalString("drawnow;");
    TRWS trws(input, dim, slayer, blayer, elayer, mask, mean, var);
    mexPrintf("Setting prior ...\n");
    mexEvalString("drawnow;");
    trws.set_prior();
    mexPrintf("Extracting surface ...\n");
    mexEvalString("drawnow;");
    trws.surface_extracting();

    plhs[0] = getMexArray(trws.result);
}
// END

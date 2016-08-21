// extract.cpp: Extract 3D surface of ice-bed layers.
// By Mingze Xu, July 2016
//
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

#include "mex.h"
#include "Instances.h"

// Dynamic smoothness
double norm_pdf(double x, double mu=MID, double s=12.0) {
    return ((20.0/(s*sqrt(2*M_PI))) * exp(-0.5*sqr((x-mu)/s)));
}

// ----------------------------------------------------------------------------------------------------
//                                              TRWS
// ----------------------------------------------------------------------------------------------------
class TRWSNode {
    public:
        size_t max_disp;
        double *prior;
        double *message[4];

        TRWSNode(size_t md) : max_disp(md) {
            prior = new double[max_disp];

            message[dir_up] = new double[max_disp];
            message[dir_down] = new double[max_disp];
            message[dir_left] = new double[max_disp];
            message[dir_right] = new double[max_disp];

            for (size_t i = 0; i < max_disp; i++) {
                message[dir_up][i] = message[dir_down][i] = message[dir_left][i] = message[dir_right][i] = 0.0;
            }
        }

        ~TRWSNode() {
            delete [] prior;
            delete [] message[dir_up];
            delete [] message[dir_down];
            delete [] message[dir_left];
            delete [] message[dir_right];
        }

        double get_msg(size_t dir, size_t disp) {
            return message[dir][disp];
        }

        void set_msg(size_t dir, size_t disp, double val) {
            message[dir][disp] = val;
        }
};

class TRWS {
    public:
        // size of template
        size_t ms;
        // rows of one slice
        size_t depth;
        // cols of one slice
        size_t height;
        // number of slices
        size_t width;

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
        vector<double> mu;
        vector<double> sigma;
        // Temporary income messages
        vector<double> incomes[5];
        // Final labels
        vector<double> results;

        TRWS(const double *input, const size_t *dim, LayerType slayer, vector<size_t> blayer, 
                PointType elayer, double *mask, double *mean, double *var, size_t _ms=11) : ms(_ms) {
            // Set dimensions
            depth = dim[0];
            height = dim[1];
            width = dim[2];

            // Set dataset
            dataset = vector<double>(depth*height*width, 0.0);
            for (size_t i = 0; i < depth*height*width; i++) {
                dataset[i] = input[i];
            }

            // Set matrix and ice mask
            for (size_t i = 0; i < height*width; i++) {
                matrix.push_back(new TRWSNode(depth-ms));
                ice_mask.push_back((int)mask[i]);
            }

            // Set ground truth
            sgt.assign(slayer.begin(), slayer.end());
            bgt.assign(blayer.begin(), blayer.end());
            egt.assign(elayer.begin(), elayer.end());

            // Init parameters Mu and Sigma
            for (size_t i = 0; i < ms; i++) {
                mu.push_back(mean[i]);
                sigma.push_back(var[i]);
            }

            // Init income messages
            incomes[dir_up] = vector<double>(depth-ms, 0.0);
            incomes[dir_down] = vector<double>(depth-ms, 0.0);
            incomes[dir_left] = vector<double>(depth-ms, 0.0);
            incomes[dir_right] = vector<double>(depth-ms, 0.0);
            incomes[dir_all] = vector<double>(depth-ms, 0.0);

            // Init results
            results = vector<double>(height*width, 0.0);
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
        void set_message(TRWSNode *nd_me, size_t dir_me, size_t h, size_t w);
        // Extract surface (MAIN)
        void surface_extracting();
};

size_t TRWS::encode(size_t h, size_t w) {
    return h + w*height;
}

size_t TRWS::encode(size_t d, size_t h, size_t w) {
    return d + h*depth + w*depth*height;
}

double TRWS::unary_cost(size_t d, size_t h, size_t w) {
    double cost = 0.0;
    size_t t = (ms-1)/2;

    // Ice mask
    if (ice_mask[encode(h, w)] == 0) {
        if (d == sgt[w][h]) {
            return cost;
        } else {
            return LARGE;
        }
    }

    // Bottom layer should be below the surface layer
    if (sgt[w][h] > 0 && sgt[w][h]+t < depth && d+t+1 < sgt[w][h]) {
        return LARGE;
    }

    // Using bottom ground truth
    if (h == MID && (d+t < bgt[w] || d+t > bgt[w]+500)) {
        return LARGE;
    }

    // Using extra ground truth
    for (size_t i = 0; i < egt.size(); i++) {
        if (w == get<0>(egt[i]) && h == get<1>(egt[i])) {
            if (d == get<2>(egt[i])) {
                return cost;
            } else {
                return LARGE;
            }
        }
    }

    // Penalty if too close to surface layer
    if (abs((int)d - (int)sgt[w][h]) < 20) {
        cost += 200;
    }

    // Template quadratic distance
    for (size_t i = 0; i < ms; i++) {
        cost += sqr(dataset[encode(d+i, h, w)] - mu[i]) / sigma[i];
    }

    return cost;
}

void TRWS::set_prior() {
    for (size_t i = 0, cp = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++, cp++) {
            for (size_t d = 0; d < depth-ms; d++) {
                matrix[cp]->prior[d] = unary_cost(d, j, i);
                // Init first loop messages according to prior
                for (size_t dir = 0; dir < 4; dir++) {
                    matrix[cp]->set_msg(dir, d, matrix[cp]->prior[d]);
                }
            }
        }
    }
}

void TRWS::set_message(TRWSNode *nd_me, size_t dir_me, size_t h, size_t w) {
    double message[depth-ms];
    double path[depth-ms];
    double temp_incomes[depth-ms];
    double beta = 1.0;
    // Different directions have different beta
    if (dir_me == dir_up || dir_me == dir_down) {
        beta = norm_pdf(h);
    } else {
        beta = norm_pdf(MID);
    }

    // First, delete message from dir_me
    for (size_t d = 0; d < depth-ms; d++) {
        message[d] = incomes[dir_all][d] - incomes[dir_me][d];
    }

    // Second, distance transform
    dt_1d(message, beta, temp_incomes, path, 0, depth-ms);

    double min_val = INFINITY;
    for (size_t d = 0; d < depth-ms; d++) {
        if (temp_incomes[d] < min_val)
            min_val = temp_incomes[d];
    }

    // Normalize messages
    for (size_t d = 0; d < depth-ms; d++) {
        nd_me->set_msg(dir_me, d, temp_incomes[d] - min_val);
    }
}

void TRWS::surface_extracting() {
    size_t loop = 0;
    size_t max_loop = 100;

    // Begin loop
    while (loop <= max_loop) {
        mexPrintf("Loop: %d\n", loop);
        mexEvalString("drawnow;");

        // Forward
        for (size_t h = 1; h < height-1; h++) {
            for (size_t w = 1; w < width-1; w++) {
                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                // For each possible state, get messages from 4 directions
                for (size_t d = 0; d < depth-ms; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d]+incomes[dir_up][d]+incomes[dir_down][d]+incomes[dir_left][d]+incomes[dir_right][d];
                    incomes[dir_all][d] *= gamma;
                }

                // Send messages to right and down
                set_message(matrix[center], dir_down, h, w);
                set_message(matrix[center], dir_right, h, w);
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

                // For each possible state, get messages from 4 directions
                for (size_t d = 0; d < depth-ms; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d]+incomes[dir_up][d]+incomes[dir_down][d]+incomes[dir_left][d]+incomes[dir_right][d];
                    if (incomes[dir_all][d] < min_val) {
                        min_val = incomes[dir_all][d];
                    }
                }

                // Normalize messages
                for (size_t d = 0; d < depth-ms; d++) {
                    incomes[dir_all][d] -= min_val;
                    incomes[dir_all][d] *= gamma;
                }

                // Send messages to left and up
                set_message(matrix[center], dir_left, h, w);
                set_message(matrix[center], dir_up, h, w);
            }
        }

        loop++;
    }

    // Compute results
    mexPrintf("Computing the solution ...\n");
    mexEvalString("drawnow;");
    double temp = 0.0;
    double min_val = INFINITY;
    int flag = -1;
    int t = (ms-1)/2;

    for (size_t h = 0; h < height; h++) {
        for (size_t w = 0; w < width; w++) {
            min_val = INFINITY;
            flag = -1;
            size_t center = encode(h, w);

            for (size_t d = sgt[w][MID]; d < depth-ms; d++) {
                temp = matrix[center]->prior[d];
                if (h > 0) {
                    size_t up = encode(h-1, w);
                    temp += matrix[up]->get_msg(dir_down, d) + abs(results[up] - (int)d);
                }

                if (h+1 < height) {
                    size_t down = encode(h+1, w);
                    temp += matrix[down]->get_msg(dir_up, d);
                }

                if (w > 0) {
                    size_t left = encode(h, w-1);
                    temp += matrix[left]->get_msg(dir_right, d) + abs(results[left] - (int)d);
                }

                if (w+1 < width) {
                    size_t right = encode(h, w+1);
                    temp += matrix[right]->get_msg(dir_left, d);
                }

                if (temp < min_val || flag == -1) {
                    min_val = temp;
                    flag = d;
                }
            }

            results[center] = flag + t;

        }
    }
}

// ----------------------------------------------------------------------------------------------------
//                                              TRWS-END
// ----------------------------------------------------------------------------------------------------

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

    plhs[0] = getMexArray(trws.results);
}
// END

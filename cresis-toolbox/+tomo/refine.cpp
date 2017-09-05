// refine.cpp: Extract 3D surface of ice-bed layers.
// By Mingze Xu, July 2016
// Correlation based mu/sigma, addition of smoothness, surface repulsion increased, input arg checks, merge with extract.cpp: John Paden
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

class Refine {
    public:
        // Model size
        size_t ms;
        // Rows of one slice
        size_t depth;
        // Cols of one slice
        size_t height;
        size_t midHeight;
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
        // Edge ground truth
        vector<size_t> edge[2];
        // Ice mask
        vector<int> ice_mask;
        // Model
        vector<double> mu;
        vector<double> sigma;
        // Shape
        double smooth_weight;
        double smooth_var;
        vector<double> smooth_slope;
        // Temporary messages
        vector<double> incomes[5];
        // Result
        vector<double> result;

        Refine(const double *input, const size_t *dim, const LayerType &slayer, const vector<size_t> &blayer, 
                const PointType &elayer, const double *_ice_mask, const double *_mu, const double *_sigma,
                double _smooth_weight, double _smooth_var, const double *_smooth_slope, const double *_edge, size_t _ms=11)
                : ms(_ms), smooth_weight(_smooth_weight), smooth_var(_smooth_var) {
            // Set dimensions
            depth = dim[0];
            height = dim[1];
            width = dim[2];
            midHeight = height/2+1;
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

            // Set Shape
            for (size_t i = 0; i < width-1; i++) {
                smooth_slope.push_back(_smooth_slope[i]);
            }

            // Set incomes
            incomes[dir_up] = vector<double>(max_disp, 0.0);
            incomes[dir_down] = vector<double>(max_disp, 0.0);
            incomes[dir_left] = vector<double>(max_disp, 0.0);
            incomes[dir_right] = vector<double>(max_disp, 0.0);
            incomes[dir_all] = vector<double>(max_disp, 0.0);

            // Set result
            result = vector<double>(height*width, 0.0);
            
            if (_edge != NULL) {
              for (size_t i = 0; i < height; i++) {
                edge[0].push_back(floor(_edge[i]));
                edge[1].push_back(floor(_edge[i+height]));
                result[i] = _edge[i];
                result[i+(width-1)*height] = _edge[i+height];
              }
            }
        }

        ~Refine() {
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
        double set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta, size_t w);
        // Set result
        void set_result();
        // Extract surface
        void surface_extracting();
};

size_t Refine::encode(size_t h, size_t w) {
    return h + w*height;
}

size_t Refine::encode(size_t d, size_t h, size_t w) {
    return d + h*depth + w*depth*height;
}

double Refine::unary_cost(size_t d, size_t h, size_t w) {
    double cost = 0.0;
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
    if (h == midHeight && (d+t < bgt[w]-20 || d+t > bgt[w]+20)) {
        return LARGE;
    }

    // Extra ground truth (uncomment these if use extra ground truth)
    for (size_t i = 0; i < egt.size(); i++) {
        if (w == get<0>(egt[i]) && h == get<1>(egt[i])) {
            cost += 2*pow(abs((int)get<2>(egt[i]) - (int)(d+t))/2.0,2);
        }
    }

    // Surface ground truth
    if (abs((int)(d+t) - (int)sgt[w][h]) < 40) {
        cost += 400 - 5*abs((int)(d+t) - (int)sgt[w][h]);
    }
//     if (abs((int)d - (int)sgt[w][h]) < 20) {
//         cost += 200;
//     }

    // Image magnitude template match
//     double tmp_cost = 0;
//     for (size_t i = 0; i < ms; i++) {
//         tmp_cost += sqr(dataset[encode(d+i,h,w)] - mu[i]) / sigma[i];
//     }
//     mexPrintf("%f\n", tmp_cost);
//     cost = cost+tmp_cost;
    
    double tmp_cost = 50;
    for (size_t i = 0; i < ms; i++) {
        tmp_cost -= dataset[encode(d+i,h,w)]*mu[i] / sigma[i];
    }
    cost = cost+tmp_cost;

    return cost;
}

void Refine::set_prior() {
    size_t t = (ms-1)/2;
    for (size_t i = 0, cp = 0; i < width; i++) {
        for (size_t j = 0; j < height; j++, cp++) {
            for (size_t d = 0; d < max_disp; d++) {
                if (edge[0].size() && (i == 0 || i == width-1)) {
                    // Running refine and on an edge slice
                    if (d+t == edge[i%width][j]) {
                        matrix[cp]->prior[d] = 0.0;
                    } else {
                        matrix[cp]->prior[d] = LARGE;
                    }
                } else {
                    // Running extract OR not on an edge slice
                    matrix[cp]->prior[d] = unary_cost(d, j, i);
                }
                for (size_t dir = 0; dir < 4; dir++) {
                    matrix[cp]->set_msg(dir, d, matrix[cp]->prior[d]);
                }
            }
        }
    }
}

double Refine::set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta, size_t h) {
    double message_in[max_disp];
    double message_out[max_disp];
    double path[max_disp];

    // First, delete message from dir_me
    for (size_t d = beg1; d < max_disp; d++) {
        message_in[d] = incomes[dir_all][d] - incomes[dir_me][d];
    }

    // beta: weighting and scaling of smoothness cost
    beta = norm_pdf((double)h, (double)midHeight, smooth_var, smooth_weight);
    
    // Second, prepare message
    dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta, smooth_slope[h-1]);
    

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

void Refine::set_result() {
    size_t t = (ms-1)/2;
    double temp = 0.0;
    double min_val = INFINITY;
    size_t flag = max_disp+1;

    for (size_t h = 0; h < height; h++) {
        for (size_t w = 1; w < width-1; w++) {
            size_t center = encode(h, w);
            min_val = INFINITY;
            flag = max_disp+1;

            // assert(sgt[w][midHeight] >= t);
            for (size_t d = sgt[w][midHeight]-t; d < max_disp; d++) {
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

void Refine::surface_extracting() {
    int loop = 0;
    int max_loop = 50;
    size_t t = (ms-1)/2;

    while (loop < max_loop) {
      if (loop > 0) {
        mexPrintf("\b\b\b\b\b\b\b\b\b\b\b", loop);
      }
        mexPrintf("loop: %04d\n", loop+1);
        mexEvalString("drawnow;");
        // Forward
        for (size_t h = 1; h < height-1; h++) {
            for (size_t w = 1; w < width-1; w++) {
                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                // assert(sgt[w][midHeight] >= t);
                for (size_t d = sgt[w][midHeight]-t; d < max_disp; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d]+incomes[dir_up][d]+incomes[dir_down][d]+incomes[dir_left][d]+incomes[dir_right][d];
                    incomes[dir_all][d] *= gamma;
                }

                double beta = 1.0;
                size_t beg1 = sgt[w][midHeight]-t;
                size_t beg2 = 0;
                // Right
                beta = norm_pdf(midHeight, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[w+1][midHeight]-t;
                set_message(matrix[center], dir_right, beg1, beg2, beta, h);
                // Down
                beta = norm_pdf(h, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[w][midHeight]-t;
                set_message(matrix[center], dir_down, beg1, beg2, beta, h);
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
                // assert(sgt[w][midHeight] >= t);
                for (size_t d = sgt[w][midHeight]-t; d < max_disp; d++) {
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
                for (size_t d = sgt[w][midHeight]-t; d < max_disp; d++) {
                    incomes[dir_all][d] -= min_val;
                    incomes[dir_all][d] *= gamma;
                }

                double beta = 1.0;
                size_t beg1 = sgt[w][midHeight]-t;
                size_t beg2 = 0;
                // Left
                beta = norm_pdf(midHeight, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[w-1][midHeight]-t;
                set_message(matrix[center], dir_left, beg1, beg2, beta, h);
                // Up
                beta = norm_pdf(h, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[w][midHeight]-t;
                set_message(matrix[center], dir_up, beg1, beg2, beta, h);
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
    if ((nrhs != 10 && nrhs != 11) || nlhs != 1) {
        mexErrMsgTxt("usage: surf = extract(input, surface, bottom, extra, mask, mean, variance, smooth_weight, smooth_var, smooth_slope, [edge])");
    }
    
    // input ==============================================================
    if (!mxIsDouble(prhs[0])) {
        mexErrMsgTxt("usage: input must be type double");
    }
    if (mxGetNumberOfDimensions(prhs[0]) != 3) {
        mexErrMsgTxt("usage: input must be a 3D matrix [rows=Nt, columns=Ndoa, slices=Nx]");
    }
    const size_t *dimInput = mxGetDimensions(prhs[0]);
    double *input = mxGetPr(prhs[0]);
    // dimInput[0]: rows of one slice
    // dimInput[1]: cols of one slice
    // dimInput[2]: number of slices

    // surface ============================================================
    if (!mxIsDouble(prhs[1])) {
        mexErrMsgTxt("usage: surface must be type double");
    }
    const size_t *dimSurface = mxGetDimensions(prhs[1]);
    if (dimSurface[0] != dimInput[1] || dimSurface[1] != dimInput[2]) {
        mexErrMsgTxt("usage: surface must have size(surface,1)=size(input,2) and size(surface,2)=size(input,3)");
    }
    double *surface = mxGetPr(prhs[1]);

    // bottom =============================================================
    if (!mxIsDouble(prhs[2])) {
        mexErrMsgTxt("usage: bottom must be type double");
    }
    if (mxGetNumberOfElements(prhs[2]) != dimInput[2]) {
        mexErrMsgTxt("usage: bottom must have numel equal to size(input,3)");
    }
    double *bottom = mxGetPr(prhs[2]);
                
    // extra ==============================================================
    if (!mxIsDouble(prhs[3])) {
      mexErrMsgTxt("usage: extra ground truth must be type double");
    }
    size_t dimExtra[2];
    if (mxGetNumberOfElements(prhs[3]) == 0) {
      dimExtra[1] = 0;
    }
    else
    {
      if (mxGetNumberOfDimensions(prhs[3]) != 2) {
        mexErrMsgTxt("usage: extra ground truth must be a 3xN array");
      }
      const size_t *tmp = mxGetDimensions(prhs[3]);
      dimExtra[0] = tmp[0];
      dimExtra[1] = tmp[1];
      if (dimExtra[0] != 3) {
        mexErrMsgTxt("usage: extra ground truth must be a 3xN array");
      }
    }
    double *extra = mxGetPr(prhs[3]);
    // dimInput[0]: 3 rows (column, x, and y for each ground truth
    // dimInput[1]: cols of extra ground truth

    // mask ===============================================================
    if (!mxIsDouble(prhs[4])) {
        mexErrMsgTxt("usage: mask must be type double");
    }
    const size_t *dimMask = mxGetDimensions(prhs[4]);
    if (dimMask[0] != dimInput[1] || dimMask[1] != dimInput[2]) {
        mexErrMsgTxt("usage: mask must have size(mask,1)=size(input,2) and size(mask,2)=size(input,3)");
    }
    double *mask = mxGetPr(prhs[4]);

    // mean ===============================================================
    if (!mxIsDouble(prhs[5])) {
        mexErrMsgTxt("usage: mean must be type double");
    }
    size_t numMean = mxGetNumberOfElements(prhs[5]);
    double *mean = mxGetPr(prhs[5]);

    // variance ===========================================================
    if (!mxIsDouble(prhs[6])) {
        mexErrMsgTxt("usage: variable must be type double");
    }
    if (numMean != mxGetNumberOfElements(prhs[6])) {
        mexErrMsgTxt("usage: variance must have numel=numel(variance)");
    }
    double *var = mxGetPr(prhs[6]);

    // smooth_weight ======================================================
    if (!mxIsDouble(prhs[7])) {
        mexErrMsgTxt("usage: smooth_weight must be type double");
    }
    if (1 != mxGetNumberOfElements(prhs[7])) {
        mexErrMsgTxt("usage: smooth_weight must be a scalar");
    }
    double *smooth_weight = mxGetPr(prhs[7]);
    if (smooth_weight[0] < 0)
        smooth_weight[0] = SCALE;

    // smooth_var =========================================================
    if (!mxIsDouble(prhs[8])) {
        mexErrMsgTxt("usage: smooth_var must be type double");
    }
    if (1 != mxGetNumberOfElements(prhs[8])) {
        mexErrMsgTxt("usage: smooth_var must be a scalar");
    }
    double *smooth_var = mxGetPr(prhs[8]);
    if (smooth_var[0] < 0)
        smooth_var[0] = SIGMA;

    // smooth_slope =======================================================
    if (!mxIsDouble(prhs[9])) {
        mexErrMsgTxt("usage: smooth_slope must be type double");
    }
    if (dimInput[1]-1 != mxGetNumberOfElements(prhs[9])) {
        mexErrMsgTxt("usage: smooth_slope must have numel=size(input,2)-1");
    }
    double *smooth_slope = mxGetPr(prhs[9]);

    double *edge;
    if (nrhs == 11) {
      // edge ===============================================================
      if (!mxIsDouble(prhs[10])) {
        mexErrMsgTxt("usage: edge must be type double");
      }
      const size_t *dimEdge = mxGetDimensions(prhs[10]);
      if (dimEdge[0] != dimInput[1] || dimEdge[1] != 2) {
        mexErrMsgTxt("usage: edge must have size(input,2) by 2");
      }
      edge = mxGetPr(prhs[10]);
    } else {
      edge = NULL;
    }

    //mexPrintf("rows of one slice (dimInput[0]): %lld\n", dimInput[0]);
    //mexPrintf("cols of one slice (dimInput[1]): %lld\n", dimInput[1]);
    //mexPrintf("number of slices (dimInput[2]): %lld\n", dimInput[2]);

    // Convert surface coordinate to integer
    LayerType slayer;
    for (size_t i = 0, cp = 0; i < dimInput[2]; i++) {
        slayer.push_back(vector<size_t>());
        for (size_t j = 0; j < dimInput[1]; j++, cp++) {
            slayer[i].push_back(floor(surface[cp]));
        }
    }

    // Convert bottom coordinate to integer
    vector<size_t> blayer;
    size_t midHeight = dimInput[1]/2+1;
    for (size_t i = 0; i < dimInput[2]; i++) {
        if (bottom[i] > 0) {
            blayer.push_back(floor(bottom[i]));
        } else {
            blayer.push_back(floor(slayer[i][midHeight]+50));
        }
    }

    // Convert extra ground truth to integer
    PointType elayer;
    for (size_t i = 0; i < dimExtra[1]; i++) {
        elayer.push_back(tuple<size_t, size_t, size_t>(floor(extra[i*3]), floor(extra[i*3+1]), floor(extra[i*3+2])));
    }

    Refine refine(input, dimInput, slayer, blayer, elayer, mask, mean, var, smooth_weight[0], smooth_var[0], smooth_slope, edge, numMean);
    refine.set_prior();
    refine.surface_extracting();

    plhs[0] = getMexArray(refine.result);
}
// END


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
        // Number of states
        double *prior;
        double *messages[4];

        TRWSNode(size_t max_disp) {
            prior = (double*)malloc(max_disp*sizeof(double));

            messages[dir_up] = (double*)malloc(max_disp*sizeof(double));
            messages[dir_down] = (double*)malloc(max_disp*sizeof(double));
            messages[dir_left] = (double*)malloc(max_disp*sizeof(double));
            messages[dir_right] = (double*)malloc(max_disp*sizeof(double));
        }
        
        ~TRWSNode() {
            free(prior);
            free(messages[dir_up]);
            free(messages[dir_down]);
            free(messages[dir_left]);
            free(messages[dir_right]);
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
        // Rows of one slice
        size_t depth;
        // Cols of one slice
        size_t height;
        size_t mid_height;
        // Number of slices
        size_t width;
        // Number of states
        size_t max_disp;

        // Dataset
        const double *image;
        // Matrix
        vector<TRWSNode *> matrix;
        // Surface ground truth
        const double *sgt;
        // Bottom ground truth
        const double *bgt;
        // Extra ground truth
        const double *egt;
        int egt_size;
        // Edge ground truth
        const double *edge;
        // Ice mask
        const double *ice_mask;
        // Model
        const double *mu;
        const double *sigma;
        const size_t ms;
        // Shape
        double smooth_weight;
        double smooth_var;
        const double *smooth_slope;
        // Temporary messages
        vector<double> incomes[5];
        // Result
        double *result;

        Refine(const double *_image, const size_t *dim_image, const double *_sgt, const double *_bgt, 
                const double *_egt, const int _egt_size, const double *_ice_mask, const double *_mu, const double *_sigma,
                const int _ms, double _smooth_weight, double _smooth_var, const double *_smooth_slope,
                const double *_edge, double *_result)
                : image(_image), sgt(_sgt), bgt(_bgt), egt(_egt), egt_size(_egt_size), ice_mask(_ice_mask),
                        mu(_mu), sigma(_sigma), ms(_ms), smooth_weight(_smooth_weight), smooth_var(_smooth_var),
                        smooth_slope(_smooth_slope), edge(_edge), result(_result) {
            // Set dimensions
            depth = dim_image[0];
            height = dim_image[1];
            width = dim_image[2];
            
            mid_height = height/2+1;
            max_disp = depth-ms;

            // Set matrix
            for (size_t i = 0; i < height*width; i++) {
                matrix.push_back(new TRWSNode(max_disp));
            }

            // Allocate incomes
            incomes[dir_up] = vector<double>(max_disp);
            incomes[dir_down] = vector<double>(max_disp);
            incomes[dir_left] = vector<double>(max_disp);
            incomes[dir_right] = vector<double>(max_disp);
            incomes[dir_all] = vector<double>(max_disp);
            
            // Assign results for edge conditions if applied
            if (edge != NULL) {
              for (size_t h = 0; h < height; h++) {
                result[h] = edge[h];
                result[h+(width-1)*height] = edge[h+height];
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
        double set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, size_t w);
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
    size_t center = encode(h,w);
    
    // Ice mask
    if (ice_mask[center] == 0 && sgt[center] > t) {
        if (d+t == sgt[center]) {
            return 0.0;
        } else {
            return LARGE;
        }
    }

    // Bottom layer should be below surface layer
    if (d+t+1 < sgt[center] && sgt[center] > t && sgt[center]+t < depth) {
        return LARGE;
    }

    // Bottom ground truth
    if (h == mid_height && (d+t < bgt[w]-20 || d+t > bgt[w]+20)) {
        return LARGE;
    }

    // Extra ground truth
    for (size_t i = 0; i < egt_size; i++) {
        if (w == egt[3*i] && h == egt[3*i+1]) {
            cost += 10*pow(abs((int)(egt[3*i+2]) - (int)(d+t))/1.0,2);
        }
    }

    // Surface ground truth
    if (abs((int)(d+t) - (int)sgt[center]) < 45) {
        cost += 450 - 10*abs((int)(d+t) - (int)sgt[center]);
    }

    // Image magnitude correlation
    double tmp_cost = 100;
    for (size_t i = 0; i < ms; i++) {
        tmp_cost -= image[encode(d+i,h,w)]*mu[i] / sigma[i];
    }
    cost = cost+tmp_cost;

    return cost;
}

void Refine::set_prior() {
    size_t t = (ms-1)/2;
    for (size_t w = 0, cp = 0; w < width; w++) {
        for (size_t h = 0; h < height; h++, cp++) {
            for (size_t d = 0; d < max_disp; d++) {
                if (edge != NULL && (w == 0 || w == width-1)) {
                    // Running refine and on an edge slice
                  size_t edge_idx;
                    if (w == 0) {
                      edge_idx = 0;
                    } else {
                      edge_idx = 1;
                    }
                    if (d+t == edge[h+height*edge_idx]) {
                        matrix[cp]->prior[d] = 0.0;
                    } else {
                        matrix[cp]->prior[d] = LARGE;
                    }
                } else {
                    // Running extract OR not on an edge slice
                    matrix[cp]->prior[d] = unary_cost(d, h, w);
                }
                for (size_t dir = 0; dir < 4; dir++) {
                    matrix[cp]->set_msg(dir, d, matrix[cp]->prior[d]);
                }
            }
        }
    }
}

double Refine::set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, size_t h) {
    double message_in[max_disp];
    double message_out[max_disp];
    double path[max_disp];

    // First, delete message from dir_me
    for (size_t d = beg1; d < max_disp; d++) {
        message_in[d] = incomes[dir_all][d] - incomes[dir_me][d];
    }

    // beta: weighting and scaling of smoothness cost
    double beta = norm_pdf((double)h, (double)mid_height, smooth_var, smooth_weight);
    
    // Second, prepare message
    if (smooth_slope == NULL) {
      dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta, beg2-beg1);
    } else {
      dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta, beg2-beg1 + smooth_slope[h]);
    }

    // Finally, normalize message so smallest message has a cost of zero
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
    size_t best_result;

    for (size_t h = 0; h < height; h++) {
        for (size_t w = 1; w < width-1; w++) {
            size_t center = encode(h, w);
            min_val = INFINITY;
            best_result = max_disp+1;

            for (size_t d = sgt[center]-t; d < max_disp; d++) {
                temp = matrix[center]->prior[d];

                if (h > 0) {
                    size_t up = encode(h-1, w);
                    temp += matrix[up]->get_msg(dir_down, d);
                }

                if (h+1 < height) {
                    size_t down = encode(h+1, w);
                    temp += matrix[down]->get_msg(dir_up, d);
                }

                if (w > 0) {
                    size_t left = encode(h, w-1);
                    temp += matrix[left]->get_msg(dir_right, d);
                }

                if (w+1 < width) {
                    size_t right = encode(h, w+1);
                    temp += matrix[right]->get_msg(dir_left, d);
                }

                if (temp < min_val || best_result > max_disp) {
                    min_val = temp;
                    best_result = d;
                }
            }

            result[center] = (int)(best_result+t);
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
        // Forward
        for (size_t h = 1; h < height-1; h++) {
            for (size_t w = 1; w < width-1; w++) {
                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                for (size_t d = sgt[center]-t; d < max_disp; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d] + incomes[dir_up][d]
                            + incomes[dir_down][d] + incomes[dir_left][d] + incomes[dir_right][d];
                }

                double beta;
                size_t beg1 = sgt[center]-t;
                size_t beg2;
                // Right
                beg2 = sgt[encode(h,w+1)]-t;
                set_message(matrix[center], dir_right, beg1, beg2, h);
                // Down
                beg2 = sgt[encode(h+1,w)]-t;
                set_message(matrix[center], dir_down, beg1, beg2, h);
            }
        }

        // Backward
        for (size_t h = height-2; h > 0; h--) {
            for (size_t w = width-2; w > 0; w--) {
                size_t center = encode(h, w);
                size_t center_mid = encode(mid_height, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                double min_val = INFINITY;
                for (size_t d = sgt[center]-t; d < max_disp; d++) {
                    incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                    incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    incomes[dir_all][d] = matrix[center]->prior[d] + incomes[dir_up][d]
                            + incomes[dir_down][d] + incomes[dir_left][d] + incomes[dir_right][d];
                }

                double beta;
                size_t beg1 = sgt[center]-t;
                size_t beg2;
                // Left
                beta = norm_pdf(mid_height, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[encode(h,w-1)]-t;
                set_message(matrix[center], dir_left, beg1, beg2, h);
                // Up
                beta = norm_pdf(h, 6.0, smooth_weight, smooth_var);
                beg2 = sgt[encode(h-1,w)]-t;
                set_message(matrix[center], dir_up, beg1, beg2, h);
            }
        }

        loop++;
    }

    set_result();
}

// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if ((nrhs < 9 && nrhs > 11) || nlhs != 1) {
        mexErrMsgTxt("Usage: surf = extract(image, sgt, bgt, egt, mask, mean, variance, smooth_weight, smooth_var, [smooth_slope], [edge])");
    }
    
    // image ==============================================================
    if (!mxIsDouble(prhs[0])) {
        mexErrMsgTxt("usage: image must be type double");
    }
    if (mxGetNumberOfDimensions(prhs[0]) != 3) {
        mexErrMsgTxt("usage: image must be a 3D matrix [rows=Nt, columns=Ndoa, slices=Nx]");
    }
    const size_t *dim_image = mxGetDimensions(prhs[0]);
    double *image = mxGetPr(prhs[0]);
    // dim_image[0]: rows of one slice
    // dim_image[1]: cols of one slice
    // dim_image[2]: number of slices

    // sgt ================================================================
    if (!mxIsDouble(prhs[1])) {
        mexErrMsgTxt("usage: sgt must be type double");
    }
    const size_t *dim_sgt = mxGetDimensions(prhs[1]);
    if (dim_sgt[0] != dim_image[1] || dim_sgt[1] != dim_image[2]) {
        mexErrMsgTxt("usage: sgt must have size(sgt,1)=size(image,2) and size(sgt,2)=size(image,3)");
    }
    double *sgt = mxGetPr(prhs[1]);

    // bgt ================================================================
    if (!mxIsDouble(prhs[2])) {
        mexErrMsgTxt("usage: bgt must be type double");
    }
    if (mxGetNumberOfElements(prhs[2]) != dim_image[2]) {
        mexErrMsgTxt("usage: bgt must have numel equal to size(image,3)");
    }
    double *bgt = mxGetPr(prhs[2]);
                
    // egt ================================================================
    if (!mxIsDouble(prhs[3])) {
      mexErrMsgTxt("usage: egt must be type double");
    }
    size_t egt_size = 0;
    if (mxGetNumberOfElements(prhs[3]) > 0) {
      if (mxGetNumberOfDimensions(prhs[3]) != 2) {
        mexErrMsgTxt("usage: egt must be a 3xN array");
      }
      const size_t *dim_egt = mxGetDimensions(prhs[3]);
      egt_size = dim_egt[1];
      if (dim_egt[0] != 3) {
        mexErrMsgTxt("usage: egt must be a 3xN array");
      }
    }
    double *egt = mxGetPr(prhs[3]);
    // dim_image[0]: 3 rows (column, x, and y for each ground truth
    // dim_image[1]: cols of extra ground truth

    // mask ===============================================================
    if (!mxIsDouble(prhs[4])) {
        mexErrMsgTxt("usage: mask must be type double");
    }
    const size_t *dim_mask = mxGetDimensions(prhs[4]);
    if (dim_mask[0] != dim_image[1] || dim_mask[1] != dim_image[2]) {
        mexErrMsgTxt("usage: mask must have size(mask,1)=size(image,2) and size(mask,2)=size(image,3)");
    }
    double *mask = mxGetPr(prhs[4]);

    // mean ===============================================================
    if (!mxIsDouble(prhs[5])) {
        mexErrMsgTxt("usage: mean must be type double");
    }
    size_t ms = mxGetNumberOfElements(prhs[5]);
    double *mean = mxGetPr(prhs[5]);

    // variance ===========================================================
    if (!mxIsDouble(prhs[6])) {
        mexErrMsgTxt("usage: variable must be type double");
    }
    if (ms != mxGetNumberOfElements(prhs[6])) {
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
    double *smooth_slope;
    if (nrhs >= 10 && mxGetNumberOfElements(prhs[9])) {
      if (!mxIsDouble(prhs[9])) {
        mexErrMsgTxt("usage: smooth_slope must be type double");
      }
      if (dim_image[1]-1 != mxGetNumberOfElements(prhs[9])) {
        mexErrMsgTxt("usage: smooth_slope must have numel=size(image,2)-1");
      }
      smooth_slope = mxGetPr(prhs[9]);
    } else {
      smooth_slope = NULL;
    }

    // edge ===============================================================
    double *edge;
    if (nrhs >= 11 && mxGetNumberOfElements(prhs[10])) {
      if (!mxIsDouble(prhs[10])) {
        mexErrMsgTxt("usage: edge must be type double");
      }
      const size_t *dim_edge = mxGetDimensions(prhs[10]);
      if (dim_edge[0] != dim_image[1] || dim_edge[1] != 2) {
        mexErrMsgTxt("usage: edge must have size(image,2) by 2");
      }
      edge = mxGetPr(prhs[10]);
    } else {
      edge = NULL;
    }

    //mexPrintf("rows of one slice (dim_image[0]): %lld\n", dim_image[0]);
    //mexPrintf("cols of one slice (dim_image[1]): %lld\n", dim_image[1]);
    //mexPrintf("number of slices (dim_image[2]): %lld\n", dim_image[2]);

    // Convert sgt coordinate to integer
    for (size_t i = 0; i < dim_image[1]*dim_image[2]; i++) {
      sgt[i] = floor(sgt[i]);
    }

    // Convert bgt coordinate to integer
    size_t mid_height = dim_image[1]/2+1;
    for (size_t w = 0; w < dim_image[2]; w++) {
        if (bgt[w] > 0) {
            bgt[w] = floor(bgt[w]);
        } else {
          // If bgt[i] == -1, then set to default value of sgt+50
            bgt[w] = sgt[w*dim_image[1] + mid_height]+50;
        }
    }

    // Convert egt to integer
    for (size_t i = 0; i < 3*egt_size; i++) {
      egt[i] = floor(egt[i]);
    }
    
    // Allocate output
    plhs[0]  = mxCreateDoubleMatrix(dim_image[1], dim_image[2], mxREAL);
    double *result = mxGetPr(plhs[0]);

    // Run refine/extract algorithm
    Refine refine(image, dim_image, sgt, bgt, egt, egt_size, mask, mean, var, ms, smooth_weight[0], smooth_var[0], smooth_slope, edge, result);
    refine.set_prior();
    refine.surface_extracting();

}

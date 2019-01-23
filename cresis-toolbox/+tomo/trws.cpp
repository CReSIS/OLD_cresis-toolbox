// trws.cpp
//
// Extract 3D surface of ice-bed layers
//
// Authors: 
//  Mingze Xu, July 2016
//  Correlation based mu/sigma, addition of smoothness, surface repulsion increased, input arg checks,
//    merge with extract.cpp, bounds, init and edge conditions updated: John Paden 2017
//  Minor style changes to comply with new C++ standards: Victor Berger 2018
//  Changes to cost function (shifted exponential decay): Victor Berger and John Paden 2018
//  Changes to cost function (geostatistical analysis): Victor Berger and John Paden 2019
//
// See also: trws.h
//
// mex -v -largeArrayDims trws.cpp

#include <iostream>
#include <cmath>
#include <stddef.h>
#include <vector>
using namespace std;
#include "mex.h"
#include "trws.h"

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

class TRWS {
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
    // Number of loops to run
    const int max_loop;
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
    size_t t;
    // Shape
    const double *smooth_weight;
    const double smooth_var;
    const double *smooth_slope;
    // 4-element vector containing start/stop index limits of the rows and
    // columns to run the algorithm on
    const ptrdiff_t *bounds;
    // Distance to ice-margin matrix
    const double *mask_dist;
    // Distance to ice-margin cost matrix
    const double *costmatrix;
    const int costmatrix_X;
    const int costmatrix_Y;
    // Mean and variance for transition model
    const double *transition_mu;
    const double *transition_sigma;
    
    // Temporary messages
    vector<double> incomes[5];
    // Result
    double *result;

    TRWS(const double *_image, const size_t *dim_image, const double *_sgt, const double *_bgt, 
         const double *_egt, const int _egt_size, const double *_ice_mask, const double *_mu, 
         const double *_sigma, const int _ms, const double *_smooth_weight, const double _smooth_var, 
         const double *_smooth_slope, const double *_edge, const int _max_loop, const ptrdiff_t *_bounds, 
         const double *_mask_dist, const double *_costmatrix, const int _costmatrix_X, const int _costmatrix_Y,
         const double *_transition_mu, const double *_transition_sigma, double *_result)
        : image(_image), sgt(_sgt), bgt(_bgt), egt(_egt), egt_size(_egt_size), ice_mask(_ice_mask), mu(_mu),
            sigma(_sigma), ms(_ms), smooth_weight(_smooth_weight), smooth_var(_smooth_var), smooth_slope(_smooth_slope),
            edge(_edge), max_loop(_max_loop), bounds(_bounds), mask_dist(_mask_dist), costmatrix(_costmatrix), 
            costmatrix_X(_costmatrix_X), costmatrix_Y(_costmatrix_Y), transition_mu(_transition_mu), 
            transition_sigma(_transition_sigma), result(_result) {
            // Set dimensions
            depth      = dim_image[0];
            height     = dim_image[1];
            width      = dim_image[2];
            mid_height = height/2;
            max_disp   = depth-ms;
            t          = (ms-1)/2;

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
    double set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta, size_t w);
    // Set result
    void set_result();
    // Extract surface
    void surface_extracting();
};

double TRWS::unary_cost(size_t d, size_t h, size_t w) {
    size_t center = encode(h, w);
    
    // Ice mask
    if (ice_mask[center] == 0 && sgt[center] > t) {
        if(d+t == sgt[center]){
            return 0.0;
        } else {
            return LARGE;
        }
    }
    
    // Set cost to large if bottom is above surface
    if (d+t+1 < sgt[center]) {
        return LARGE;
    }
        
    // Set cost to large if far from center ground truth (if present)
    if ((bgt[w] != -1) && (h == mid_height) && (d+t < bgt[w] - 20 || d+t > bgt[w] + 20)) {
        return LARGE;
    }

    double cost = 0;
    
    // Increase cost if far from extra ground truth
    for (size_t i = 0; i < egt_size; i++) {
        if (h == egt[3*i+1] && abs(w - egt[3*i]) < 4) {
            cost += 10*sqr(abs((int)(egt[3*i+2]) - (int)(d+t))/(1.0 + 0.5*abs(w - egt[3*i])));
        }
    }
   
    if (mask_dist && costmatrix) {
        // Distance from nearest ice-mask
        // Probabilistic measurement
        int DIM = (std::isinf(mask_dist[center]) || mask_dist[center] >=  costmatrix_Y) ? costmatrix_Y - 1 : mask_dist[center];
        cost += costmatrix[costmatrix_X*DIM + d + t  + 1 - (int)sgt[center]];
    }

    // Image magnitude correlation
    double tmp_cost = 0;
    for (size_t i = 0; i < ms; i++) {
        cost -= image[encode(d+i,h,w)]*mu[i] / sigma[i];
    }
   
    return cost;
}

size_t TRWS::encode(size_t h, size_t w) {
    return h + w*height;
}

size_t TRWS::encode(size_t d, size_t h, size_t w) {
    return d + h*depth + w*depth*height;
}

void TRWS::set_prior() {
    for (size_t w = 0, cp = 0; w < width; w++) {
        for (size_t h = 0; h < height; h++, cp++) {
            for (size_t d = 0; d < max_disp; d++) {
                if (h < bounds[0] || h > bounds[1] || w < bounds[2] || w > bounds[3]) {
                    matrix[cp]->prior[d] = 0.0;
                } 
                else if (edge != NULL && w == 0 && edge[h] > 0) {
                    if (d+t == edge[h]) {
                        matrix[cp]->prior[d] = 0.0;
                    } 
                    else {
                        matrix[cp]->prior[d] = LARGE;
                    }
                } 
                else if (edge != NULL && w == width-1 && edge[h+height] > 0) {
                    if (d+t == edge[h+height]) {
                        matrix[cp]->prior[d] = 0.0;
                    } else {
                        matrix[cp]->prior[d] = LARGE;
                    }
                } 
                else {
                    // No edge boundary conditions OR not on an edge slice
                    matrix[cp]->prior[d] = unary_cost(d, h, w);
                }
                for (size_t dir = 0; dir < 4; dir++) {
                    matrix[cp]->set_msg(dir, d, matrix[cp]->prior[d]);
                }
            }
        }
    }
}

double TRWS::set_message(TRWSNode *nd_me, size_t dir_me, size_t beg1, size_t beg2, double beta, size_t h) {
    double message_in[max_disp];
    double message_out[max_disp];
    double path[max_disp];

    // First, delete message from dir_me
    for (size_t d = beg1; d < max_disp; d++) {
        message_in[d] = incomes[dir_all][d] - incomes[dir_me][d];
    }

    // Second, prepare message
    if (smooth_slope == NULL) {
        dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta, beg1-beg2);
    } 
    else {
        dt(message_in, message_out, path, beg1, max_disp-1, beg2, max_disp-1, beta, beg1-beg2 + smooth_slope[h]);
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

void TRWS::set_result() {
    double temp = 0.0;
    double min_val = INFINITY;
    size_t best_result;

    for (size_t h = bounds[0]; h <= bounds[1]; h++) {
        for (size_t w = bounds[2]; w <= bounds[3]; w++) {
            size_t center = encode(h, w);
            min_val = INFINITY;
            best_result = max_disp+1;

            for (size_t d = max(1.0,sgt[center]-t); d < max_disp; d++) {
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

                // Check to see if this d is the minimum
                if (temp < min_val || best_result > max_disp) {
                    min_val = temp;
                    best_result = d;
                }
            }

            result[center] = (int)(best_result+t);
        }
    }
}

void TRWS::surface_extracting() { 
    // Propagate h/height from the midpoint out (mid_height) and let it be
    // biased by the current loop's results. This means the starting point
    // can have a large affect on the result because propagation away from
    // the starting point is much more effective and far reaching. In our
    // case, the midpoint always has ground truth so it is a good.
    int loop = 0;
    double heights_array[height];
    for (int i=0; i<=mid_height; i++) {
        heights_array[i] = mid_height-i; // From center to left
    }
    for (int i=mid_height+1; i<height; i++) {
        heights_array[i] = i; // From center to right
    }
    double forward_incomes[depth];

    // Propagate w/width depending on whether or not edge conditions exist
    int edge_mode = 0;
    if (edge != NULL) {
        if (edge[0] > 0) {
            edge_mode = edge_mode + 1;
        }
        if (edge[height] > 0) {
            edge_mode = edge_mode + 2;
        }
    }
    // Propagates left to right
    double width_array_left[width];
    for (int i=0; i<width; i++) {
        width_array_left[i] = i;
    }
    // Propagates right to left
    double width_array_right[width];
    for (int i=0; i<width; i++) {
        width_array_right[i] = width-1-i;
    }

    while (loop < max_loop) {
        if (loop > 0) {
            mexPrintf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", loop);
        }
        mexPrintf("  loop: %2d of %2d\n", loop+1, max_loop);
        // Forward
        for (size_t h_idx = 0; h_idx < height; h_idx++) {
            size_t h = heights_array[h_idx];
            mexEvalString("drawnow;");
            for (size_t w_idx = 0; w_idx < width; w_idx++) {
                size_t w;
                if (edge_mode == 2 || edge_mode == 3 && loop%2) {
                    w = width_array_right[w_idx];
                } 
                else {
                    w = width_array_left[w_idx];
                }

                if (h<bounds[0] || h>bounds[1] || w<bounds[2] || w>bounds[3]) {
                    continue;
                }

                size_t center = encode(h, w);
                size_t up = encode(h-1, w);
                size_t down = encode(h+1, w);
                size_t left = encode(h, w-1);
                size_t right = encode(h, w+1);

                for (size_t d = max(1.0,sgt[center]-t); d < max_disp; d++) {
                    if (h > 1)
                    {
                        incomes[dir_up][d] = matrix[up]->get_msg(dir_down, d);
                    }
                    else {
                        incomes[dir_up][d] = 0;
                    }

                    if (h < height-1)
                    {
                        incomes[dir_down][d] = matrix[down]->get_msg(dir_up, d);
                    }
                    else {
                        incomes[dir_down][d] = 0;
                    }

                    // In the case where we want unbiased propagation (no edge
                    // conditions or edge_mode = 0), we use the result in
                    // forward_incomes which is the result from the previous loop.
                    // Since the loop is starting from the left, we only need to do
                    // this for the left income since the right income is already
                    // unbiased (from the last-loop).
                    if (w > 1) {
                        if (edge_mode == 0) {
                            incomes[dir_left][d] = forward_incomes[d];
                        } 
                        else {
                            incomes[dir_left][d] = matrix[left]->get_msg(dir_right, d);
                        }
                    } 
                    else {
                        incomes[dir_left][d] = 0;
                    }

                    if (w < width-1) {
                        incomes[dir_right][d] = matrix[right]->get_msg(dir_left, d);
                    }
                    else {
                        incomes[dir_right][d] = 0;
                    }

                    incomes[dir_all][d] = matrix[center]->prior[d] + incomes[dir_up][d]
                            + incomes[dir_down][d] + incomes[dir_left][d] + incomes[dir_right][d];
                }
                if (edge_mode == 0) {
                    for (size_t d = max(1.0,sgt[right]-t); d < max_disp; d++) {
                        forward_incomes[d] = matrix[center]->get_msg(dir_right, d);
                    }
                }

                double beta;
                size_t beg1 = max(1.0,sgt[center]-t);
                size_t beg2;

                if (!transition_mu || !transition_sigma) {
                    beta = norm_pdf((double)h, (double)mid_height, smooth_var, smooth_weight[0]);
                }
                else {
                    beta = norm_pdf((double)h, transition_mu[w], transition_sigma[w], smooth_weight[0]);
                }
                    
                if (w < width-1) {
                    // Right
                    beg2 = max(1.0,sgt[encode(h,w+1)]-t);
                    set_message(matrix[center], dir_right, beg1, beg2, beta, h);
                }

                if (w > 1) {
                    // Left
                    beg2 = max(1.0,sgt[encode(h,w-1)]-t);
                    set_message(matrix[center], dir_left, beg1, beg2, beta, h);
                }

                beta = norm_pdf((double)h, (double)mid_height, smooth_var, smooth_weight[1]);
                    
                if (h < height-1) {
                    // Down
                    beg2 = max(1.0,sgt[encode(h+1,w)]-t);
                    set_message(matrix[center], dir_down, beg1, beg2, beta, h);
                }

                if (h > 1) {
                    // Up
                    beg2 = max(1.0,sgt[encode(h-1,w)]-t);
                    set_message(matrix[center], dir_up, beg1, beg2, beta, h);
                }
            }
        }
        loop++;
    }
    set_result();
}

// MATLAB FUNCTION START
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if ((nrhs < 9 && nrhs > 17) || nlhs != 1) {
        mexErrMsgTxt("Usage: labels = trws(image, sgt, bgt, egt, mask, mean, variance, smooth_weight, smooth_var, [smooth_slope], [edge], [max_loop], [bounds], [mask_dist], [costmatrix], [transition_mu], [transition_sigma])");
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
    // dim_image[0]: 3 rows (column, x, and y for each ground truth)
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
    if (mxGetNumberOfElements(prhs[7]) != 2) {
        mexErrMsgTxt("usage: smooth_weight must be a 2 element vector");
    }
    double *smooth_weight = mxGetPr(prhs[7]);
    if (smooth_weight[0] < 0) {
        smooth_weight[0] = SCALE;
    }
    if (smooth_weight[1] < 0) {
        smooth_weight[1] = SCALE;
    }

    // smooth_var =========================================================
    if (!mxIsDouble(prhs[8])) {
        mexErrMsgTxt("usage: smooth_var must be type double");
    }
    if (mxGetNumberOfElements(prhs[8]) != 1) {
        mexErrMsgTxt("usage: smooth_var must be a scalar");
    }
    double *smooth_var = mxGetPr(prhs[8]);
    if (smooth_var[0] < 0) {
        smooth_var[0] = SIGMA;
    }

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
    } 
    else {
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
    } 
    else {
        edge = NULL;
    }

    // max_loop ===========================================================
    double max_loop;
    if (nrhs >= 12 && mxGetNumberOfElements(prhs[11])) {
        if (!mxIsDouble(prhs[11])) {
            mexErrMsgTxt("usage: max_loop must be type double");
        }
        if (mxGetNumberOfElements(prhs[11]) != 1) {
            mexErrMsgTxt("usage: max_loop must be a scalar");
        }
        max_loop = floor(mxGetPr(prhs[11])[0]);
    } 
    else {
        max_loop = MAX_LOOP;
    }

    // bounds =============================================================
    ptrdiff_t bounds[4];
    if (nrhs >= 13 && mxGetNumberOfElements(prhs[12])) {
        if (!mxIsInt64(prhs[12])) {
            mexErrMsgTxt("usage: bounds must be type int64");
        }
        if (mxGetNumberOfElements(prhs[12]) != 4) {
            mexErrMsgTxt("usage: bounds must be a 4 element vector");
        }
        ptrdiff_t *tmp = (ptrdiff_t*)mxGetPr(prhs[12]);
        bounds[0] = tmp[0];
        bounds[1] = tmp[1];
        bounds[2] = tmp[2];
        bounds[3] = tmp[3];
        if (bounds[0] < 0) {
            bounds[0] = 0;
        }
        if (bounds[1] < 0) {
            bounds[1] = dim_image[1]-1;
        }
        if (bounds[2] < 0) {
            bounds[2] = 0;
        }
        if (bounds[3] < 0) {
            bounds[3] = dim_image[2]-1; 
        }
        if (bounds[0] >= dim_image[1]) {
            mexErrMsgTxt("usage: bounds[0] < size(input,2)");
        }
        if (bounds[1] >= dim_image[1]) {
            mexErrMsgTxt("usage: bounds[1] < size(input,2)");
        }
        if (bounds[1] < bounds[0]) {
            mexErrMsgTxt("usage: bounds[1] must be greater than bounds[0]");
        }
        if (bounds[2] >= dim_image[2]) {
            mexErrMsgTxt("usage: bounds[2] < size(input,3)");
        }
        if (bounds[3] >= dim_image[2]) {
            mexErrMsgTxt("usage: bounds[3] < size(input,3)");
        }
        if (bounds[3] < bounds[2]) {
            mexErrMsgTxt("usage: bounds[3] must be greater than bounds[2]");
        }
    } 
    else {
        // Default settings is to process all rows and all columns
        bounds[0] = 0;
        bounds[1] = dim_image[1]-1;
        bounds[2] = 0;
        bounds[3] = dim_image[2]-1;
    }

    // mask_dist ==========================================================
    const double *mask_dist;
    if (nrhs >= 14 && mxGetNumberOfElements(prhs[13])) {
        if (!mxIsDouble(prhs[13])) {
            mexErrMsgTxt("usage: mask_dist must be type double");
        }
        const size_t *dim_mask_dist = mxGetDimensions(prhs[13]);
        if (dim_mask_dist[0] != dim_image[1] || dim_mask_dist[1] != dim_image[2]) {
            mexErrMsgTxt("usage: mask_dist must have size(mask_dist,1)=size(image,2) and size(mask_dist,2)=size(image,3)");
        }
        mask_dist = mxGetPr(prhs[13]);
    }
    else {
        mask_dist = NULL;
    }
    
    // costmatrix =========================================================
    const double *costmatrix;
    int costmatrix_X, costmatrix_Y;
    if (nrhs >= 15 && mxGetNumberOfElements(prhs[14])) {
        if (!mxIsDouble(prhs[14])) {
            mexErrMsgTxt("usage: costmatrix must be type double");
        }
        costmatrix = mxGetPr(prhs[14]);
        costmatrix_X = mxGetM(prhs[14]);
        costmatrix_Y = mxGetN(prhs[14]);
    }
    else {
        costmatrix = NULL;
        costmatrix_X = -1;
        costmatrix_Y = -1;
    }
    
    // transition_mu ======================================================
    const double *transition_mu;
    if (nrhs >= 16 && mxGetNumberOfElements(prhs[15])) {
        if (!mxIsDouble(prhs[15])) {
            mexErrMsgTxt("usage: transition_mu must be type double");
        }
        if (mxGetNumberOfElements(prhs[15]) != dim_image[1]) {
            mexErrMsgTxt("usage: transition_mu must have size(transition_mu,1)=size(image,1)");
        }   
        transition_mu = mxGetPr(prhs[15]);
    }
    else {
        transition_mu = NULL;
    }
    
    // transition_sigma ===================================================
    const double *transition_sigma;
    if (nrhs >= 17 && mxGetNumberOfElements(prhs[16])) {
        if (!mxIsDouble(prhs[16])) {
            mexErrMsgTxt("usage: transition_sigma must be type double");
        }
        if (mxGetNumberOfElements(prhs[16]) != dim_image[1]) {
            mexErrMsgTxt("usage: transition_sigma must have size(transition_sigma,1)=size(image,1)");
        }   
        transition_sigma = mxGetPr(prhs[16]);
    }
    else {
        transition_sigma = NULL;
    }
    
    // ====================================================================
    
    // Convert sgt coordinate to integer
    for (size_t i = 0; i < dim_image[1]*dim_image[2]; i++) {
        sgt[i] = floor(sgt[i]);
    }

    // Convert bgt coordinate to integer
    size_t mid_height = dim_image[1]/2;
    for (size_t w = 0; w < dim_image[2]; w++) {
        if (bgt[w] > 0) {
            bgt[w] = floor(bgt[w]);
        } 
        else {
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

    // Run TRWS algorithm
    TRWS obj(image, dim_image, sgt, bgt, egt, egt_size, mask, mean, var, ms, smooth_weight, 
              smooth_var[0], smooth_slope, edge, max_loop, bounds, mask_dist, costmatrix, 
              costmatrix_X, costmatrix_Y, transition_mu, transition_sigma, result);

    obj.set_prior();
    obj.surface_extracting();
}


#ifndef __SAR_PROC_TASK_H_INCLUDED__
#define __SAR_PROC_TASK_H_INCLUDED__

#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;

#include "mex.h"
#include "mat.h"

#define sqr(x) ((x)*(x))

const double C = 299792458.0003452;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
void sar_pixel_CMFR(double &out_real, double &out_imag, size_t pixel_idx, size_t *support_limits, double *angle_limits);
void sar_pixel_RMFR(double &out_real, double &out_imag, size_t pixel_idx, size_t *support_limits, double *angle_limits);
double interp(double * x_lims,double *y, size_t N, double xq);
double getRefractionRange(size_t,size_t,double &);
double getRange(size_t,size_t,double &);
double getSurfaceIntersectionNewton(double,double,double);
double getSurfaceIntersectionGolden(double,double,double);
size_t getAlongTrackIdx(double &);

#endif
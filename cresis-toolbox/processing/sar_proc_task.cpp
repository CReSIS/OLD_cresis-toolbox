
#include "sar_proc_task.h"

// Real and Imaginary Parts in Pulse Compressed Radar Data
double *data_real;
double *data_imag;

//
double *phase_center;
double *along_track;
// Should be passed in using c++ indices
double *output_along_track;
double *output_pos;
double *kx_support_limits;
double *pixels;
double *matched_sig_lib_real;
double *matched_sig_lib_imag;
double dt;
double dr;
double fc;
double *fcs_x;
double *k_window;
double t0;

double surf_max;

double *surf_along_track;
double *surf_poly;
double *surf_line;

double sigma_r;

double trig_arg_mult;
size_t N_k;

size_t l_ft;
size_t l_st;

size_t pixel_rows;
size_t pixel_cols;

size_t lib_entries;
size_t lib_length;
size_t mid_idx;

size_t poly_n;

double n0;
double n1;

const double tau = 0.618033988749895;

double (*range_funct)(size_t,size_t,double &);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  
  void (*sar_proc_fh)(double &, double &, size_t , size_t *, double *);
  
  if (nrhs != 2) {
    cerr << "nrhs: " << nrhs << endl;
    mexErrMsgTxt("usage: sar_proc_task(param,data).");
  }
  
  // ADD CHECKS FOR EXISTANCE OF FIELDS
  
	mexPrintf("\tIn SAR Processor...\n");
	
	mxArray * field;
	
	field = mxGetField(prhs[0],0,"phase_center");
	if(field==NULL){
		mexErrMsgTxt("No ""phase_center"" field in param struct.");
	}
	phase_center = mxGetPr(field);
  
	field = mxGetField(prhs[0],0,"along_track");
	if(field==NULL){
		mexErrMsgTxt("No ""along_track"" field in param struct.");
	}
	along_track = mxGetPr(field);
	
	field = mxGetField(prhs[0],0,"output_along_track");
	if(field==NULL){
		mexErrMsgTxt("No ""output_along_track"" field in param struct.");
	}
	output_along_track = mxGetPr(field);
  
	field = mxGetField(prhs[0],0,"output_pos");
	if(field==NULL){
		mexErrMsgTxt("No ""output_pos"" field in param struct.");
	}
	output_pos = mxGetPr(field);
	
	field = mxGetField(prhs[0],0,"pixel");
	if(field==NULL){
		mexErrMsgTxt("No ""pixel"" field in param struct.");
	}
	pixels = mxGetPr(field);
	
	field = mxGetField(prhs[0],0,"matched_sig_lib");
	if(field==NULL){
		mexErrMsgTxt("No ""matched_sig_lib"" field in param struct.");
	}  
	// size of matched filter response library
	lib_length = mxGetM(mxGetField(prhs[0],0,"matched_sig_lib"));
	lib_entries = mxGetN(mxGetField(prhs[0],0,"matched_sig_lib"));
	matched_sig_lib_real = mxGetPr(field);
	if(mxIsComplex(field)){
    sar_proc_fh = &sar_pixel_CMFR;
		matched_sig_lib_imag = mxGetPi(field);
	}else{
    sar_proc_fh = &sar_pixel_RMFR;
	}
	
	field = mxGetField(prhs[0],0,"k_window");
	if(field==NULL){
		mexErrMsgTxt("No ""k_window"" field in param struct.");
	}
	k_window = mxGetPr(field);
	
	field = mxGetField(prhs[0],0,"fcs_x");
	if(field==NULL){
		mexErrMsgTxt("No ""fcs_x"" field in param struct.");
	}
	fcs_x = mxGetPr(field);
  
	field = mxGetField(prhs[0],0,"surf_along_track");
	if(field==NULL){
		mexErrMsgTxt("No ""surf_along_track"" field in param struct.");
	}
  surf_along_track = mxGetPr(field);
  
	field = mxGetField(prhs[0],0,"surf_poly");
	if(field==NULL){
		mexErrMsgTxt("No ""surf_poly"" field in param struct.");
	}
	surf_poly = mxGetPr(field);
	poly_n = mxGetM(field)-1;
	
	if(poly_n<2){
		mexErrMsgTxt("Polynomial at least be 1st order");
	}
  
	field = mxGetField(prhs[0],0,"surf_line");
	if(field==NULL){
		mexErrMsgTxt("No ""surf_poly"" field in param struct.");
	}
	surf_line = mxGetPr(field);
	  
	field = mxGetField(prhs[0],0,"kx_support_limits");
	if(field==NULL){
		mexErrMsgTxt("No ""kx_support_limits"" field in param struct.");
	}
	kx_support_limits = mxGetPr(field);
	
	field = mxGetField(prhs[0],0,"fc");
	if(field==NULL){
		mexErrMsgTxt("No ""fc"" field in param struct.");
	}
	fc = mxGetScalar(field);
	
	field = mxGetField(prhs[0],0,"dt");
	if(field==NULL){
		mexErrMsgTxt("No ""dt"" field in param struct.");
	}
	dt = mxGetScalar(field);
	dr = C*dt/2;
	
	field = mxGetField(prhs[0],0,"t0");
	if(field==NULL){
		mexErrMsgTxt("No ""t0"" field in param struct.");
	}
	t0 = mxGetScalar(field);
	
	field = mxGetField(prhs[0],0,"n0");
	if(field==NULL){
		mexErrMsgTxt("No ""n0"" field in param struct.");
	}
	n0 = mxGetScalar(field);
	
	field = mxGetField(prhs[0],0,"n1");
	if(field==NULL){
		mexErrMsgTxt("No ""n1"" field in param struct.");
	}
	n1 = mxGetScalar(field);
	
	field = mxGetField(prhs[0],0,"refraction_flag");
	if(field==NULL){
		range_funct = &getRefractionRange;
	}else if(mxGetScalar(field)){
		range_funct = &getRefractionRange;
	}else{
		range_funct = &getRange;
	}
  
  data_real = mxGetPr(prhs[1]);
  data_imag = mxGetPi(prhs[1]);
  
  // size of data (slow time/fast time)
  l_ft = mxGetM(prhs[1]);
  l_st = mxGetN(prhs[1]);
  
  
  // size of SAR output map
  const size_t *dimensions = mxGetDimensions(mxGetField(prhs[0],0,"pixel"));
  size_t n_dims = mxGetNumberOfDimensions(mxGetField(prhs[0],0,"pixel"));
  if(n_dims<2){
    mexErrMsgTxt("pixels input should have at least 2 dimensions (3 x columns).");
  }
  
  if(*dimensions!=3){
    mexErrMsgTxt("First dimension of pixel array must have 3 elements");
  }
  
  pixel_rows = *(dimensions+1);
  size_t p_c = 1;
  if(n_dims==3){
    p_c = *(dimensions+2);
  }
  pixel_cols = p_c;
  
  if(l_st!=mxGetN(mxGetField(prhs[0],0,"phase_center"))){
    mexErrMsgTxt("phase_center must have same number of columns as data");
  }
	
	if(l_st!=mxGetN(mxGetField(prhs[0],0,"along_track"))){
		mexErrMsgTxt("along_track must have same number of columns as data");
	}
  
	if(pixel_cols!=mxGetN(mxGetField(prhs[0],0,"output_along_track"))){
		mexErrMsgTxt("pixel must have same number of columns as output_along_track");
	}
  
	if(pixel_cols!=mxGetN(mxGetField(prhs[0],0,"output_pos"))){
		mexErrMsgTxt("output_pos must have same number of columns as output_along_track");
	}
	
	if(l_st!=mxGetN(mxGetField(prhs[0],0,"surf_along_track"))){
		mexErrMsgTxt("surf_along_track must have same number of columns along_track");
	}
	
	field = mxGetField(prhs[0],0,"surf_max");
	if(field==NULL){
		surf_max = surf_along_track[0];
		for(size_t surf_idx = 1 ; surf_idx<l_st ; surf_idx++){
			if(surf_along_track[surf_idx]>surf_max){
				surf_max = surf_along_track[surf_idx];
			}
		}
	}else{
		surf_max = mxGetScalar(field);
	}
	
	if(l_st-1!=mxGetN(mxGetField(prhs[0],0,"surf_poly"))){
		mexErrMsgTxt("surf_poly must have same number of columns along_track-1");
	}
  
  if(mxGetM(mxGetField(prhs[0],0,"phase_center"))!=3){
    mexErrMsgTxt("phase_center should have 3 columns");
  }
  
  if(l_st!=mxGetN(mxGetField(prhs[0],0,"fcs_x"))){
    mexErrMsgTxt("fcs_x must have same number of columns as phase_center");
  }
  
  if(mxGetM(mxGetField(prhs[0],0,"fcs_x"))!=3){
    mexErrMsgTxt("fcs_x should have 3 columns");
  }
	
	if(along_track[0]!=0){
		mexErrMsgTxt("First element of along_track vector should be 0");
	}
	
	for(int i = 1;i<l_st;i++){
		if(along_track[i]<along_track[i-1]){
			mexErrMsgTxt("along_track vector should be monotonically increasing");
		}
	}
  
	for(int i = 1;i<pixel_cols;i++){
		if(output_along_track[i]<output_along_track[i-1]){
			mexErrMsgTxt("output_along_track vector should be monotonically increasing");
		}
	}
  
  // calculate trig argument to avoid redundant calculation
  trig_arg_mult = 2*M_PI*fc;
  
  // index of center of library entry
  mid_idx = (lib_length-1)/2;
  
  // length of kx-domain window
  N_k = mxGetM(mxGetField(prhs[0],0,"k_window"));
  
	sigma_r = C*dt/(2*lib_length);
	
  // SAR output map
  mxArray * SAR_output = mxCreateDoubleMatrix(pixel_rows, pixel_cols, mxCOMPLEX);
  double * SAR_output_real = mxGetPr(SAR_output);
  double * SAR_output_imag = mxGetPi(SAR_output);
	
  size_t last_max_support_limits[2] = {0,l_st-1};
  // vector containing slow time limit pair for image pixels
  size_t support_limits[2] = {0,l_st-1};
  // slow time limits of last evaluated pixel in column
  size_t last_lim[2] = {0,l_st-1};
  
  // loop through pixel columns
  for(size_t st_idx = 0;st_idx<pixel_cols;st_idx++){
		
		mexPrintf("\tPixel column %d of %d\n",st_idx+1,pixel_cols);
		mexEvalString("pause(.001);");
    
    double kx_limits[2] = {0,0};
    
    double d[3];
    // top pixel coordinates of column
    size_t pixel_idx = pixel_rows*st_idx;
    double * pixel_ptr = pixels+3*pixel_idx;
    double pixel_x = *pixel_ptr;
    double pixel_y = *(pixel_ptr+1);
    double pixel_z = *(pixel_ptr+2);
    
    size_t phase_center_idx = 0;
            
		// ptr to x coordinate of phase_center over pixel
		double *phase_center_coor_ptr = phase_center + last_max_support_limits[0]*3;
    // ptr to z component of flight nadir vector phase_center over pixel
    double *fcs_x_ptr = fcs_x + last_max_support_limits[0]*3;
		
		double kx = 0;
    
    double kx_last = 0;
    // search from 1st phase center in support window of top pixel of last column to last phase center
    // assumes limits of this top pixel are further along phase center path than top pixel of last column
    // IF GOES PAST L_ST?????
    for(phase_center_idx = last_max_support_limits[0];phase_center_idx<l_st;phase_center_idx++){
      // z comp. of distance vector to pixel (post-dec to y phase_center coord.)
      d[0] = (*phase_center_coor_ptr++) - pixel_x;
      // y comp. of distance vector to pixel (post-dec to x phase_center coord.)
      d[1] = (*phase_center_coor_ptr++) - pixel_y;
      // x comp. of distance vector to pixel (post-dec to z phase_center coord. of next in phase_center)
      d[2] = (*phase_center_coor_ptr++) - pixel_z;
      
      // find inner product of distance vector with slant vector
      // find distance from pixel to phase_center
      double z = 0;
      double norm_d = 0;
      for(int i = 0;i<3;i++){
        z += d[i]* (*fcs_x_ptr++);
        norm_d += d[i]*d[i];
      }
      // kx from slant vector to phase_center-pixel vector
      kx = -z/sqrt(norm_d);
      
      // if phase center is within pixel's support kx
      if(kx<kx_support_limits[1]){
        if(phase_center_idx == last_max_support_limits[0]){
          support_limits[0] = last_max_support_limits[0];
          kx_limits[1] = kx;
        }else{
          // pre-increment support limit index
          support_limits[0] = phase_center_idx-1;
          last_lim[0] = support_limits[0];
          last_max_support_limits[0] = support_limits[0];
          kx_limits[1] = kx_last;
        }
        break;
      }
      kx_last = kx;
    }
    
		// If no phase center within support limits was found for top pixel????? Shouldn't process pixel ==== support_limits[0] = l_st???
    if(phase_center_idx==l_st){
      kx_limits[1] = kx;
      support_limits[0] = l_st-1;
      last_lim[0] = support_limits[0];
      last_max_support_limits[0] = support_limits[0];
    }
		
		// Limit search for very first pixel
		if(st_idx==0){
			last_max_support_limits[1] = last_max_support_limits[0];
		}
    
    // ptr to x coordinate of phase center over pixel
    phase_center_coor_ptr = phase_center + last_max_support_limits[1]*3;
    
    // ptr to x component of flight nadir vector of phase center over pixel
    fcs_x_ptr = fcs_x + last_max_support_limits[1]*3;
    
    // search from last phase center in support window of top pixel of last column to last phase center
    // assumes limits of this top pixel are further along phase center path than top pixel of last column
    // IF GOES PAST L_ST?????
    for(phase_center_idx = last_max_support_limits[1];phase_center_idx<l_st;phase_center_idx++){
      // x comp. of vector to pixel (post-inc to y coord.)
      d[0] = (*phase_center_coor_ptr++) - pixel_x;
      // y comp. of vector to pixel (post-inc to z coord.)
      d[1] = (*phase_center_coor_ptr++) - pixel_y;
      // z comp. of vector to pixel (post-inc to x coord. of next phase center)
      d[2] = (*phase_center_coor_ptr++) - pixel_z;
      
      // find inner product of distance vector with slant vector
      // find distance from pixel to phase center
      double z = 0;
      double norm_d = 0;
      for(int i = 0;i<3;i++){
        z += d[i]*(*fcs_x_ptr++);
        norm_d += d[i]*d[i];
      }
      // kx from slant vector to phase_center-pixel vector
      kx = -z/sqrt(norm_d);
      
      // if phase center is on or outside pixel's support kx limit
      if(kx<=kx_support_limits[0]){
        // pre-increment support limit index
        support_limits[1] = phase_center_idx;
        last_lim[1] = support_limits[1];
        last_max_support_limits[1] = support_limits[1];
        kx_limits[0] = kx;
        break;
      }
    }
    
		// If no phase center was found outside the support limits
    if(phase_center_idx==l_st){
      kx_limits[0] = kx;
    }
		
    sar_proc_fh(SAR_output_real[pixel_idx],SAR_output_imag[pixel_idx],pixel_idx,support_limits,kx_limits);
    
		pixel_ptr = pixels+3*(pixel_rows*st_idx+1);
		
    // loop through all other pixels in column to find limits
    // top -> bottom
    for(size_t pixel_idx = pixel_rows*st_idx+1;pixel_idx<(pixel_rows*(st_idx+1));pixel_idx++){
      // x coordinate of pixel
      pixel_x = *pixel_ptr++;
      // y coordinate of pixel
      pixel_y = *(pixel_ptr++);
      // z coordinate of pixel
      pixel_z = *(pixel_ptr++);
      
      // ptr to x coordinate of phase center right limit of pixel below
      phase_center_coor_ptr = phase_center+last_lim[1]*3;
      // ptr to x component of flight nadir vector of phase center right limit of pixel below
      fcs_x_ptr = fcs_x + last_lim[1]*3;
      
			size_t phase_center_idx = 0;
      //
      for(phase_center_idx = last_lim[1];phase_center_idx<l_st;phase_center_idx++){
        
        // z comp. of vector to pixel (post-inc to y coord.)
        d[0] = (*phase_center_coor_ptr++) - pixel_x;
        // y comp. of vector to pixel (post-inc to x coord.)
        d[1] = (*phase_center_coor_ptr++) - pixel_y;
        // x comp. of vector to pixel (post-inc to z coord. of next phase center)
        d[2] = (*phase_center_coor_ptr++) - pixel_z;
        
        // find inner product of distance vector with slant vector
        // find distance from pixel to phase center
        double z = 0;
        double norm_d = 0;
        for(int i = 0;i<3;i++){
          z += d[i]*(*fcs_x_ptr++);
          norm_d += d[i]*d[i];
        }
        // kx from slant vector to phase_center-pixel vector
        kx = -z/sqrt(norm_d);
        
        // if phase center is on or outside pixel's support kx limit
        if(kx<=kx_support_limits[0]){
          // pre-increment support limit index
          support_limits[1] = phase_center_idx;
          last_lim[1] = phase_center_idx;
          kx_limits[0] = kx;
          break;
        }
      }
			
			if(phase_center_idx==l_st){
				support_limits[1] =  l_st-1;
				last_lim[1] = support_limits[1];
				kx_limits[0] = kx;
			}
      
      // ptr to z coordinate of left limit of pixel below
      phase_center_coor_ptr = phase_center+last_lim[0]*3+2;
      // ptr to z component of flight nadir vector of phase center left limit of pixel below
      fcs_x_ptr = fcs_x + last_lim[0]*3+2;
      
      //
      kx_last = 0;
      for(size_t phase_center_idx = last_lim[0]+1;phase_center_idx-->0;){
        
        // x comp. of vector to pixel (post-inc to y coord.)
        d[2] = (*phase_center_coor_ptr--) - pixel_z;
        // y comp. of vector to pixel (post-inc to z coord.)
        d[1] = (*phase_center_coor_ptr--) - pixel_y;
        // z comp. of vector to pixel (post-inc to x coord. of next phase center)
        d[0] = (*phase_center_coor_ptr--) - pixel_x; 
        
        // find inner product of distance vector with slant vector
        // find distance from pixel to phase center
        double z = 0;
        double norm_d = 0;
        for(int i = 2;i>=0;i--){
          z += d[i]*(*fcs_x_ptr--);
          norm_d += d[i]*d[i];
        }
        // kx from slant vector to phase_center-pixel vector
        kx = -z/sqrt(norm_d);
        
        // if phase center is on or outside pixel's support kx limit
        if(kx>=kx_support_limits[1]){
          // pre-increment support limit index
          support_limits[0] = phase_center_idx;
          last_lim[0] = phase_center_idx;
          kx_limits[1] = kx_last;
          break;
        }
        kx_last = kx;
      }
      
      if(phase_center_idx==	UINT_MAX){
				support_limits[0] = 0;
				last_lim[0] = 0;
				kx_limits[1] = kx;
      }
      
      sar_proc_fh(SAR_output_real[pixel_idx],SAR_output_imag[pixel_idx],pixel_idx,support_limits,kx_limits);
    }		
  }
  plhs[0] = SAR_output;
}


void sar_pixel_CMFR(double &out_real, double &out_imag, size_t pixel_idx, size_t *support_limits, double *kx_limits){
  
  double *pixel_ptr = pixels+3*pixel_idx;
  double pixel_x = *pixel_ptr;
  double pixel_y = *(pixel_ptr+1);
  double pixel_z = *(pixel_ptr+2);
  
  double sum_k = 0;
  double inner_prod_real = 0;
  double inner_prod_imag = 0;
  
  // ptr to x coordinate of left-most phase center for pixel
  double *phase_center_coor_ptr = phase_center+3*support_limits[0];
  // ptr to x element of nadir vector for right-most phase center for pixel
  double *fcs_x_ptr = fcs_x+3*support_limits[0];
    
	double xi_initial = NAN;
	
  // find inner product of rx image with target at pixe
  // loop from left phase center limit -> right phase center limit
  for(size_t phase_center_idx = support_limits[0];phase_center_idx<=support_limits[1];phase_center_idx++){
    
    double fcs_vect_x = *fcs_x_ptr++;
    double fcs_vect_y = *fcs_x_ptr++;
    double fcs_vect_z = *fcs_x_ptr++;
    
    //// post-increment ptr to y coord.
    double r_x = (*phase_center_coor_ptr++) - pixel_x;
    //// post-inc ptr to z coord.
    double r_y = (*phase_center_coor_ptr++) - pixel_y;
    //// post-inc ptr to x coord. of next phase center
    double r_z = (*phase_center_coor_ptr++) - pixel_z;
    //
    //// range of pixel to phase center
    //double r = r_x*r_x+r_y*r_y+r_z*r_z;
    //r = sqrt(r);
		
		double r = range_funct(pixel_idx,phase_center_idx,xi_initial);
    
    // angle = acos( dot(range_vect,nadir_vect)/abs(range_vect) )
    double kx = -(r_z*fcs_vect_z + r_y*fcs_vect_y + r_x*fcs_vect_x)/r;
    
    // weight applied by kx domain window
    double k_weight = interp(kx_limits,k_window,N_k,kx);
    
    // time delay from phase_center to pixel and back
    double t_delay = r*2/C;
		
    // row index of peak matched filter response for simulated target at pixel
		double time_idx_double = (t_delay/dt-t0/dt);
    size_t time_idx = (size_t)time_idx_double;
		
    // index of library entree associated with sub-time-step delay
    size_t lib_idx = round(((time_idx_double)-(double)(time_idx))*lib_entries);
    if(lib_idx==lib_entries){
      lib_idx = 0;
			time_idx++;
    }
		
    // index of first time bin included in inner product
    size_t idx1 = time_idx-mid_idx;
    if(time_idx<mid_idx){idx1=0;}
		else if(idx1>l_ft-1){continue;}
    // index of last time bin included in inner product
    size_t idx2 = time_idx+mid_idx;
    if(idx2>(l_ft-1)){idx2=l_ft-1;}
		else if(idx2<0){continue;}
    
    // ptr to first time bin included in inner product
		size_t data_ptr_offset = l_ft*phase_center_idx+idx1;
    double *data_ptr_real = data_real+data_ptr_offset;
    double *data_ptr_imag = data_imag+data_ptr_offset;
    // prt to first library time bin included in inner product
		size_t lib_ptr_offset = lib_length*lib_idx+mid_idx-(time_idx-idx1);
    double *lib_ptr_real = matched_sig_lib_real+lib_ptr_offset;
    double *lib_ptr_imag = matched_sig_lib_imag+lib_ptr_offset;
    
    // inner product of column
    double col_inner_prod_real = 0;
    double col_inner_prod_imag = 0;
		
		// trig_arg_mult = 2*M_PI*fc
		double trig_arg = trig_arg_mult*t_delay;
		// Real part of exp(-i*2*M_PI*fc*td)
		double carrier_real = cos(trig_arg);
		// Imaginary part of exp(-i*2*M_PI*fc*td)
		double carrier_imag = -sin(trig_arg);
		
    // loop through time bins in column
		// multiply signal by complex conjugate of matched filter
    for(size_t td_bin = idx1;td_bin<=idx2;td_bin++){
      
      // complex multiplication with complex conj. --> C_i = R0*R1-I0*(-I1) + i*(R0*(-I1)+R1*I0)
      // add product to running column inner product (InnProd[0:i] = sum(C[0]:C[i-1]) + C[i])
      double prod_real = (*data_ptr_real) * *lib_ptr_real;
      prod_real += (*data_ptr_imag) * *lib_ptr_imag;
      double prod_imag = -(*data_ptr_real++) * *lib_ptr_imag++;
      prod_imag += (*data_ptr_imag++) * *lib_ptr_real++;
      
      col_inner_prod_real += prod_real;
      col_inner_prod_imag += prod_imag;
    }
		// multiply by complex conjugate of carrier phase delay
		double col_inner_prod_m_carrier_real = col_inner_prod_real * carrier_real + col_inner_prod_imag * carrier_imag;
		double col_inner_prod_m_carrier_imag = col_inner_prod_imag * carrier_real - col_inner_prod_real * carrier_imag;
    
    // add weighted column inner product total inner product
    inner_prod_real += col_inner_prod_m_carrier_real*k_weight;
    inner_prod_imag += col_inner_prod_m_carrier_imag*k_weight;
    // add kx domain window weight for phase_center-pixel combo to running sum
    sum_k += k_weight;
  }
  
  // divide by sum of the k_weights used on phase centers in inner product
  inner_prod_real /= sum_k;
  inner_prod_imag /= sum_k;
  
  out_real = inner_prod_real;
  out_imag = inner_prod_imag;
}

void sar_pixel_RMFR(double &out_real, double &out_imag, size_t pixel_idx, size_t *support_limits, double *kx_limits){
  
  double *pixel_ptr = pixels+3*pixel_idx;
  double pixel_x = *pixel_ptr;
  double pixel_y = *(pixel_ptr+1);
  double pixel_z = *(pixel_ptr+2);
  
  double sum_k = 0;
  double inner_prod_real = 0;
  double inner_prod_imag = 0;
  
  // ptr to x coordinate of left-most phase center for pixel
  double *phase_center_coor_ptr = phase_center+3*support_limits[0];
  // ptr to x element of nadir vector for right-most phase center for pixel
  double *fcs_x_ptr = fcs_x+3*support_limits[0];
  
	double xi_initial = -1.0;
	
  // find inner product of rx image with target at pixe
  // loop from left phase center limit -> right phase center limit
  for(size_t phase_center_idx = support_limits[0];phase_center_idx<=support_limits[1];phase_center_idx++){
    
    double fcs_vect_x = *fcs_x_ptr++;
    double fcs_vect_y = *fcs_x_ptr++;
    double fcs_vect_z = *fcs_x_ptr++;
    
    //// post-increment ptr to y coord.
    double r_x = (*phase_center_coor_ptr++) - pixel_x;
    //// post-inc ptr to z coord.
    double r_y = (*phase_center_coor_ptr++) - pixel_y;
    //// post-inc ptr to x coord. of next phase center
    double r_z = (*phase_center_coor_ptr++) - pixel_z;
    //
    //// range of pixel to phase center
    //double r = r_x*r_x+r_y*r_y+r_z*r_z;
    //r = sqrt(r);
		
		double r = range_funct(pixel_idx,phase_center_idx,xi_initial);
    
    // angle = acos( dot(range_vect,nadir_vect)/abs(range_vect) )
    double kx = -(r_z*fcs_vect_z + r_y*fcs_vect_y + r_x*fcs_vect_x)/r;
    
    // weight applied by kx domain window
    double k_weight = interp(kx_limits,k_window,N_k,kx);
    
    // time delay from phase_center to pixel and back
    double t_delay = r*2/C;
		
    // row index of peak matched filter response for simulated target at pixel
		double time_idx_double = (t_delay/dt-t0/dt);
    size_t time_idx = (size_t)time_idx_double;
		
    // index of library entry associated with sub-time-step delay
    size_t lib_idx = round(((time_idx_double)-(double)(time_idx))*lib_entries);
    if(lib_idx==lib_entries){
      lib_idx = 0;
			time_idx++;
    }
		
    // index of first time bin included in inner product
    size_t idx1 = time_idx-mid_idx;
    if(time_idx<mid_idx){idx1=0;}
		else if(idx1>l_ft-1){continue;}
    // index of last time bin included in inner product
    size_t idx2 = time_idx+mid_idx;
    if(idx2>(l_ft-1)){idx2=l_ft-1;}
		else if(idx2<0){continue;}
    
    // ptr to first time bin included in inner product
		size_t data_ptr_offset = l_ft*phase_center_idx+idx1;
    double *data_ptr_real = data_real+data_ptr_offset;
    double *data_ptr_imag = data_imag+data_ptr_offset;
    // prt to first library time bin included in inner product
    double *lib_ptr_real = matched_sig_lib_real+lib_length*lib_idx+mid_idx-(time_idx-idx1);
    
    // inner product of column
    double col_inner_prod_real = 0;
    double col_inner_prod_imag = 0;
		
		// trig_arg_mult = 2*M_PI*fc
		double trig_arg = trig_arg_mult*t_delay;
		// Real part of exp(-i*2*M_PI*fc*td)
		double carrier_real = cos(trig_arg);
		// Imaginary part of exp(-i*2*M_PI*fc*td)
		double carrier_imag = -sin(trig_arg);
		
    // loop through time bins in column
		// multiply signal by complex conjugate of matched filter
    for(size_t td_bin = idx1;td_bin<=idx2;td_bin++){
      
      // complex multiplication with complex conj. --> C_i = R0*R1-I0*(-I1) + i*(R0*(-I1)+R1*I0)
      // add product to running column inner product (InnProd[0:i] = sum(C[0]:C[i-1]) + C[i])
      double prod_real = (*data_ptr_real++) * *lib_ptr_real;
      double prod_imag = (*data_ptr_imag++) * *lib_ptr_real++;
      
      col_inner_prod_real += prod_real;
      col_inner_prod_imag += prod_imag;
    }
		// multiply by complex conjugate of carrier phase delay
		double col_inner_prod_m_carrier_real = col_inner_prod_real * carrier_real + col_inner_prod_imag * carrier_imag;
		double col_inner_prod_m_carrier_imag = col_inner_prod_imag * carrier_real - col_inner_prod_real * carrier_imag;
    
    // add weighted column inner product total inner product
    inner_prod_real += col_inner_prod_m_carrier_real*k_weight;
    inner_prod_imag += col_inner_prod_m_carrier_imag*k_weight;
    // add kx domain window weight for phase_center-pixel combo to running sum
    sum_k += k_weight;
  }
  
  // divide by sum of the k_weights used on phase centers in inner product
  inner_prod_real /= sum_k;
  inner_prod_imag /= sum_k;
  
  out_real = inner_prod_real;
  out_imag = inner_prod_imag;
}

// interpolates point y(xq) from function y(x)
// axis x exist from x_lim[0]<=x<=x_lim[1] and has N equally spaced points
double interp(double * x_lim,double *y, size_t N, double xq){
  
  double N_dx = (N-1)*(xq-x_lim[0])/(x_lim[1]-x_lim[0]);
  size_t x_idx_1 = (size_t)N_dx;
  double x_remainder = N_dx-double(x_idx_1);
  //double m = (y[x_idx_1+1]-y[x_idx_1])/dx;
  return y[x_idx_1] + x_remainder*(y[x_idx_1+1]-y[x_idx_1]);
}

double getRange(size_t pixel_idx,size_t phase_center_idx,double & xi_initial){
  double *pixel_ptr = pixels+3*pixel_idx;
  double x_pixel = *pixel_ptr;
  double y_pixel = *(pixel_ptr+1);
  double z_pixel = *(pixel_ptr+2);
  
	double * phase_center_ptr = phase_center+3*phase_center_idx;
  double x_pc = *phase_center_ptr;
  double y_pc = *(phase_center_ptr+1);
  double z_pc = *(phase_center_ptr+2);
	
	return sqrt(sqr(x_pixel-x_pc)+sqr(y_pixel-y_pc)+sqr(z_pixel-z_pc));
}

double getRefractionRange(size_t pixel_idx,size_t phase_center_idx,double & xi_initial){
	
  double *pixel_ptr = pixels+3*pixel_idx;
  double x_pixel = *pixel_ptr;
  double y_pixel = *(pixel_ptr+1);
  double z_pixel = *(pixel_ptr+2);
  
	double * phase_center_ptr = phase_center+3*phase_center_idx;
  double x_pc = *phase_center_ptr;
  double y_pc = *(phase_center_ptr+1);
  double z_pc = *(phase_center_ptr+2);
	
	// slow time column of pixel
	size_t output_along_track_idx = pixel_idx/pixel_rows;
	// along track coordinate of pixel
	double xt = output_along_track[output_along_track_idx];
	
	size_t pc_idx = 3*output_along_track_idx;
	
	// Could be named much better
	double pc_x = output_pos[pc_idx];
	double pc_y = output_pos[pc_idx+1];
	double pc_z = output_pos[pc_idx+2];
	// distance from coordesponding output along track phase center position and pixel position
	double yt = -sqrt(sqr(pc_x-x_pixel)+sqr(pc_y-y_pixel)+sqr(pc_z-z_pixel));
	
	// along track coordinate of phase center
	double xr = along_track[phase_center_idx];
	
	// distance from measuring phase center and phase center cooresponding to pixel
	double pc_dist = sqrt(sqr(pc_x-x_pc)+sqr(pc_y-y_pc)+sqr(pc_z-z_pc));
	double along_track_dist = abs(xt - xr);
	
	// ratio between actual dist. and along-track dist.
	double alpha = pc_dist/along_track_dist;
	
	if(alpha==0 || pc_dist==0){
		alpha = 1;
	}
  
  // initial estimate for point of intersection
  double x_surface;
	if(!isnan(xi_initial)){
		if(xi_initial < 0){
			// NaN is returned if there is no intersection with surface or if surface intersection is below pixel
			// In either case, it is assumed that all following phase centers evaluated for this pixel would also return NaN
			x_surface = getSurfaceIntersectionNewton(xr,xt,yt);
		}else{
			x_surface = xi_initial;
		}
	}
	
	if(isnan(x_surface)){
		return sqrt(sqr(x_pixel-x_pc)+sqr(y_pixel-y_pc)+sqr(z_pixel-z_pc));
	}
	
	double r = NAN;
	double r_last = NAN;
	
	for(int i = 0;i<15;i++){
	
		size_t x_surface_idx = getAlongTrackIdx(x_surface);
		// solve for surface y(x)
		// ptr to 0th coefficient of section
		double xd = x_surface-along_track[x_surface_idx];
		double *surface_coeff_ptr = surf_poly+(poly_n+1)*(x_surface_idx+1);
		// y = a0
		double y_surface = *--surface_coeff_ptr;
		// y = a0 + a1*x
		y_surface += *--surface_coeff_ptr*xd;
		// dy/dx = a1;
		double surface_der = *surface_coeff_ptr;
		double surface_der_2 = 0;
		
		double xs[poly_n+1];
		xs[0] = 0;
		xs[1] = xd;
		
		double factorial = 1;
		for(int n = 2;n<=poly_n;n++){
			xs[n] = xs[n-1]*xd;
			y_surface += *--surface_coeff_ptr * xs[n];
			surface_der += (double)n * *surface_coeff_ptr * xs[n-1];
			factorial *= (double)n;
			surface_der_2 += factorial * *surface_coeff_ptr * xs[n-2];
		}
		
		double vr[2] = {xr - x_surface,-y_surface};
		double norm_vr = sqrt(sqr(vr[0])+sqr(vr[1]));
		double vt[2] = {xt - x_surface,yt - y_surface};
		double norm_vt = sqrt(sqr(vt[0])+sqr(vt[1]));
		
		r_last = r;
		r = norm_vr+norm_vt;
		
		if(abs(r-r_last)<sigma_r){
			// save x coordinate intersection
			xi_initial = x_surface;
			// scale the surface intersection x coordinate
      x_surface = alpha*x_surface;
      double r_r = sqrt(sqr(xr-x_surface)+sqr(y_surface));
			double r_t = sqrt(sqr(xt-x_surface)+sqr(vt[1]));
			double r = n0*r_r+n1*r_t;
			return r;
		}
		
		double norm_n = sqrt(1+sqr(surface_der));
		
		double Dr = (-surface_der*vr[0]+vr[1])/norm_n/norm_vr;
		double Dt = (surface_der*vt[0]-vt[1])/norm_n/norm_vt;
				
		double sin_angr = 0;
		double sin_angt = 0;
		if(Dr<1.0 && Dr>-1.0){
			sin_angr = sqrt(1-sqr(Dr));
		}
		if(Dt<1.0 && Dt>-1.0){
			sin_angt = sqrt(1-sqr(Dt));
		}		
		
		if((vr[0]+vr[1]*surface_der)>0){
			if(vt[0]+vt[1]*surface_der>0){
				if(norm_vr>norm_vt){
					sin_angt = -sin_angt;
				}else{
					sin_angr = -sin_angr;
				}
			}
		}else if(vt[0]+vt[1]*surface_der<0){
			if(norm_vr>norm_vt){
				sin_angt = -sin_angt;
			}else{
				sin_angr = -sin_angr;
			}
		}
		
		double Cost = n0*sin_angr-n1*sin_angt;
		
		if(Cost==0){
			// save x coordinate intersection
			xi_initial = x_surface;
			// scale the surface intersection x coordinate
      x_surface = alpha*x_surface;
      double r_r = sqrt(sqr(xr-x_surface)+sqr(y_surface));
			double r_t = sqrt(sqr(xt-x_surface)+sqr(vt[1]));
			double r = n0*r_r+n1*r_t;
			return r;
		}
		
		double der_2_over_norm_n = surface_der_2/norm_n;
		//double arg = surface_der*surface_der_2/sqr(norm_n);
		double arg = surface_der*der_2_over_norm_n/norm_n;
		//double dCdx = n0*Dr/sin_angr*(vr[0]*surface_der_2/norm_n/norm_vr + Dr*(surface_der*surface_der_2/sqr(norm_n)-(vr[0]+surface_der*vr[1])/sqr(norm_vr)))
		//           - n1*Dt/sin_angt*(-vt[0]*surface_der_2/norm_n/norm_vt + Dt*(surface_der*surface_der_2/sqr(norm_n)-(vt[0]+surface_der*vt[1])/sqr(norm_vt)));
		
		double dCdx;
		if(sin_angr==0){
			dCdx = n0*norm_n/norm_vr
					- n1*Dt/sin_angt*(-vt[0]*der_2_over_norm_n/norm_vt + Dt*(arg-(vt[0]+surface_der*vt[1])/sqr(norm_vt)));
		
		}else if(sin_angt==0){
			dCdx = n0*Dr/sin_angr*(vr[0]*der_2_over_norm_n/norm_vr + Dr*(arg-(vr[0]+surface_der*vr[1])/sqr(norm_vr)))
					+ n1*norm_n/norm_vt;
		
		}else{
			dCdx = n0*Dr/sin_angr*(vr[0]*der_2_over_norm_n/norm_vr + Dr*(arg-(vr[0]+surface_der*vr[1])/sqr(norm_vr)))
					- n1*Dt/sin_angt*(-vt[0]*der_2_over_norm_n/norm_vt + Dt*(arg-(vt[0]+surface_der*vt[1])/sqr(norm_vt)));
		}
		//double dCdx = n0*Dr/sin_angr*(vr[0]*der_2_over_norm_n/norm_vr + Dr*(arg-norm_n*sin_angr/norm_vr))
		//           - n1*Dt/sin_angt*(-vt[0]*der_2_over_norm_n/norm_vt + Dt*(arg-norm_n*sin_angt/norm_vt));
					
		double dx = -Cost/dCdx;
		
		if(x_surface==along_track[0]){
			if(dx<0){
				xi_initial = NAN;
				return NAN;
			}
		}else if(x_surface==along_track[l_st-1]){
			if(dx>0){
				xi_initial = NAN;
				return NAN;
			}
		}
		
		//if(Cost>0){
		//	if(to_left){
		//		lims[0] = x_surface;
		//	}else{
		//		lims[1] = x_surface;
		//	}
		//}else{
		//	if(to_left){
		//		lims[1] = x_surface;
		//	}else{
		//		lims[0] = x_surface;
		//	}
		//}
		
		x_surface += dx;
		
		//if(x_surface>lims[1] || x_surface<lims[2]){
		//	x_surface = (lims[0]+lims[1])*.5;
		//}
	}
	xi_initial = NAN;
	return NAN;
}

double getSurfaceIntersectionNewton(double xr,double xt,double yt){
  
	if(yt>=surf_max){
		return NAN;
	}
	
	if(xr==xt){
		if(surf_along_track[getAlongTrackIdx(xr)]<yt){
			return NAN;
		}
		return xr;
	}
	
	// sloped of line from phase center to pixel
  double m = yt/(xt-xr);
	// y-intesect  of line from phase center to pixel
  double b = -m*xr;
  
	// limits on search
  size_t a_lim = 0;
  size_t b_lim = l_st-1;
  
  // initial estimate with line surface estimation
  //double x_surface = (surf_line[1]-b)/(m-surf_line[0]);
	double x_surface = xt;
  
	// should ensure initially passes if statement
  double x_surface_old = x_surface+x_surface;
	
	double d_along_track_est = (along_track[l_st-1] - along_track[0])/(l_st-1);
	double threshold = d_along_track_est*.001;
	
	double y_surface;
	
	// stopping condition?????
  while(abs(x_surface_old-x_surface)>threshold){
		
		size_t x_surface_idx = getAlongTrackIdx(x_surface);
    // solve for surface y(x)
		// ptr to 0th coefficient of section
		double xd = x_surface-along_track[x_surface_idx];
    double *surface_coeff_ptr = surf_poly+(poly_n+1)*(x_surface_idx+1);
		// y = a0
    y_surface = *--surface_coeff_ptr;
		// y = a0 + a1*x
		y_surface += *--surface_coeff_ptr*(xd);
		// dy/dx = a1;
    double surface_der = *surface_coeff_ptr;
		
		double xs[poly_n+1];
		xs[1] = xd;
		
    for(int i = 2;i<=poly_n;i++){
			xs[i] = xs[i-1]*xd;
      y_surface += *--surface_coeff_ptr * xs[i];
			surface_der += (double)i * *surface_coeff_ptr * xs[i-1];
    }
    
		// Cost at x
    double c = y_surface - m*x_surface - b;
		
		if(c==0.0){
			return x_surface;
		}
    
    // cost derivative at x
    double dcdx = surface_der-m;
    
    // Calculate step size from x_surface;
    double dx = -c/dcdx;
    
		//double x_surface_new = x_surface+dx;
		//if(x_surface_new == x_surface_old){
		//	x_surface_new = (x_surface_old+x_surface)*0.5;
		//}
		
		if(abs(surface_der-m)<abs(m)*0.1 || surface_der/m>1.0){
			if(abs(dx)<d_along_track_est){
				break;
			}
			dx = d_along_track_est;
			if(m<0){
				if(c>0){
					//left
					dx = -dx;
				}
			}else{
				if(c<0){
					//left
					dx = -dx;
				}
			}
		}
		
		x_surface_old = x_surface;
		// New estimate
		//x_surface = x_surface_new;
		x_surface = x_surface+dx;
		
    // If new estimate, x, is right of x_surface
    if(dx>0){
		//
		//	// No intersection with surface
			if(x_surface==along_track[l_st-1]){
				return NAN;
			}
		//	
    //  // left limit becomes x_surface
    //  a_lim = x_surface;
		//	
    //  // If Newton search has exceeded known limits
    //  if(x_surface>b_lim){
    //    x_surface = (a_lim+b_lim)*0.5;
    //  }
    //  
    //  // If new estimate, x, is left of x_surface
    }else if(dx<0){
		//
		//	// No intersection with surface
			if(x_surface==along_track[0]){
				return NAN;
			}
			
      // right limit becomes x_surface
    //  b_lim = x_surface;
      
      // If Newton search has exceeded known limits
    //  if(x_surface<a_lim){
    //    x_surface = (a_lim+b_lim)*0.5;
    //  }
    }		
  }
	// if intersection point is below target
	if(y_surface<yt){
		return NAN;
	}
	return x_surface;
}

// getOutputAlongTrackIdx?????
size_t getAlongTrackIdx(double &a){
	
	if(isnan(a)){
    mexErrMsgTxt("Input to getAlongTrackIdx is NAN.");
	}

	size_t lims[2] = {0,l_st-1};

	if(a<0){
		a = along_track[0];
		return 0;
	}else if(a>=along_track[l_st-1]){
		a = along_track[l_st-1];
		return l_st-2;
	}

	double m_est = (along_track[l_st-1] - along_track[0])/(l_st-1);
	size_t a_idx = a/m_est;
	
	if(along_track[a_idx]<=a){
		if(along_track[a_idx+1]>a){
			return a_idx;
		}
		lims[0] = a_idx;
	}else{
		lims[1] = a_idx;
	}
	
	size_t a_idx_previous = l_st;
	while(a_idx!=a_idx_previous){
		a_idx_previous = a_idx;
		
		double da_est = along_track[a_idx+1] - along_track[a_idx];
		a_idx += (a-along_track[a_idx])/da_est;
		
		if(a_idx<=lims[0] || a_idx>=lims[1]){
			a_idx = (lims[0]+lims[1])*0.5;
		}
		
		if(along_track[a_idx]<=a){
			if(along_track[a_idx+1]>a){
				return a_idx;
			}
			lims[0] = a_idx;
		}else{
			lims[1] = a_idx;
		}
	}
	return a_idx;
}

double getSurfaceIntersectionGolden(double xr,double xt,double yt){
  
	if(yt>=surf_max){
		return NAN;
	}
	
	if(xr==xt){
		if(surf_along_track[getAlongTrackIdx(xr)]<yt){
			return NAN;
		}
		return xr;
	}
	
	double d_along_track_est = (along_track[l_st-1] - along_track[0])/(l_st-1);
	double threshold = d_along_track_est*.001;
	
	// sloped of line from phase center to pixel
  double m = yt/(xt-xr);
	// y-intesect  of line from phase center to pixel
  double b_int = -m*xr;
	
	double a = xr;
	double b = xt;
	if(xr>xt){
		a = xt;
		b = xr;
	}
	
	double w = b-a;
	
	double y = b-w*tau;
	
	
	size_t y_idx = getAlongTrackIdx(y);
	// solve for surface y(x)
	// ptr to 0th coefficient of section
	double yd = y-along_track[y_idx];
	double *y_coeff_ptr = surf_poly+(poly_n+1)*(y_idx+1);
	// y = a0
	double yy = *--y_coeff_ptr;
	// y = a0 + a1*x
	yy += *--y_coeff_ptr*(yd);
	// dy/dx = a1;
	double y_der = *y_coeff_ptr;
	
	double ys[poly_n+1];
	ys[1] = yd;
	
	for(int i = 2;i<=poly_n;i++){
		ys[i] = ys[i-1]*yd;
		yy += *--y_coeff_ptr * ys[i];
		y_der += (double)i * *y_coeff_ptr * ys[i-1];
	}
	
	double cy = yy - m*y - b_int;
	cy = cy*cy;
			
	int from_right = 1;
	
	while(w>threshold){
		double x;
		if(from_right){
			x = a+w*tau;
		}else{
			x = b-w*tau;
		}
		
		size_t x_idx = getAlongTrackIdx(x);
		// solve for surface y(x)
		// ptr to 0th coefficient of section
		double xd = x-along_track[x_idx];
		double *x_coeff_ptr = surf_poly+(poly_n+1)*(x_idx+1);
		// y = a0
		double xy = *--x_coeff_ptr;
		// y = a0 + a1*x
		xy += *--x_coeff_ptr*(xd);
		// dy/dx = a1;
		double x_der = *x_coeff_ptr;
		
		double xs[poly_n+1];
		xs[1] = xd;
		
		for(int i = 2;i<=poly_n;i++){
			xs[i] = xs[i-1]*xd;
			xy += *--x_coeff_ptr * xs[i];
			x_der += (double)i * *x_coeff_ptr * xs[i-1];
		}
			
		double cx = xy - m*x - b_int;
		cx = cx*cx;
		
		if(cx<cy){
			y = x;
			cy = cx;
			if(from_right){
				a = y;
			}else{
				b = y;
			}
		}else{
			if(from_right){
				b = x;
				from_right = 0;
			}else{
				a = x;
				from_right = 1;
			}
		}
		
		w = b-a;
	}
	return (b+a)/2;
}
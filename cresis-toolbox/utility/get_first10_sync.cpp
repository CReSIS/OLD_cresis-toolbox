/*
 * Compile with:
 *   mex get_first10_sync.c
 */

#include "mex.h"
#include "stdio.h"

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  const int SYNCS_TO_LOAD = 10;
  const int DATA_BLOCK_SIZE = 65536;
  
  /* Check for proper number of arguments. */
  if(nrhs!=2) {
    mexErrMsgTxt("Two inputs required (MCoRDS filename and first byte).");
  } else if(nlhs>1) {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  mwSize mrows,ncols;
  
  /* The first input must be a string (row vector of chars).*/
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsChar(prhs[0]) || !(mrows==1) ) {
    mexErrMsgTxt("First input must be a string.");
  }
  
  /* The second input must be a scalar .*/
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])
      || !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Second input must be a real scalar double.");
  }
  
  // =======================================================================
  // Get filename
  char *filename;
  /* copy the string data from prhs[0] into a C string filename.    */
  filename = mxArrayToString(prhs[0]);
  if (filename == NULL) 
  {
    mexErrMsgTxt("Could not convert input to string.");
  }
  
  // =======================================================================
  // Get first byte to load in
  int first_byte;
  first_byte = (int)*mxGetPr(prhs[1]);
  //mexPrintf("%d\n", first_byte);
  
  // =======================================================================
  // Get output pointer
  plhs[0] = mxCreateDoubleMatrix(1, SYNCS_TO_LOAD, mxREAL);
  double *sync_pos;
  sync_pos = mxGetPr(plhs[0]);
  
  // =======================================================================
  // Open file
  FILE *fid;
  fid = fopen(filename, "rb");
  if (fid == NULL)
  {
    mexPrintf("Tried to open %s\n", filename);
    mexErrMsgTxt("Failed to open file.");
  }
  // =======================================================================
  // Skip past bad bytes if user requests
  fseek(fid, first_byte, SEEK_SET);
  
  // =======================================================================
  // Look for the first 10 (SYNCS_TO_LOAD) 0xDEADBEEF markers
  unsigned char data[DATA_BLOCK_SIZE];
  int out_idx = 0;
  int num_read = DATA_BLOCK_SIZE;
  while (~feof(fid) && num_read == DATA_BLOCK_SIZE && out_idx < SYNCS_TO_LOAD)
  {      
    int cur_pos = ftell(fid);
    num_read = fread(data, sizeof(char), DATA_BLOCK_SIZE, fid);
    for (int idx = 0; idx < num_read-3; idx++)
    {
      if (data[idx] == 0xDE && data[idx+1] == 0xAD && data[idx+2] == 0xBE && data[idx+3] == 0xEF)
      {
        sync_pos[out_idx] = cur_pos + idx;
        out_idx++;
        if (out_idx == SYNCS_TO_LOAD)
        {
          break;
        }
      }
    }
    fseek(fid, -3, SEEK_CUR);
  }
  
  fclose(fid);
  
}

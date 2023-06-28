/*
 * file : delvar.c
 *
 * delvar MAT-file variable
 *
 * Downloaded from Matlab Central.
 */

#include "mat.h"
#include "mex.h"

/*
 * usage of MEX-file
 */
void
printUsage()
{
  mexPrintf("Usage: %s MAT-file variable\n",mexFunctionName());
}

void
mexFunction( int nlhs,
             mxArray *plhs[],
             int nrhs,
             const mxArray *prhs[]
           )
{
  MATFile *mfp;
  //mxArray *status;
  char *filename;
  char *variable;
  int buffersize;
  int status;

  /*
   * the value returned by the left-hand side is the status. return 1 if
   * there is a failure. return 0 if there is success (variable was
   * removed from MAT-file.
   */
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  mxGetPr(plhs[0])[0] = 1;

  /*
   * error checking for input arguments
   */
  if (nrhs!=2) {
    printUsage();
    return;
  }

  if ( (!mxIsChar(prhs[0])) ||
       (!mxIsChar(prhs[1])) ) {
    // Usage of function incorrect
    printUsage();
    mxGetPr(plhs[0])[0]=1;
    return;
  }

  /*
   * get filename to open
   */
  buffersize=mxGetM(prhs[0])*mxGetN(prhs[0])+1;
  filename=(char*)mxCalloc(buffersize,sizeof(char));
  mxGetString(prhs[0],filename,buffersize);

  /*
   * open MAT-file
   */
  mfp=matOpen(filename,"u+");
  if (mfp==(MATFile *)NULL) {
    // Failed to open file
    mxFree(filename);
    mxGetPr(plhs[0])[0]=2;
    return;
  }

  /*
   * get variable to delete
   */
  buffersize=mxGetM(prhs[1])*mxGetN(prhs[1])+1;
  variable=(char*)mxCalloc(buffersize,sizeof(char));
  mxGetString(prhs[1],variable,buffersize);
    
  /*
   * delete variable from MAT-file
   */
  status=matDeleteVariable(mfp,variable);
  if (status!=0) {
    // Failed to delete variable
    mxFree(filename);
    mxFree(variable);
    matClose(mfp);
    mxGetPr(plhs[0])[0]=10+status;
    return;
  }

  /*
   * cleanup
   */
  mxFree(filename);
  mxFree(variable);
  matClose(mfp);

  /*
   * change return status to success
   */
  mxGetPr(plhs[0])[0]=0;
}

// Compile with mex -v -largeArrayDims basic_load_hdr_mex.cpp
//   -v: verbose
//   -largeArrayDims: required for 64 bit
//
// Assumption is that files are 2^31 bytes so the 32 bit integer can be
// used for record offsets

#include <stdio.h>
#include <stdlib.h>
#include "mat.h"
#include "mex.h"

/*
 * usage of MEX-file
 */
void
printUsage()
{
  mexPrintf("[file_size offset field1 ... fieldN] = %s(fn, frame_sync, field_offsets, field_type);\n", mexFunctionName());
  mexPrintf("\n");
  mexPrintf(" Loads header fields from radar data files that use frame sync. Function\n");
  mexPrintf(" used by create_segment_raw_file_list_v2.m. Assumes ieee-be file byte order and ieee-le host byte order\n");
  mexPrintf("\n");
  mexPrintf(" fn = string containing filename to load\n");
  mexPrintf(" frame_sync = One of the following frame sync values:\n");
  mexPrintf("    hex2dec('DEADBEEF'), hex2dec('1ACFFC1D'), hex2dec('BADA55E5')\n");
  mexPrintf(" field_offsets = row vector containing the offsets to each 32 bit sized header\n");
  mexPrintf("    field that will be returned. Size of varargout is equal to the length of\n");
  mexPrintf("    this vector.  For example, [1 2 3], would return the 3 32 bit values\n");
  mexPrintf("    preceding each frame sync.\n");
  mexPrintf(" file_type is cell array of types equal in length to field_offsets\n");
  mexPrintf("\n");
  mexPrintf("file_size is in bytes, offset is vector containing byte offset to each frame sync\n");
  mexPrintf("field1 through fieldN are from the file\n");
  mexPrintf("\n");
  mexPrintf("[file_size offset seconds fractions] = %s(fn, hex2dec('DEADBEEF'), [2 3], {uint32(1) uint32(1)})\n", mexFunctionName());

}

unsigned short swap_bytes_16bit(unsigned short val)
{
  return val>>8 | val<<8;
}

unsigned int swap_bytes_32bit(unsigned int val)
{
  return val>>24 | val>>8&0x0000FF00 | val<<8&0x00FF0000 | val<<24;
}

unsigned long long swap_bytes_64bit(unsigned long long val)
{
  return val>>56 | val>>40&0x000000000000FF00 | val>>24&0x0000000000FF0000 | val>>8&0x00000000FF000000 
    | val<<8&0x000000FF00000000 | val<<24&0x0000FF0000000000 | val<<40&0x00FF000000000000 | val<<56;
}

void
mexFunction( int nlhs,
             mxArray *plhs[],
             int nrhs,
             const mxArray *prhs[]
           )
{
  char *filename;
  ptrdiff_t buffersize;

  if (nlhs<1) {
    printUsage();
    return;
  }

  /*
   * the value returned by the left-hand side is the status. return 1 if
   * there is a failure. return 0 if there is success (variable was
   * removed from MAT-file.
   */
  plhs[0]=mxCreateDoubleMatrix(1,1,mxREAL);
  mxGetPr(plhs[0])[0] = 0;

  /*
   * error checking for input arguments
   */
  if (nrhs!=4) {
    printUsage();
    return;
  }

  if ( (!mxIsChar(prhs[0])) ||
       (!mxIsClass(prhs[1], "uint32")) ||
       mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1 ||
       (!mxIsClass(prhs[2], "uint32")) ||
       (!mxIsClass(prhs[3], "cell")) ) {
    // Usage of function incorrect
    printUsage();
    return;
  }

  ptrdiff_t num_rows = mxGetM(prhs[2]);
  ptrdiff_t num_cols = mxGetN(prhs[2]);
  ptrdiff_t num_outputs;
  if (num_rows == 0) {
    num_outputs = 0;
  } else if (num_rows == 1) {
    num_outputs = num_cols;
  } else {
    // Usage of function incorrect
    printUsage();
    return;
  }

  if (num_outputs != nlhs - 2 ||
    mxGetM(prhs[2]) !=  mxGetM(prhs[3]) ||
    mxGetN(prhs[2]) !=  mxGetN(prhs[3]) )
  {
    // Usage of function incorrect
    printUsage();
    return;
  }

  // Allocate later
  plhs[1]=mxCreateNumericMatrix(0,0,mxUINT32_CLASS,mxREAL);
  for (int idx=2; idx <= num_outputs+1; idx++)
  {
    mxArray *class_type = mxGetCell(prhs[3], idx-2);

    if (class_type == NULL)
    {
      printUsage();
      return;
    }

    mxClassID class_type_ID = mxGetClassID(class_type);

    bool good_type = false;
    switch (class_type_ID)  {
      case mxINT8_CLASS:
        good_type = true;
        break; 
      case mxUINT8_CLASS:
        good_type = true;
        break;
      case mxINT16_CLASS:
        good_type = true;
        break;
      case mxUINT16_CLASS:
        good_type = true;
        break;
      case mxINT32_CLASS:
        good_type = true;
        break;
      case mxUINT32_CLASS:
        good_type = true;
        break;
      case mxINT64_CLASS:
        good_type = true;
        break;
      case mxUINT64_CLASS:
        good_type = true;
        break;
      case mxSINGLE_CLASS:
        good_type = true;
        break; 
      case mxDOUBLE_CLASS:
        good_type = true;
        break;
      default:
        break;
    }
    if (!good_type)
    {
      printUsage();
      return;
    }

    plhs[idx]=mxCreateNumericMatrix(0,0,class_type_ID,mxREAL);
  }

  /*
   * get filename to open
   */
  buffersize=mxGetM(prhs[0])*mxGetN(prhs[0])+1;
  filename=(char*)mxCalloc(buffersize,sizeof(char));
  mxGetString(prhs[0],filename,buffersize);

  // Get frame_sync
  unsigned int frame_sync = ((unsigned int*)mxGetPr(prhs[1]))[0];

  // ======================================================================
  // Try to keep the last record in the file unless it is too short to include
  // all the header information needed by this function. To do this, we simply
  // avoid using these last "min_complete_record_size" bytes.  A "complete"
  // record for this function is one that contains all the dynamic output
  // variables that are requested (i.e. not all the data in the record...
  // just the outputs we need right now).
  ptrdiff_t min_complete_record_size = 4; // 4 for the frame sync
  
  // Get dynamic output variable offsets and sizes, update 
  // min_complete_record_size if needed
  for (ptrdiff_t out_idx = 0; out_idx < num_outputs; out_idx++)
  {
    mxArray *class_type = mxGetCell(prhs[3], out_idx);
    mxClassID class_type_ID = mxGetClassID(class_type);
    ptrdiff_t field_offset = ((unsigned int *)mxGetPr(prhs[2]))[out_idx];

    ptrdiff_t field_offset_and_size;

    switch (class_type_ID)  {
      case mxINT8_CLASS:
        field_offset_and_size = field_offset + sizeof(char);
        break; 
      case mxUINT8_CLASS:
        field_offset_and_size = field_offset + sizeof(unsigned char);
        break;
      case mxINT16_CLASS:
        field_offset_and_size = field_offset + sizeof(short);
        break;
      case mxUINT16_CLASS:
        field_offset_and_size = field_offset + sizeof(unsigned short);
        break;
      case mxINT32_CLASS:
        field_offset_and_size = field_offset + sizeof(int);
        break;
      case mxUINT32_CLASS:
        field_offset_and_size = field_offset + sizeof(unsigned int);
        break;
      case mxINT64_CLASS:
        field_offset_and_size = field_offset + sizeof(long long);
        break;
      case mxUINT64_CLASS:
        field_offset_and_size = field_offset + sizeof(unsigned long long);
        break;
      case mxSINGLE_CLASS:
        field_offset_and_size = field_offset + sizeof(float);
        break; 
      case mxDOUBLE_CLASS:
        field_offset_and_size = field_offset + sizeof(double);
        break;
    }
    
    if (4+field_offset_and_size > min_complete_record_size)
    {
      // Dynamic variable requires bigger minimum complete record size
      min_complete_record_size = 4+field_offset_and_size;
    }
  }
  // mexPrintf("min_complete_record_size = %d\n", min_complete_record_size); // DEBUG
  
  // ======================================================================
  // Open the data file and read all the contents
  FILE *fptr;
  fptr = fopen(filename,"r");
  if (fptr == 0)
  {
    mxFree(filename);
    return;
  }

  fseek(fptr, 0, SEEK_END);
  ptrdiff_t file_size = ftell(fptr);
  fseek(fptr, 0, SEEK_SET);

  unsigned char *data;
  data = (unsigned char *)mxMalloc(file_size);

  fread(data, 1, file_size, fptr);

  fclose(fptr);

  // ======================================================================
  // Search for all the frame syncs in the file and record their byte offset
  int *offset;
  ptrdiff_t offset_size = 20000;
  offset = (int *)mxRealloc(NULL,offset_size * sizeof(int));

  unsigned char fs1 = frame_sync >> 24;
  unsigned char fs2 = (frame_sync >> 16) % 256;
  unsigned char fs3 = (frame_sync >> 8) % 256;
  unsigned char fs4 = frame_sync % 256;
  //mexPrintf("%u\n", frame_sync); // DEBUG
  //mexPrintf("%x %x %x %x\n", fs1, fs2, fs3, fs4); // DEBUG
  mwSize offset_idx = 0;
  for (ptrdiff_t idx = 0; idx < file_size-min_complete_record_size; idx++)
  {
    if (data[idx] == fs1 && data[idx+1] == fs2 && data[idx+2] == fs3 && data[idx+3] == fs4)
    {
      if (offset_idx >= offset_size)
      {
        offset_size = 2*offset_size;
        offset = (int *)mxRealloc(offset,offset_size*sizeof(int));
      }
      offset[offset_idx] = idx;
      offset_idx++;
    }
  }
  
  const mwSize num_records = offset_idx;

  // ======================================================================
  // Get dynamic output variables in each of the records that was found
  for (ptrdiff_t out_idx = 0; out_idx < num_outputs; out_idx++)
  {
    mxArray *class_type = mxGetCell(prhs[3], out_idx);
    mxClassID class_type_ID = mxGetClassID(class_type);

    void *outvar;

    switch (class_type_ID)  {
      case mxINT8_CLASS:
        outvar = mxMalloc(num_records*sizeof(char));
        break; 
      case mxUINT8_CLASS:
        outvar = mxMalloc(num_records*sizeof(unsigned char));
        break;
      case mxINT16_CLASS:
        outvar = mxMalloc(num_records*sizeof(short));
        break;
      case mxUINT16_CLASS:
        outvar = mxMalloc(num_records*sizeof(unsigned short));
        break;
      case mxINT32_CLASS:
        outvar = mxMalloc(num_records*sizeof(int));
        break;
      case mxUINT32_CLASS:
        outvar = mxMalloc(num_records*sizeof(unsigned int));
        break;
      case mxINT64_CLASS:
        outvar = mxMalloc(num_records*sizeof(long long));
        break;
      case mxUINT64_CLASS:
        outvar = mxMalloc(num_records*sizeof(unsigned long long));
        break;
      case mxSINGLE_CLASS:
        outvar = mxMalloc(num_records*sizeof(float));
        break; 
      case mxDOUBLE_CLASS:
        outvar = mxMalloc(num_records*sizeof(double));
        break;
    }
    ptrdiff_t field_offset = ((unsigned int *)mxGetPr(prhs[2]))[out_idx];

    for (ptrdiff_t idx = 0; idx < num_records; idx++)
    {
      switch (class_type_ID)  {
        case mxINT8_CLASS:
          ((char*)outvar)[idx] = ((char*)(data + offset[idx] + field_offset))[0];
          break; 
        case mxUINT8_CLASS:
          ((unsigned char*)outvar)[idx] = ((unsigned char*)(data + offset[idx] + field_offset))[0];
          break;
        case mxINT16_CLASS:
          ((short*)outvar)[idx] = swap_bytes_16bit(((short*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxUINT16_CLASS:
          ((unsigned short*)outvar)[idx] = swap_bytes_16bit(((unsigned short*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxINT32_CLASS:
          ((int*)outvar)[idx] = swap_bytes_32bit(((int*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxUINT32_CLASS:
          ((unsigned int*)outvar)[idx] = swap_bytes_32bit(((unsigned int*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxINT64_CLASS:
          ((long long*)outvar)[idx] = swap_bytes_64bit(((long long*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxUINT64_CLASS:
          ((unsigned long long*)outvar)[idx] = swap_bytes_64bit(((unsigned long long*)(data + offset[idx] + field_offset))[0]);
          break;
        case mxSINGLE_CLASS:
          ((float*)outvar)[idx] = (float)swap_bytes_32bit((int)((float*)(data + offset[idx] + field_offset))[0]);
          break; 
        case mxDOUBLE_CLASS:
          ((double*)outvar)[idx] = (double)swap_bytes_64bit((long long)((double*)(data + offset[idx] + field_offset))[0]);
          break;
      }
    }
    // Assign this dynamic output variable
    mxSetData(plhs[2+out_idx], outvar);
    mxSetM(plhs[2+out_idx], 1);
    mxSetN(plhs[2+out_idx], num_records);
  }

  // ======================================================================
  // cleanup
  mxFree(filename);
  mxFree(data);

  // ======================================================================
  // Assign offset output variable
  mxSetData(plhs[1], offset);
  mxSetM(plhs[1], 1);
  mxSetN(plhs[1], num_records);

  // ======================================================================
  // change return status to success
  mxGetPr(plhs[0])[0]=file_size;
}


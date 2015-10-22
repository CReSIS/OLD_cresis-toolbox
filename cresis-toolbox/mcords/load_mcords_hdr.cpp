// [hdr wfs status] = load_mcords_hdr(fn, param);
//
// Loads MCORDS 2009 Antarctica DC-8 headers.
//   - The function reads one file at a time
//   - It returns enough information through the status structure so that
//     multiple files can
//     easily be concatenated (see create_records_mcords.m for example)
//
// Compile with (-O == optimization, should be on by default):
//   mex -O load_mcords_hdr.cpp
// Your active directory must be the same directory as the file to compile it.
//
// Using .cpp file ending to allow "//" commenting and to allow
// declarations anywhere in the code. Otherwise it is C-code.
//
// INPUTS:
//   fn = filename string
//   param = parameter structure
//     first_byte = first byte to read in (usually 0, except for the
//       first file index, then see example below)
//     fs = sampling frequency of radar, used to interpret counts/tics
//       in header field variables
//     rec_size = record size (see example for how to compute)
//     sync_lock = is the sync locked (i.e. does first_byte point to
//       a valid record)
//     last_bytes = bytes to prepend to the data, first_byte is relative
//       to these
//
// OUTPUTS:   
//   hdr.ver = num_recx1 vector of uint8
//   hdr.seconds = num_recx1 vector of uint32
//   hdr.fractions = num_recx1 vector of uint32
//   hdr.epri = num_recx1 vector of uint32
//   hdr.offset = num_recx1 vector of int32
//     The first element may be negative which means that the start of
//     this record is in the previous file.
//
//   wfs.num_sam = hdr.num_wf x 1 vector of double
//   wfs.which_bits = hdr.num_wf x 1 vector of double
//   wfs.t0 = hdr.num_wf x 1 vector of double
//   wfs.presums = hdr.num_wf x 1 vector of double
//
//   status.last_bytes = last bytes in the file
//   status.sync_lock = current status of the sync
//
// Example to run:
//    syncs = get_first10_sync(fn, 2^25 for file index 0 otherwise 0);
//    rec_size = mean(diff(syncs));
//    first_byte = find(diff(syncs) == rec_size, 1);
//    param.first_byte = syncs(first_byte);
//    param.rec_size = median(diff(syncs));
//    param.sync_lock = 1;
//    param.last_bytes = '';
//    param.fs = 1e9/9;
//    [hdr wfs status] = load_mcords_hdr(fn, param);
//
// Authors: John Paden

//HEADER FILE FORMAT
//Header Size: 4*8 + 16*8 = 160 bytes
//struct header
//{
//  uint32 frame_sync; // aka magic
//  uint32 version;    // file type version (e.g. May 2010 is 2)
//  uint32 seconds;    // # of seconds
//  uint32 fraction;   // # of fs/2 clock cycles
//  uint32 epri;       // effective PRI number
//  uint32 num_wf;     // number of waveforms
//  uint32 which_bits; // ???
//  uint32 dec_cfg;    // decimation config settings
//  struct wf wfs[16];
//}
//struct wf
//{
//  uint32(31:14) tmp;
//  uint32(13:0) num_sam;    // number of samples
//  uint32(31:28) tmp;
//  uint32(27:24) which_bits // whichbits for waveform
//  uint32(23:10) t0 // sample delay in fs clock cycles + 11 us
//  uint32(9:0) presums      // number of presums - 1
//}

#include "mex.h"
#include "stdio.h"
#include "string.h"

// Declare support functions which are defined at the bottom of this file
inline unsigned int big2littleint32(unsigned char *buffer);
int find_next_sync(FILE *fid, int bytes_left);

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  
  // Check for proper number of arguments
  if (nrhs != 2)
  {
    mexErrMsgTxt("Two inputs required (MCoRDS filename, param structure).");
  }
  else if (nlhs != 3)
  {
    mexErrMsgTxt("There must be 3 output arguments.");
  }
  
  mwSize mrows,ncols;
  
  // The first input must be a string (row vector of chars)
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsChar(prhs[0]) || !(mrows==1) ) {
    mexErrMsgTxt("First input must be a string.");
  }
  
  // The second input must be a 1x1 struct
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsStruct(prhs[1]) || !(mrows==1 && ncols==1) ) {
    mexErrMsgTxt("Second input must be a struct.");
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
  
  // Second argument can have first_byte double scalar field
  mxArray *mx_first_byte;
  int first_byte;
  mx_first_byte = mxGetField(prhs[1], 0, "first_byte");
  if (mx_first_byte == NULL || mxIsEmpty(mx_first_byte))
  {
    first_byte = 0;
  }
  else
  {
    first_byte = (int)mxGetPr(mx_first_byte)[0];
  }
  //mexPrintf("%d\n", first_byte);
  
  // Second argument MUST have rec_size double scalar field
  mxArray *mx_rec_size;
  int rec_size;
  mx_rec_size = mxGetField(prhs[1], 0, "rec_size");
  if (mx_rec_size == NULL || mxIsEmpty(mx_rec_size))
  {
    mexErrMsgTxt("Record size must be set.");
  }
  else
  {
    rec_size = (int)mxGetPr(mx_rec_size)[0];
  }
  //mexPrintf("%d\n", rec_size);
  
  // Second argument can have sync_lock double scalar field
  mxArray *mx_sync_lock;
  int sync_lock;
  mx_sync_lock = mxGetField(prhs[1], 0, "sync_lock");
  if (mx_sync_lock == NULL || mxIsEmpty(mx_sync_lock))
  {
    sync_lock = 0;
  }
  else
  {
    sync_lock = (int)mxGetPr(mx_sync_lock)[0];
  }
  //mexPrintf("%d\n", sync_lock);
  
  // Second argument MUST have fs double scalar field
  mxArray *mx_fs;
  double fs;
  mx_fs = mxGetField(prhs[1], 0, "fs");
  if (mx_fs == NULL || mxIsEmpty(mx_fs))
  {
    mexErrMsgTxt("fs must be set.");
  }
  else
  {
    fs = (int)mxGetPr(mx_fs)[0];
  }
  //mexPrintf("%f\n", fs);
  
  // Second argument can have last_bytes uint8 vector field
  mxArray *mx_last_bytes;
  unsigned char *last_bytes;
  int last_bytes_len;
  mx_last_bytes = mxGetField(prhs[1], 0, "last_bytes");
  if (mx_last_bytes == NULL || mxIsEmpty(mx_last_bytes))
  {
    last_bytes = 0;
    last_bytes_len = 0;
  }
  else
  {
    last_bytes = (unsigned char *)mxGetPr(mx_last_bytes);
    last_bytes_len = mxGetNumberOfElements(mx_last_bytes);
  }
  //mexPrintf("%p\n", last_bytes);
  
  // =======================================================================
  // Open file
  FILE *fid;
  //mexPrintf("Opening file %s\n", filename);
  fid = fopen(filename, "rb");
  if (fid == NULL)
  {
    mexErrMsgTxt("Failed to open file.");
  }
  
  // =======================================================================
  // Determine the maximum number of complete records in the file
  fseek(fid, 0, SEEK_END);
  int file_size = ftell(fid);
  int num_rec = (last_bytes_len + file_size - first_byte)/rec_size;
  //mexPrintf("Reading %d records\n", num_rec); // DEBUG
  
  // =======================================================================
  // Get output pointers
  
  //    hdr.ver = uint8
  //    hdr.seconds = uint32
  //    hdr.fractions = uint32
  //    hdr.epri = uint32
  //    hdr.num_wf = uint8
  int num_hdr_fields= 5;
  const char *hdr_fields[] = {"ver","seconds","fractions","epri","offset"};
  plhs[0] = mxCreateStructMatrix(1, 1, num_hdr_fields, hdr_fields);

  mxArray *mx_ver = mxCreateNumericArray(1, &num_rec, mxUINT8_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[0], 0, 0, mx_ver);
  unsigned char *ver = (unsigned char *)mxGetPr(mx_ver);

  mxArray *mx_seconds = mxCreateNumericArray(1, &num_rec, mxUINT32_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[0], 0, 1, mx_seconds);
  unsigned int *seconds = (unsigned int *)mxGetPr(mx_seconds);
  
  mxArray *mx_fractions = mxCreateNumericArray(1, &num_rec, mxUINT32_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[0], 0, 2, mx_fractions);
  unsigned int *fractions = (unsigned int *)mxGetPr(mx_fractions);
  
  mxArray *mx_epri = mxCreateNumericArray(1, &num_rec, mxUINT32_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[0], 0, 3, mx_epri);
  unsigned int *epri = (unsigned int *)mxGetPr(mx_epri);

  mxArray *mx_offset = mxCreateNumericArray(1, &num_rec, mxINT32_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[0], 0, 4, mx_offset);
  int *offset = (int *)mxGetPr(mx_offset);
  
  // last_bytes = uint8 (the bytes of the last incomplete record)
  // sync_lock = double (status of the sync lock)
  int num_status_fields = 2;
  const char *status_fields[] = {"last_bytes","sync_lock"};
  plhs[2] = mxCreateStructMatrix(1, 1, num_status_fields, status_fields);
  mx_sync_lock = mxCreateDoubleScalar(sync_lock);
  mxSetFieldByNumber(plhs[2], 0, 1, mx_sync_lock);
  
  // =======================================================================
  // Go to first record
  fseek(fid, first_byte, SEEK_SET);

  // =======================================================================
  // If sync is lost, find it before reading in records
  if (!sync_lock)
  {
    sync_lock = find_next_sync(fid, file_size - first_byte - rec_size);
    mxGetPr(mx_sync_lock)[0] = sync_lock;
    if (sync_lock)
    {
      mexPrintf("No good records in this file.");
      return;
    }
  }

  // =======================================================================
  // =======================================================================
  // Read in headers
  // =======================================================================
  // =======================================================================
  
  // mexCallMATLAB(0,NULL,0,NULL, "toc"); // DEBUG
  unsigned char *buffer;
  
  // =======================================================================
  // Read in NUM_REC_BLOCK records at a time to speed up loading
  unsigned int NUM_REC_BLOCK = 512;
  buffer = (unsigned char*)mxCalloc(rec_size*NUM_REC_BLOCK, sizeof(char));
  if (buffer == NULL)
  {
    mexErrMsgTxt("Failed to allocate buffer");
  }
  
  int rec;
  for (rec = 0; rec < num_rec && ftell(fid) <= file_size-rec_size 
          && sync_lock;)
  {
    int num_read;
    // Handle end of file case when there are not NUM_REC_BLOCK records left
    if (num_rec - rec > NUM_REC_BLOCK)
    {
      num_read = NUM_REC_BLOCK;
    }
    else
    {
      num_read = num_rec - rec;
    }
    
    // Handle the beginning case when the user may have passed in some
    // bytes from the previous file and we need to prepend those to our
    // result.
    int bytes_read;
    if (rec == 0 && last_bytes_len > 0)
    {
      memcpy(buffer, last_bytes, last_bytes_len);
      bytes_read = fread(buffer + last_bytes_len,1,rec_size*num_read - last_bytes_len,fid);
    }
    else
    {
      bytes_read = fread(buffer,1,rec_size*num_read,fid);
    }

    // Go through each of the num_read blocks and parse their headers
    for (int sub_rec = 0; sub_rec < num_read; sub_rec++, rec++)
    {
      // ==============================================================
      // Check for bad sync marker
      if (big2littleint32(&buffer[sub_rec*rec_size+0*4]) != 0xDEADBEEF)
      {
        // Check to see if 0xDEADBEEF: if not, then find the next one
        // and fseek to it and break
        fseek(fid,-bytes_read+sub_rec*rec_size, SEEK_CUR);
        //mexPrintf("  ==> Lost sync lock (0xDEADBEEF frame sync not found, %d)\n", ftell(fid));
        int bytes_left = file_size - ftell(fid) - rec_size;
        sync_lock = find_next_sync(fid, bytes_left);
        if (!sync_lock)
        {
          //mexPrintf("  ==> No more locks in this file\n");
          mxGetPr(mx_sync_lock)[0] = 1;
        }
        break;
      }
      
      // ==============================================================
      // Read, parse, and interpret header fields
      //  uint32 frame_sync; // aka magic
      //  uint32 version;    // file type version (e.g. May 2010 is 2)
      //  uint32 seconds;    // # of seconds
      //  uint32 fraction;   // # of fs/2 clock cycles
      //  uint32 epri;       // effective PRI number
      //  uint32 num_wf;     // number of waveforms
      //  uint32 which_bits; // ???
      //  uint32 dec_cfg;    // decimation config settings
      // ==============================================================
      // mexPrintf("0x%08x\n", big2littleint32(&buffer[sub_rec*rec_size+0*4])); // DEBUG
      
      offset[rec] = ftell(fid) - bytes_read + sub_rec*rec_size - last_bytes_len;
      if (offset[rec] > file_size-rec_size)
      {
        //mexPrintf("  ==> Too short %d\n", offset[rec]);
        fseek(fid, offset[rec], SEEK_SET);
        break;
      }
      ver[rec] = buffer[sub_rec*rec_size+4+3];
      seconds[rec] = big2littleint32(&buffer[sub_rec*rec_size+2*4]);
      fractions[rec] = big2littleint32(&buffer[sub_rec*rec_size+3*4]);
      epri[rec] = big2littleint32(&buffer[sub_rec*rec_size+4*4]);
      int num_wfs = buffer[sub_rec*rec_size+5*4+3];
      
      // ==============================================================
      // Read in waveform header
      //  uint32(31:14) tmp;
      //  uint32(13:0) num_sam;    // number of samples
      //  uint32(31:28) tmp;
      //  uint32(27:24) which_bits // whichbits for waveform
      //  uint32(23:10) t0 // sample delay in fs clock cycles
      //  uint32(9:0) presums      // number of presums - 1
      if (rec == 0)
      {
        int num_wfs_fields = 4;
        const char *wfs_fields[] = {"num_sam", "which_bits", "t0", "presums"};
        plhs[1] = mxCreateStructMatrix(1, 1, num_wfs_fields, wfs_fields);

        mxArray *mx_num_sam = mxCreateDoubleMatrix(num_wfs, 1, mxREAL);
        mxSetFieldByNumber(plhs[1], 0, 0, mx_num_sam);
        double *num_sam = mxGetPr(mx_num_sam);

        mxArray *mx_which_bits = mxCreateDoubleMatrix(num_wfs, 1, mxREAL);
        mxSetFieldByNumber(plhs[1], 0, 1, mx_which_bits);
        double *which_bits = mxGetPr(mx_which_bits);

        mxArray *mx_t0 = mxCreateDoubleMatrix(num_wfs, 1, mxREAL);
        mxSetFieldByNumber(plhs[1], 0, 2, mx_t0);
        double *t0 = mxGetPr(mx_t0);

        mxArray *mx_presums = mxCreateDoubleMatrix(num_wfs, 1, mxREAL);
        mxSetFieldByNumber(plhs[1], 0, 3, mx_presums);
        double *presums = mxGetPr(mx_presums);
 
        for (int wf = 0; wf < num_wfs; wf++)
        {
          num_sam[wf] = ((buffer[sub_rec*rec_size + 4*8 + 8*wf + 2] & 0x3F) << 8)
                  + (buffer[sub_rec*rec_size + 4*8 + 8*wf + 3] & 0xFF);
          
          which_bits[wf] = (buffer[sub_rec*rec_size + 4*8 + 8*wf + 4]) & 0X0F;
          
          t0[wf] = (((buffer[sub_rec*rec_size + 4*8 + 8*wf + 5] & 0xFF) << 6)
                  + ((buffer[sub_rec*rec_size + 4*8 + 8*wf + 6] & 0xFC) >> 2)) / fs - 11.0e-6;
          
          presums[wf] = 1 + (buffer[sub_rec*rec_size + 4*8 + 8*wf + 6] & 0x03 << 8)
                  + buffer[sub_rec*rec_size + 4*8 + 8*wf + 7] & 0xFF;
          // mexPrintf("%d of %d %f\n", wf, num_wfs, presums[wf]); // DEBUG
        }
      }
    }
    last_bytes_len = 0;
    
    //mexCallMATLAB(0,NULL,0,NULL, "pause"); // DEBUG
  }
  mxFree(buffer);
  
  // mexPrintf("Read %d records of %d\n", rec, num_rec); // DEBUG
  if (rec != num_rec)
  {
    mxSetM(mx_ver,rec);
    mxSetM(mx_seconds,rec);
    mxSetM(mx_fractions,rec);
    mxSetM(mx_epri,rec);
    mxSetM(mx_offset,rec);
  }
  
  // =======================================================================
  // Go to last incomplete record and copy it
  int last_bytes_out_len = (file_size - ftell(fid));
  mx_last_bytes = mxCreateNumericArray(1, &last_bytes_out_len, mxUINT8_CLASS, mxREAL);
  mxSetFieldByNumber(plhs[2], 0, 0, mx_last_bytes);
  
  unsigned char *last_bytes_out = (unsigned char *)mxGetPr(mx_last_bytes);

  fread(last_bytes_out, 1, last_bytes_out_len, fid);
  
  fclose(fid);
  return;
}

inline unsigned int big2littleint32(unsigned char *buffer)
{
  unsigned int new_val = (buffer[0] << 24)
              + (buffer[1] << 16)
              + (buffer[2] << 8)
              + buffer[3];
  return(new_val);
}


inline int find_next_sync(FILE *fid, int bytes_left)
{
  const int DATA_BLOCK_SIZE = 65536;
  int sync_lock;
  
  // =======================================================================
  // Find next 0xDEADBEEF marker
  unsigned char data[DATA_BLOCK_SIZE];
  int total_read = 0;
  while (total_read < bytes_left)
  {
    int num_to_read;
    if (bytes_left-total_read > DATA_BLOCK_SIZE)
    {
      num_to_read = DATA_BLOCK_SIZE;
    }
    else
    {
      num_to_read = bytes_left-total_read;
    }
    int num_read = fread(data, sizeof(char), num_to_read, fid);
    total_read = total_read + num_read;
    for (int idx = 0; idx < num_read-3; idx++)
    {
      if (data[idx] == 0xDE && data[idx+1] == 0xAD && data[idx+2] == 0xBE && data[idx+3] == 0xEF)
      {
        fseek(fid, -num_read + idx, SEEK_CUR);
        //mexPrintf(" ==> Relocked (%d)\n", ftell(fid));
        sync_lock = 1;
        return(sync_lock);
      }
    }
    fseek(fid, -3, SEEK_CUR);
  }
  
  sync_lock = 0;
  return(sync_lock);
}


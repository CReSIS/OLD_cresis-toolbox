// Compile with mex -v -largeArrayDims arena_packet_strip_mex.cpp
//   -v: verbose
//   -largeArrayDims: required for 64 bit
//   Note: Matlab automatically deallocates memory from mx*alloc
//
// Strips data out of network packet headers from Remote Sensing Solutions
// Arena ADC/DAQ files and stores data records into separate output file.
// Data record byte offset and headers are extracted and returned.
//
// Assumptions:
// 1. Files are <2^31 bytes so the 32 bit integer can be used for record offsets
// 2. Packets are always longer than the record header size
//
// Loading steps:
// 1. Read all packets out of the file and write records into a separate file
// 2. Search for radar sync markers and extract record byte offsets and
//    timing header fields into hdr
// 3. Return hdr which contains file locations for markers and timing header fields
// 4. Return state information for the next call to arena_packet_strip_mex
//
// Run with no arguments to see syntax.
//
// Author: John Paden

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mat.h"
#include "mex.h"

const ptrdiff_t ALLOC_INCREMENT = 128;

/*
 * usage of MEX-file
 */
void
printUsage(int line_number)
{
  mexPrintf("Failed arena_packet_strip_mex on line %d\n", line_number);
  mexPrintf("[hdr,last_bytes_len,num_expected,pkt_counter] \n");
  mexPrintf("   = %s(fn,out_fn,last_bytes,last_bytes_len,num_expected,\n", mexFunctionName());
  mexPrintf("        pkt_counter,min_num_expected,max_num_expected,default_num_expected,\n");
  mexPrintf("        num_header_fields,length_field_offset);\n");
  mexPrintf("\n");
  mexPrintf(" fn: char containing input filename\n");
  mexPrintf(" out_fn: char containing output filename, can be the same as input\n");
  mexPrintf(" last_bytes: uint8 vector (MUST BE PREALLOCATED TO BE num_header_fields*8)\n");
  mexPrintf(" last_bytes_len: int32 containing length of last_bytes, set to -1 on first call\n");
  mexPrintf(" num_expected: int32 containing number of expected bytes before the next sync,\n");
  mexPrintf("   -1 means unknown/unlocked and searching for frame sync, set to -1 on first call\n");
  mexPrintf(" pkt_counter: int32 containing last packet counter value, set to -1 on first call\n");
  mexPrintf(" min_num_expected: minimum record size in bytes. Set to avoid bad header length field causing large data loss.\n");
  mexPrintf(" max_num_expected: maximum record size in bytes. Set to avoid bad header length field causing large data loss.\n");
  mexPrintf(" default_num_expected: default record size in bytes, used when min/max record size is violated\n");
  mexPrintf("   Usually set to min_num_expected.\n");
  mexPrintf(" num_header_fields: int32 containing the number of 8 byte blocks in the header\n");
  mexPrintf("   Header must include the record length field to guarantee files are ready properly.\n");
  mexPrintf(" length_field_offset: int32 containing byte offset of length field from sync\n");
  mexPrintf("\n");
  mexPrintf(" hdr: N by Nx uint64 matrix containing header fields for each record. First 4 bytes of header\n");
  mexPrintf("   sync field is overwritten with the offset to the record in the output file.\n");
  mexPrintf(" last_bytes_len: int32 containing length of valid data stored in last_bytes\n");
  mexPrintf(" num_expected: int32 containing number of expected bytes before the next sync, -1 means unlocked and searching for frame sync\n");
  mexPrintf(" pkt_counter: int32 containing last packet counter value\n");
  mexPrintf("\n");
  mexPrintf(" last_bytes is modified in place\n");
  mexPrintf("\n");
  mexPrintf("When loading a sequence of files, the files should be loaded in the order\n");
  mexPrintf("collected. The first call should have last_bytes_len, num_expected, and\n");
  mexPrintf("pkt_counter all set to -1. Subsequent calls should use the returned values\n");
  mexPrintf("of these fields.\n");
  mexPrintf("\n");
  mexPrintf("Example:\n");
  mexPrintf(" See run_arena_packet_strip.m\n");
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

  if (nlhs!=4) {
    printUsage(__LINE__);
    return;
  }

  /*
   * error checking for input arguments
   */
  if (nrhs!=11) {
    printUsage(__LINE__);
    return;
  }

  if ( !mxIsChar(prhs[0]) )
  {
    mexPrintf("%d: First argument must be input filename\n", __LINE__);
    printUsage(__LINE__);
    return;
  }
  if ( !mxIsChar(prhs[1]) )
  {
    mexPrintf("%d: Second argument must be output filename\n", __LINE__);
    printUsage(__LINE__);
    return;
  }
  
  if ( !mxIsClass(prhs[2], "uint8") )
  {
    mexPrintf("%d: Third argument must be a preallocated uint8 vector\n", __LINE__);
    printUsage(__LINE__);
    return;
  }
  
  for (int input_idx = 3; input_idx <= 10; input_idx++) 
  {
    if ( !mxIsClass(prhs[input_idx], "int32") 
            || mxGetM(prhs[input_idx]) != 1 
            || mxGetN(prhs[input_idx]) != 1)
    {
      mexPrintf("%d: Argument %d must be int32 scalar\n", __LINE__,input_idx+1);
      printUsage(__LINE__);
      return;
    }
  }
  
  // Declare inputs
  char *last_bytes = (char*)mxGetPr(prhs[2]);
  int *last_bytes_len = (int*)mxGetPr(prhs[3]);
  int *num_expected = (int*)mxGetPr(prhs[4]);
  int *pkt_counter = (int*)mxGetPr(prhs[5]);
  int *min_num_expected = (int*)mxGetPr(prhs[6]);
  int *max_num_expected = (int*)mxGetPr(prhs[7]);
  int *default_num_expected = (int*)mxGetPr(prhs[8]);
  int *num_header_fields = (int*)mxGetPr(prhs[9]);
  int *length_field_offset = (int*)mxGetPr(prhs[10]);
  
  ptrdiff_t header_size = *num_header_fields * 8;
  
  // Create variable outputs (allocating the memory later)
  plhs[0]=mxCreateNumericMatrix(0,0,mxUINT64_CLASS,mxREAL);
  plhs[1]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  plhs[2]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  plhs[3]=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);

  /*
   * get filename to read from
   */
  buffersize=mxGetM(prhs[0])*mxGetN(prhs[0])+1;
  filename=(char*)mxCalloc(buffersize,sizeof(char));
  if (filename==NULL)
  {
    mexPrintf("%d: Memory allocation failure\n", __LINE__);
    return;
  }
  mxGetString(prhs[0],filename,buffersize);

  /*
   * read in data file
   */
  FILE *fptr;
  size_t file_size;
  char *data;
  
  fptr = fopen(filename,"rb"); // "b" required for windows to avoid \r issues
  mxFree(filename);
  if (fptr == 0)
  {
    mexPrintf("%d: Input file failed to open\n", __LINE__);
    return;
  }

  fseek(fptr, 0, SEEK_END);
  file_size = ftell(fptr);
  fseek(fptr, 0, SEEK_SET);
  
  data = (char *)mxMalloc(file_size);
  if (data==NULL)
  {
    mexPrintf("%d: Memory allocation failure\n", __LINE__);
    return;
  }
  fread(data, 1, file_size, fptr);

  fclose(fptr);
  
  /*
   * Get output filename to write to
   */
  buffersize=mxGetM(prhs[1])*mxGetN(prhs[1])+1;
  filename=(char*)mxCalloc(buffersize,sizeof(char));
  if (filename==NULL)
  {
    mexPrintf("%d: Memory allocation failure\n", __LINE__);
    return;
  }
  mxGetString(prhs[1],filename,buffersize);
  
  /*
   * Open output file
   */
  FILE *out_fptr;
  
  out_fptr = fopen(filename,"wb"); // "b" required for windows to avoid \r issues
  mxFree(filename);
  if (out_fptr == 0)
  {
    mexPrintf("%d: Output file failed to open\n", __LINE__);
    return;
  }
  
  /*
   * Preallocate hdr array:
   *   This is dynamically reallocated as required during parsing.
   */
  char *hdr;
  ptrdiff_t num_hdr = 0;
  ptrdiff_t alloc_hdr = ALLOC_INCREMENT;
  hdr = (char *)mxMalloc(header_size * alloc_hdr);
  if (hdr==NULL)
  {
    mexPrintf("%d: Memory allocation failure\n", __LINE__);
    fclose(out_fptr);
    return;
  }

  /*
   * Parse data and write out packets into header
   */
  
  // Setup for parsing
  int num_expected_original = -1;
  int last_record_offset = -1;
  int out_offset = 0; // Offset in the output file
  
  ptrdiff_t idx = 0;
  ptrdiff_t out_idx = 0;
  while (idx < file_size)
  {
    // Read in type and length of payload
    int payload_length, new_pkt_counter;
    idx += 12;
    payload_length = ((int*)(data + idx))[0] - 16;
    idx += 16;
    new_pkt_counter = ((int*)(data + idx))[0];
    idx += 4;
    //mexPrintf("%d: %lld %u %d %u\n", __LINE__, idx, payload_length, *pkt_counter, new_pkt_counter);
    
    if (*pkt_counter == -1)
    {
      // First time calling
      mexPrintf("%d: Found sync (%d)\n", __LINE__, *pkt_counter, new_pkt_counter);
      *num_expected = -1;
    }
    else if (new_pkt_counter != (*pkt_counter)+1)
    {
      // Packet dropped so switch to sync lost state
      mexPrintf("%d: Dropped packet last-pkt(%d) new-pkt(%d)\n", __LINE__, *pkt_counter, new_pkt_counter);
      *num_expected = -1;
    }
    *pkt_counter = new_pkt_counter;
    
    // Search for sync and record headers if found
    char *start = (char*)(data+idx); // Start of current packet's payload
    ptrdiff_t offset; // Offset in current packet from start (i.e. the offset into the payload)
    
    // If last_bytes_len > 0, then a sync has already been found in the
    // previous packet (could be from the previous file) and the rest of
    // the header should be right at the start of this packet.
    if (*last_bytes_len > 0)
    {
      // 1. Copy new header to hdr
      if (num_hdr >= alloc_hdr)
      {
        alloc_hdr += ALLOC_INCREMENT;
        char *hdr_tmp;
        hdr_tmp = (char *)mxRealloc(hdr, header_size * alloc_hdr);
        if (hdr_tmp==NULL)
        {
          mexPrintf("%d: Memory allocation failure\n", __LINE__);
          fclose(out_fptr);
          return;
        }
        hdr = hdr_tmp;
      }
      memcpy(hdr+num_hdr*header_size,last_bytes,*last_bytes_len);
      memcpy(hdr+num_hdr*header_size+*last_bytes_len,start,header_size-*last_bytes_len);
      //mexPrintf("%d: %lld %d 0x%016llx\n", __LINE__, idx, num_hdr+1, ((unsigned long long int *)(hdr+num_hdr*header_size))[0]);
      // 2. Write output file offset for this record into hdr (overwriting the first 32 bits of the current header)
      //    offset can be negative and this means that the record starts in the previous file and the offset is measured from the end of the file
      ((unsigned int *)(hdr+num_hdr*header_size))[0] = out_offset - *last_bytes_len;
      num_hdr++;
      // 3a. Check length_field_offset
      ptrdiff_t radar_header_length_offset = 12-*last_bytes_len; // Offset of radar header length field from start
      int expected_length_field_offset;
      if (radar_header_length_offset < 0)
      {
        // Profile data format and length fields were in the last bytes from the previous packet
        expected_length_field_offset = 20+(*(int *)(last_bytes+12));
      }
      else
      {
        // Profile data format and length fields are in this packet
        expected_length_field_offset = 20+(*(int *)(start+radar_header_length_offset));
      }
      if (*length_field_offset != expected_length_field_offset)
      {
        mexPrintf("%d: length_field_offset may be wrong:\n  length_field_offset=%d, expected_length_field_offset=%d\n", __LINE__, *length_field_offset, expected_length_field_offset); // PADEN
      }
      // 3b. Determine number of bytes expected in this record
      offset = *length_field_offset-*last_bytes_len;
      *last_bytes_len = 0;
      int profile_data_format;
      if (offset < 0)
      {
        // Profile data format and length fields were in the last bytes from the previous packet
        profile_data_format = (*(int *)(last_bytes+*length_field_offset-4));
        *num_expected = (*(int *)(last_bytes+*length_field_offset));
        offset = 0;
      }
      else
      {
        // Profile data format and length fields are in this packet
        profile_data_format = (*(int *)(start+offset-4));
        *num_expected = (*(int *)(start+offset));
        offset += 4;
      }
      num_expected_original = *num_expected;
      //mexPrintf("%d: %d\n", __LINE__, *num_expected);
      int num_expected_bins = *num_expected;
      switch (profile_data_format)
      {
        case 0x00000000:
          // num_expected is in units of bytes
          // num_expected_bins needs to be adjusted for 16 bit range bins
          //  >>1 = /2, 2 byte samples
          num_expected_bins = num_expected_bins >> 1;
          break;
        case 0x00010000:
          // num_expected is in units of bytes
          // num_expected_bins needs to be adjusted for 16 bit IQ range bins
          //  >>2 = /4, 4 byte samples
          num_expected_bins = num_expected_bins >> 2;
          //*num_expected = *num_expected >> 1; // These three lines replaced the above line for temporary hack;
          //*num_expected = *num_expected - 1536*4; // con't from previous line.
          //num_expected_bins = num_expected_bins >> 3; // con't from previous line.
          break;
        case 0x00020000:
        case 0x00030000:
          // num_expected is in units of bytes
          // num_expected_bins needs to be adjusted for 8 byte bins:
          //   >>3 = /8, 2 IQ channels, 4 byte samples or 2*4 = 8
          num_expected_bins = num_expected_bins >> 3;
          //*num_expected = *num_expected << 3; // This line replaced the above line for temporary hack
          break;
      }
      if (num_expected_bins < *min_num_expected || num_expected_bins > *max_num_expected)
      {
        mexPrintf("%d: %d Bad Record Length %d\n", __LINE__, ((unsigned int *)(hdr+num_hdr*header_size))[0], *num_expected);
        *num_expected = *default_num_expected;
      }
    }
    else
    {
      offset = 0;
    }

    // Start parsing the new packet data
    while (offset < payload_length)
    {
      if (*num_expected >= payload_length-offset)
      {
        // Skip to end of packet since num_expected bytes in record is more than the packet length
        *num_expected -= (payload_length-offset);
        offset = payload_length;
        //mexPrintf("%d: %d\n", __LINE__, *num_expected);
      }
      else
      {
        if (*num_expected >= 0)
        {
          // Skip to where the sync should be
          offset += *num_expected;
          *num_expected = 0;
        }
        // Search for sync (if sync locked, it should be the first uint64)
        while (offset < payload_length)
        {
          const long long c_sync = 0x7F80000080000000;
          if (((unsigned long long *)(start+offset))[0] == c_sync)
          {
            if (*num_expected == -1)
            {
              mexPrintf("%d: Sync Found At %d\n", 
                      __LINE__, idx+offset);
            }
            last_record_offset = idx+offset;
            // Found sync
            if (offset + header_size > payload_length)
            {
              // Header is incomplete
              // 1. Fill last_bytes
              memcpy(last_bytes,start+offset,payload_length-offset);
              *last_bytes_len = payload_length-offset;
              // 2. Skip to end
              offset = payload_length;
              break;
            }
            else
            {
              // Header is complete
              // 1. Copy header to hdr
              if (num_hdr >= alloc_hdr)
              {
                alloc_hdr += ALLOC_INCREMENT;
                char *hdr_tmp;
                hdr_tmp = (char *)mxRealloc(hdr, header_size * alloc_hdr);
                if (hdr_tmp==NULL)
                {
                  mexPrintf("%d: Memory allocation failure\n", __LINE__);
                  fclose(out_fptr);
                  return;
                }
                hdr = hdr_tmp;
              }
              memcpy(hdr+num_hdr*header_size,start+offset,header_size);
              //mexPrintf("%d: %lld %d 0x%016llx\n", __LINE__, idx+offset, num_hdr+1, ((unsigned long long int *)(hdr+num_hdr*header_size))[0]);
              // 2. Write output file offset into header (overwriting the first 32 bits of the header)
              ((unsigned int *)(hdr+num_hdr*header_size))[0] = out_offset + offset;
              num_hdr++;
              // 3a. Check length_field_offset
              ptrdiff_t radar_header_length_offset = offset + 12; // Offset of radar header length field from start
              int expected_length_field_offset;
              expected_length_field_offset = 20+(*(int *)(start+radar_header_length_offset));
              if (*length_field_offset != expected_length_field_offset)
              {
                mexPrintf("%d: length_field_offset may be wrong:\n  length_field_offset=%d, expected_length_field_offset=%d\n", __LINE__, *length_field_offset, expected_length_field_offset); // PADEN
              }
              // 3b. Determine number of bytes expected in this record
              offset += *length_field_offset;
              *last_bytes_len = 0;
              int profile_data_format;
              profile_data_format = (*(int *)(start+offset-4));
              *num_expected = (*(int *)(start+offset));
              num_expected_original = *num_expected;
              offset += 4;
              //mexPrintf("%d: %d\n", __LINE__, *num_expected);
              int num_expected_bins = *num_expected;
              switch (profile_data_format)
              {
                case 0x00000000:
                  // num_expected is in units of bytes
                  // num_expected_bins needs to be adjusted for 16 bit range bins
                  //  >>1 = /2, 2 byte samples
                  num_expected_bins = num_expected_bins >> 1;
                  break;
                case 0x00010000:
                  // num_expected is in units of bytes
                  // num_expected_bins needs to be adjusted for 16 bit IQ range bins
                  //  >>2 = /4, 4 byte samples
                  num_expected_bins = num_expected_bins >> 2;
                  //*num_expected = *num_expected >> 1; // These three lines replaced the above line for temporary hack;
                  //*num_expected = *num_expected - 1536*4; // con't from previous line.
                  //num_expected_bins = num_expected_bins >> 3; // con't from previous line.
                  break;
                case 0x00020000:
                case 0x00030000:
                  // num_expected is in units of bytes
                  // num_expected_bins needs to be adjusted for 8 byte bins:
                  //   >>3 = /8, 2 IQ channels, 4 byte samples or 2*4 = 8
                  num_expected_bins = num_expected_bins >> 3;
                  //*num_expected = *num_expected << 3; // This line replaced the above line for temporary hack
                  break;
              }
              if (num_expected_bins < *min_num_expected || num_expected_bins > *max_num_expected)
              {
                mexPrintf("%d: %d Bad Record Length %d\n", __LINE__, ((unsigned int *)(hdr+num_hdr*header_size))[0], num_expected_bins);
                *num_expected = *default_num_expected;
              }
              break;
            }
          }
          else
          {
            if (*num_expected == 0)
            {
              // Since num_expected != -1, we were locked and have now become unlocked
              mexPrintf("%d: Sync Lost LastSync(%d) ExpectedSync(%d) LastRecordSize(%d)\n",
                      __LINE__, last_record_offset, idx+offset, num_expected_original);
//               mexPrintf("  %d %d %d\n", payload_length, out_offset+offset,((unsigned int *)(hdr+(num_hdr-1)*header_size))[0]);
//               for (int debug_idx=0; debug_idx<1; debug_idx++)
//               {
//                 mexPrintf("  0x%016llx\n", ((unsigned long long *)(start+offset-64))[debug_idx]);
//               }
              if (num_hdr > 0)
              {
                // Print out last header to help with debugging
                for (int header_idx=0; header_idx<*num_header_fields; header_idx++)
                {
                  mexPrintf("  %d 0x%016llx\n", header_idx, ((unsigned long long int *)(hdr+(num_hdr-1)*header_size))[header_idx]);
                }
              }
              else
              {
                mexPrintf("  First record so no record header to print\n");
              }
              *num_expected = -1;
            }
            offset += 8;
          }
        }
      }
    }
    // Write out new packet
    size_t num_bytes;
    num_bytes = fwrite(data+idx, 1, payload_length, out_fptr);
    if (num_bytes != payload_length)
    {
      mexPrintf("%d: File write failure\n", __LINE__);
      fclose(out_fptr);
      return;
    }
    out_offset += payload_length;
    idx += payload_length;
  }
  fclose(out_fptr);
  
  // Copy hdr to output
  /* Point mxArray to dynamicData */
  mxSetData(plhs[0], hdr);
  mxSetM(plhs[0], *num_header_fields);
  mxSetN(plhs[0], num_hdr);
  
  ((int *)mxGetPr(plhs[1]))[0] = *last_bytes_len;
  ((int *)mxGetPr(plhs[2]))[0] = *num_expected;
  ((int *)mxGetPr(plhs[3]))[0] = *pkt_counter;

  return;
}

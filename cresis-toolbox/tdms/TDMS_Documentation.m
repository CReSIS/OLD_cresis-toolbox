%                           TDMS_Documentation
%
%   INTRO:
%   =======================================================================
%       The main file for reading TDMS files is called TDMS_readTMDSfile.
%
%   The goal of this function was to simply serve in reading the data from
%   the TDMS file.  I have explicitly chosen to not include any batch
%   processing, file selection, or post processing (except for some minimal
%   things) into this file.  I hope this simplifies the file, and that any
%   changing of the way the data is packaged, saving, batch processing,
%   etc, can easily be written by a wrapper function.
%
%   FUNCTIONS USER CAN CALL/USE:
%   =======================================================================
%   see list in TDMS_readTDMSfile.m
%
%   OPTIONAL PARAMTERS:
%   =======================================================================
%       There are a few set of optional parameters which are defined fully
%   in the TDMS_readTDMSfile.m & TDMS_retrievingSubsets.m
%
%   ACTIVE ISSUES:
%   =======================================================================
%   See TDMS_TODO_List
%
%   IMPORTANT USAGE NOTES
%   =======================================================================
%   1) Properties are not cast to double, this can cause unexpected results 
%   when performing math with mixed data types i.e. you can get 3/4 => 1
%
%   2) Only the latest property value is kept, one could make an effort to
%   read the property value as it updates, but since one can not ensure the
%   size of the write segments, i.e. 1 write call ~= 1 segment, then the 
%   time at which the property changed is a bit ambiguous (see
%   TDMS_fileFormatNotes.m)
%
%   3) Timestamps have a rediculous setup for accuracy.  In some cases this
%   accuracy may be real, but I have decided to simplify this result down
%   to a double, thus if you really need the accuracy that the timestamp
%   allows for, you may be out of luck with this code
%
%   4) Timestamps are in UTC format.  They can be converted to local time by
%   using of the optional parameters (see above).
%
%   5) Timestamp data is kept in double format.  It is processed so that
%   the function datestr() will convert the double value to a human
%   readable form.

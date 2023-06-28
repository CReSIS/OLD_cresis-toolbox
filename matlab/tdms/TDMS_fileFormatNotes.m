%TDMS_fileFormatNotes

%1) TDMS is a binary format for quickly dumping data to disk, along with
%additional properties

%2) There is no compression in a TDMS file

%3) TDMS files consist of meta data sections, and raw data sections. In
%general the bulk of the processing comes from interpreting the meta data
%sections.

%4) To reduce read time, you can defragment the TDMS file, which rewrites
%the file so as to minimize the # of meta data sections

%5) TDMS does not currently support 2 dimensional arrays.  To save a 2d
%array the array first needs to be linearized and then saved.  Upon
%reading, the data can be reshaped to the original 2d matrix size.

%6) National Instruments does not specify how chunks of data relate to
%write commands.  In other words, each write does not necessarily correspond 
%to a data chunk.

%7) You can change the property value in a TDMS file, once it has
%been written, but only the last property value is considered valid.  In
%other words, let's say with every write, you save a boolean as to whether that
%data is good or bad, in an associated write.  Thus, for every data write,
%you have a property write.  First, there is no guarantee that each
%changing of the property value will actually show up.  Second, there is no
%guarantee that you will have the same # of property value changes as you
%will # of data sections, to later match things up.  Thus, although it
%might be possible to retain the property at different points in the file,
%it most likely will not make sense, and thus only the last value is
%returned.

%8) In general, the more you write to TDMS files in Labview, the more
%fragemented the data will become (more meta sections).  This is especially
%true when the type of data is changing between writes such as writing a
%set of integers, then a set of doubles, then a set of integers again. This
%can be fixed later by defragmenting but can also be fixed by hold more
%values in memory and writing less often.

%9) TDMS only has 3 different types of objects, a root, which holds
%everything, a group object, which holds channel objects, and channel
%objects, which hold data.  TDMS does not support storing raw data in the
%group or root objects.  All channel objects must belong to a group object.

%10 ASCII Illustration of file structure
%
%   SEGMENT 1            SEGMENT 2    SEGMENT 3
%---------------------- ------------ --------------------------------
%META   RAW             META RAW     META    RAW
%mmmmm rrrrrrrrrrrrrrrr mmmm rrrrrrr mmmmmmm rrrrrrrrrrrrrrrrrrrrrrrrrrrr
%
%Let's look at "RAW" for segment 3
%
%       Chunk 1              Chunk 2
% -------------------  -------------------
%Chan 1 Chan2   Chan3 Chan 1 Chan2   Chan3
%111111122222222333333111111222222222333333
%
%   NOTE: The # of values for each channel can not be changed for each chunk





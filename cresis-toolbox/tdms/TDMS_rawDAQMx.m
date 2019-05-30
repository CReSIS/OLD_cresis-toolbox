%{

TDMS files can store many different types of data inside of them, such as 
uint8, int16, doubles, etc. A full list is outside the scope of this document.

One of these data types is rawDAQmx. I think this data type is actually the
bit representation captured by the DAQ (hence the name raw). The presence
of this data type in the meta for a TDMS file refers to the fact that the
data has a true IEEE like data type underneath, as well as properties that
are stored with the file that inform the user how to scale the data from
the the bit-level representation to a true number.

The primary reason one might use this format is because it gets rid of the
intermediate step of converting the data to a known data type before
saving. This could be advantageous because it could save computational
time and memory from not doing the conversion, as well as time and disk space
from writing a smaller file to disk.

This data type however is not an open specification. In the TDMS file the
header before the data will specify that upcoming data has this data type,
allowing one to skip over that data section.

I have contacted the NI about this and received the following response:
%==========================================================================
I apologize but since most of the TDMS format  is a National Instrument
proprietary file structure I can really only help you out with any of the
information that is already provided.  I have included a link to a
KnowledgeBase below that walks through TDMS in a lot more detail.  Please
have a look through and hopefully it will answer some of your questions.

TDMS File Format Internal Structure
http://zone.ni.com/devzone/cda/tut/p/id/5696
%==========================================================================

That being said I have been able to obtain the data from this format in
some example files that I created. However, since these example files are
relatively few, and since I have no idea how to generalize my code, I have
not formally included support for this data type. If I get some more
example files with rawDAQmx data I would likely include a mode in which one
could attempt to retrieve data from this data type.

Finally, National Instruments does provide some old drivers. I think they
only work for 32 bit Windows but they may support the rawDAQmx data type
...

%}
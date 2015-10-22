Install Visual C++

unzip json-c-master.zip

1) Open json-c.vcproj with Visual C++

2) Modify json-util.c file as following:
.
.
.

line 159:
static void sscanf_is_broken_test()
{
	int ret_errno;
    int64_t num64;
    int is_int64_min;
    int ret_errno2;
    int is_int64_max;


	(void)sscanf(" -01234567890123456789012345", "%" SCNd64, &num64);
	/*int*/ ret_errno = errno;
	/*int*/ is_int64_min = (num64 == INT64_MIN);

	(void)sscanf(" 01234567890123456789012345", "%" SCNd64, &num64);
	/*int*/ ret_errno2 = errno;
	/*int*/ is_int64_max = (num64 == INT64_MAX);

	if (ret_errno != ERANGE || !is_int64_min ||
	    ret_errno2 != ERANGE || !is_int64_max)
	{
		MC_DEBUG("sscanf_is_broken_test failed, enabling workaround code\n");
		sscanf_is_broken = 1;
	}
}
.
.
.




3) Modify json_objet.c  as following

555	 -  return sprintbuf(pb, "%f", jso->o.c_double);
555	 +  return sprintbuf(pb, "%.15g", jso->o.c_double);



4) Compile the library to: 
H:\scripts\matlab\json_mex\json-c-master\Release


5) Change "...\matlab-json-master\tojson.c(1)":

#include <json/json.h>    ==>   #include <json.h>



6) In Matlab under working directory of "matlab-json-master":

mex -ljson -g LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:MSVSRT.lib /NODEFAULTLIB:LIBCMT.lib" -IH:\scripts\matlab\json_mex\json-c-master\ -LH:\scripts\matlab\json_mex\json-c-master\Release tojson.c

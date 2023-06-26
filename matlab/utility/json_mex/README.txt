
Install:
autoconf, automake, libtool

TO MAKE STATIC JSON-C:

mkdir /users/hwang/bin/

INSIDE /users/hwang/scripts/matlab/json_mex/json-c-master DO THESE TASKS:

Run:
./configure --prefix=/users/hwang/bin/

Then modify the Makefile by adding "-fPIC" to the compiler flags:
CC = gcc
CCDEPMODE = depmode=gcc3
CFLAGS = -g -O2 -fPIC
CPP = gcc -E
CPPFLAGS =  -fPIC

Then run:
make
make install


TO MAKE STATIC fromjson AND tojson:

Move the .so files out of the /users/hwang/bin/lib/ directory so that the
compiler will not see them.

Copy mexopt.sh and add -fPIC to compiler flags

        glnxa64)
#----------------------------------------------------------------------------
            RPATH="-Wl,-rpath-link,$TMW_ROOT/bin/$Arch"
            # StorageVersion: 1.0
            # CkeyName: GNU C
            # CkeyManufacturer: GNU
            # CkeyLanguage: C
            # CkeyVersion:
            CC='gcc'
            CFLAGS='-ansi -D_GNU_SOURCE'
            CFLAGS="$CFLAGS  -fexceptions"
            CFLAGS="$CFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CLIBS="$RPATH $MLIBS -lm"
            COPTIMFLAGS='-O -DNDEBUG'
            CDEBUGFLAGS='-g'
            CLIBS="$CLIBS -lstdc++"
#
            # C++keyName: GNU C++
            # C++keyManufacturer: GNU
            # C++keyLanguage: C++
            # C++keyVersion: 
            CXX='g++'
            CXXFLAGS='-ansi -D_GNU_SOURCE'
            CXXFLAGS="$CXXFLAGS -fPIC -fno-omit-frame-pointer -pthread"
            CXXLIBS="$RPATH $MLIBS -lm"
            CXXOPTIMFLAGS='-O -DNDEBUG'
            CXXDEBUGFLAGS='-g'
...
            LD="$COMPILER"
            LDEXTENSION='.mexa64'
            LDFLAGS="-pthread -shared
-Wl,--version-script,$TMW_ROOT/extern/lib/$Arch/$MAPFILE -Wl,--no-undefined"
            LDOPTIMFLAGS='-O'
            LDDEBUGFLAGS='-g'



mex -v -L/users/hwang/bin/lib/ -f mexopts.sh -ljson-c -I/users/hwang/bin/include/json  -g fromjson.c

mex -v -L/users/hwang/bin/lib/ -f mexopts.sh -ljson-c -I/users/hwang/bin/include/json  -g tojson.c

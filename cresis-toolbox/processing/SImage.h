#ifndef __SIMAGE_H__
#define __SIMAGE_H__

#include "DTwoDimArray.h"
#include <string.h>

// A very simple image class.
//
class SDoublePlane : public _DTwoDimArray<double>
{
 public:
  SDoublePlane() { }
  SDoublePlane(int _rows, int _cols)  : _DTwoDimArray<double>(_rows, _cols)
    { 
      // be nice and initialize plane to all 0's
      memset(data_ptr(), 0, sizeof(double) * rows() * cols());
    }

 SDoublePlane(int _rows, int _cols, double *ptr)  : _DTwoDimArray<double>(_rows, _cols, ptr)
    { 
    }


};

#endif

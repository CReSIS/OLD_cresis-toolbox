// This class is just a basic, templated 2-D array.
//
// David Crandall, 2003-
// crandall@cs.cornell.edu
//

#ifndef __DTWODIMARRAY_H__
#define __DTWODIMARRAY_H__
#include <assert.h>

template<class T>
class _DTwoDimArray
{
 public:

  ///////////////////////////////////////////////////////
  // constructors
  _DTwoDimArray()
    {
      data = 0;
      data_area = 0;
      _rows = _cols = 0;
    }

  _DTwoDimArray(int __rows, int __cols)
  {
    _rows = __rows;
    _cols = __cols;
    
    data = 0;
    data_area = 0;
    
    initialize_storage();
  }

  _DTwoDimArray(int __rows, int __cols, const T *array)
    {
      _rows = __rows;
      _cols = __cols;
      
      data = 0;
      data_area = 0;
      
      initialize_storage();
      
      memcpy(data_area, array, _rows * _cols * sizeof(T));
    }

  _DTwoDimArray(const _DTwoDimArray<T> &other)
    {
      assert(this != &other);
      
      data = 0;
      data_area = 0;
      
      *this = other;
    }


  ///////////////////////////////////////////////////////
  // destructor
  ~_DTwoDimArray()
    {
      deallocate_storage();
    }


  ///////////////////////////////////////////////////////
  // assignment operator
  _DTwoDimArray<T> &operator=(const _DTwoDimArray<T> &other)
    {
      if(this == &other)
	return *this;
      
      // profiler->begin(4);
      if(!data || _rows != other.rows() || _cols != other.cols())
	{
	  _rows = other.rows();
	  _cols = other.cols();
	  
	  initialize_storage();
	}
      
      memcpy(data_area, other.data_area, _rows * _cols * sizeof(T));
      
      // profiler->end(4);
      return *this;
    }


  ///////////////////////////////////////////////////////
  // element access
  inline T *operator[](int row) const { return data[row]; }


  ///////////////////////////////////////////////////////
  // informational 
  inline int rows() const { return _rows; }
  inline int cols() const { return _cols; }
  inline T *data_ptr() const { return data_area; }
  inline T **row_pointers() const { return data; }

 protected:

  void deallocate_storage()
    {
      if(data)
	{
	  delete[] data;
	  delete[] data_area;
	  
	  data = 0;
	  data_area = 0;
	}
    }

  void initialize_storage()
  {
    //  profiler->begin(2);
    if(data)
      deallocate_storage();
    
    if(_rows > 0)
      {
	data = new T *[_rows];
	data_area = new T[_rows * _cols];
      }
    else
      {
	data = 0;
	data_area = 0;
      }
    
    T *cp = data_area;
    for(int i=0; i<_rows; i++, cp+=_cols)
      data[i] = cp;
    //  profiler->end(2);
  }


  /// Array of pointers to the beginning of each matrix row in data_area
  T **data;

  /// Pointer to actual matrix data
  T *data_area;

  /// Size of matrix
  int _rows, _cols;
};


#endif

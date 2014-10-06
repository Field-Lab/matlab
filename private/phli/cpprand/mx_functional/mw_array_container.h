/* mw_array_container.h
 * Bring in some STL container functionality for mwArray.
 * Derived from: http://jeffmatherphotography.com/dispatches/2008/08/using-c-iterators-on-matlab-mxarrays/
 * Ver 0.4
 * Peter H. Li 2011 FreeBSD License
 */

#ifndef MW_ARRAY_CONTAINER
#define MW_ARRAY_CONTAINER

#include <mex.h>
#include "skip_ptr.h"


// Treat an mxArray like a container.
// (See Josuttis. "The C++ Standard Library." 1999. p. 219-220)
template <class T> class MWArrayContainer 
{
public:
  typedef T                          value_type;
  typedef T*                         pointer;
  typedef const T*                   const_pointer;
  typedef skip_ptr<T>                iterator;
  typedef const skip_ptr<T, const T> const_iterator;
  typedef T&                         reference;
  typedef const T&                   const_reference;
  typedef mwIndex                    size_type;
  typedef ptrdiff_t                  difference_type;
  
  typedef std::reverse_iterator<iterator>       reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  MWArrayContainer() {}
  
  MWArrayContainer(const mxArray *inarr) {
    _data = static_cast<pointer>(mxGetData(inarr));
    _nrows = mxGetM(inarr);
    _ncols = mxGetN(inarr);
    _nelem = mxGetNumberOfElements(inarr);
    _ndim  = mxGetNumberOfDimensions(inarr);
    _dims = mxGetDimensions(inarr);
    _step = 1;
  }
  
  
  // Derived structs
  MWArrayContainer get_col(mwIndex col) const {
    const size_type dims[] = {_nrows, 1};
    return MWArrayContainer(begin() + col*_nrows, 2, dims, _step);
  }
  
  MWArrayContainer get_row(mwIndex row) const {
    const size_type dims[] = {1, _ncols};
    return MWArrayContainer(begin() + row, 2, dims, _nrows*_step);
  }
  
  MWArrayContainer get_block(int blocknum, int nblocks) const {
    difference_type blocksize = _nelem / nblocks;
    const_iterator blockstart = begin() + blocknum*blocksize;
  
    size_type dims[2];
    if (blocknum == (nblocks-1)) dims[0] = end() - blockstart;
    else dims[0] = blocksize;
    dims[1] = 1;
    
    return MWArrayContainer(blockstart, 2, dims, _step);
  }
  
  
  // Iterator support
  iterator       begin()       { return iterator(_data, _step); }
  const_iterator begin() const { return const_iterator(_data, _step); }
  iterator       end()       { return begin() + _nelem; }
  const_iterator end() const { return begin() + _nelem; }
  
  // Direct element access
  reference       operator[](size_t ind)       { return *(begin() + ind); }
  const_reference operator[](size_t ind) const { return *(begin() + ind); }
  
  // Size is constant
  size_type size() const { return _nelem; }
  size_type max_size() const { return _nelem; }
  
  
  
  // Additional iterators
  // Need to think hard about steps, etc. here
  iterator rowbegin(size_type row) { return iterator(_data + row*_step, _nrows*_step); }
  iterator rowend(size_type row) { return rowbegin(row) + _ncols; }
  iterator colbegin(size_type col) { return iterator(_data + col*_nrows*_step, _step); }
  iterator colend(size_type col) { return colbegin(col) + _nrows; }
  
  
  
  // Getters
  size_type nrows() const { return _nrows; }
  size_type ncols() const { return _ncols; }
  size_type nelem() const { return _nelem; }
  size_type ndim()  const { return _ndim;  }
  const size_type* const dims() const { return _dims; } // Needs testing
  
  
private:
  MWArrayContainer(const_iterator begin, size_type ndim, const size_type *dims, difference_type step) :
    _data(begin.get()), _ndim(ndim), _dims(dims), _step(step), _nrows(dims[0]), _ncols(dims[1]), _nelem(1) {
      for (int i = 0; i < ndim; i++) _nelem *= dims[i];
  }
    
  pointer _data;
  size_type _nrows, _ncols, _nelem, _ndim; 
  const size_type *_dims;
  difference_type _step;
};


#endif
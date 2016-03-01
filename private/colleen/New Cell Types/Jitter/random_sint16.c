#include <matrix.h>
#include <mex.h>   
#include "Random.h"

#define NDIMS 2


SInt16 RandJavaShort(RandJavaState state) {
  *state = (*state * 0x5DEECE66DLL + 0xBLL) & 0xFFFFFFFFFFFFLL;
  return (SInt16) (*state >> 32LL);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{

    SInt64 *state;
    SInt16 *draw;
    const mwSize dims[]={1,1};
    
    state = (SInt64 *)mxGetData( prhs[0] );
  
    //associate outputs
    plhs[0] = mxCreateNumericArray(NDIMS,dims,mxINT16_CLASS,mxREAL);
    draw = (SInt16 *)mxGetData( plhs[0] );
    
    
    //Actually do something
    *draw = RandJavaShort(state); 

    return;
}
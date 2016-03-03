
#include "math.h"
#include <matrix.h>
#include <mex.h>
#include "Random.h"
#include "stdio.h"
#include <stdlib.h>
#include <string.h>

#define NDIMS_A 2
#define NDIMS_B 3
#define RANDJAVA_SCALE (4294967296.0)
//#define RANDJAVA_SCALE (65535.0f)

SInt16 RandJavaNbit(RandJavaState state) {
    *state = (*state * 0x5DEECE66DLL + 0xBLL) & 0xFFFFFFFFFFFFLL;

        return state;

}

SInt32 RandJavaLong(RandJavaState state) {
    *state = (*state * 0x5DEECE66DLL + 0xBLL) & 0xFFFFFFFFFFFFLL;
    return (SInt32) (*state >> 16LL);
}

float RandJavaFloat(RandJavaState state)
{
    return(((float)(UInt16)(RandJavaLong(state)>>16LL))/65535.0f);
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    SInt16 draw;
    SInt64 *state;
    
    mxArray *lut_in_mxarray;
    unsigned char *lut;
            
    mxArray *map_in_mxarray;
    unsigned short *map;
    
    mxArray *backrgb_in_mxarray;
    unsigned char *backrgb;
    
    unsigned char *image_pattern, *prefilled_seq;
    int w, h, width, height, image_index, map_index, lut_index, noise_type, n_bits, cnt, map_value,  m_width, m_height;
    float probability;
    
    //get seed
    state = (SInt64 *)mxGetData( prhs[0] );
    
    // array size
    width = mxGetScalar(prhs[1]);
    height = mxGetScalar(prhs[2]);
    
    // lut
    lut_in_mxarray = mxDuplicateArray(prhs[3]);
    lut = (unsigned char *)mxGetData( lut_in_mxarray );
    
    // map
    map_in_mxarray = mxDuplicateArray(prhs[4]);
    map = (unsigned short *)mxGetData( map_in_mxarray );
    
    // background [R G B] for map and sparse
    backrgb_in_mxarray = mxDuplicateArray(prhs[5]);
    backrgb = (unsigned char *)mxGetData( backrgb_in_mxarray );
    
    // map size
    m_width = mxGetScalar(prhs[6]);
    m_height = mxGetScalar(prhs[7]);
    
    //RGB/BW, Gaussian/Binary, sparse probability flags
    noise_type = mxGetScalar(prhs[8]);
    n_bits = mxGetScalar(prhs[9]);
    probability = mxGetScalar(prhs[10]);

    //associate output
    const mwSize dims_B[]= {4,width,height};
    plhs[0] = mxCreateNumericArray(NDIMS_B,dims_B,mxUINT8_CLASS,mxREAL);
    draw = (unsigned char *)mxGetData( plhs[0] );
    
    
    int if_fill = 1;
    
    for( h=0; h < m_height; h++) {
        draw = RandJavaNbit( state ));
    };
        

    return;
}  // mex file

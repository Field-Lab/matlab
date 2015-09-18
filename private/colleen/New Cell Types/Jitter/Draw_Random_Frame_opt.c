
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

SInt16 RandJavaNbit(RandJavaState state, int n_bit) {
    *state = (*state * 0x5DEECE66DLL + 0xBLL) & 0xFFFFFFFFFFFFLL;
    
    if (n_bit == 1){
        return (SInt16) (*state >> 47LL);}
    else if (n_bit == 3){
        return (SInt16) (*state >> 45LL);}
    else if (n_bit == 8){
        return (SInt16) (*state >> 40LL); }
    else{
        return (SInt16) (*state >> 32LL); }
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
    unsigned char *lut
            
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
    image_pattern = (unsigned char *)mxGetData( plhs[0] );
    
    int map_lined[4][m_width][m_height];
    prefilled_seq = (unsigned char *)(map_lined);
    
    int if_fill = 1;
    
    for( h=0; h < m_height; h++) {
        
        image_index = 4 * h * m_width;       
        
        for( w=m_width; w != 0 ; w--) {
                        
            if (probability != 1.0)
                if_fill =  RandJavaFloat(state) < probability;
            
            if (if_fill){ // fill color values                
                for (cnt = 0; cnt<3; cnt++){
                    if (noise_type==3 || cnt==0) //Gauss RGB - 2 additional draws
                        lut_index = (int)(RandJavaNbit( state, n_bits )*3);
                    prefilled_seq[image_index++] = lut[lut_index + cnt];
                }
            }
            else{ // fill background values
                for (cnt=0; cnt<3; cnt++){
                    prefilled_seq[image_index++] = backrgb[cnt];
                }
            }
            prefilled_seq[image_index++] = (unsigned char)0xff; // alpha channel
            
        }
    }
    
    if (m_height==1){ // map based noise, transform. Make sure in matlab class that height is 1, width is not!
        
        for( h=0; h < height; h++) {
            image_index = 4 * h * width;
            map_index = h * width;
            for( w=0; w < width; w++) {
                
                map_value = (int) map[map_index++];
                
                if (map_value>0){ //'cone'
                    cnt = map_value * 4;
                    image_pattern[image_index++] = prefilled_seq[cnt++];
                    image_pattern[image_index++] = prefilled_seq[cnt++];
                    image_pattern[image_index++] = prefilled_seq[cnt];
                }
                else{ // intercone space
                    image_pattern[image_index++] = backrgb[0];  //  R
                    image_pattern[image_index++] = backrgb[1];  // G
                    image_pattern[image_index++] = backrgb[2];  // B
                }
                image_pattern[image_index++] = (unsigned char)0xff;
            }
            
        }
    }
    
    else{ // simply copy stuff
        memcpy(image_pattern, prefilled_seq, sizeof(unsigned char)*4*width*height);
    }

    return;
}  // mex file

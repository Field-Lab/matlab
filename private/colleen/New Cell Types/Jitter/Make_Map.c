#include <matrix.h>
#include <mex.h>   
#include "stdio.h"

#define NDIMS_A 2
#define NDIMS_B 3

// inputs
// width, height, lut in mxarray, map_in_mxarray, backrgb vect (3 element) 



// No checking... should there be safety checking? like for lut_index beyo1nd bound?what about proper num inputs?

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{
    
    unsigned char *image_mat; // unsigned char is uint8
    
    mxArray *lut_in_mxarray;
    unsigned char *lut;
    
    mxArray *map_in_mxarray;
    unsigned char *map;
    
    mxArray *backrgb_in_mxarray;
    unsigned char *backrgb;

    
    int w, h;
    //int width, height, r_index, g_index, b_index, lut_index;
    int width, height, map_index, image_index, lut_index, map_value;
    
    // Retrieve input scalar size values
    width = mxGetScalar(prhs[0]);
    height = mxGetScalar(prhs[1]);
    
    // Retrieve lut 
    lut_in_mxarray = mxDuplicateArray(prhs[2]);  
    lut = (unsigned char *)mxGetData( lut_in_mxarray );
    
    // Retrieve map
    map_in_mxarray = mxDuplicateArray(prhs[3]);  
    map = (unsigned char *)mxGetData( map_in_mxarray );

    // Retrieve back rgb vect
    backrgb_in_mxarray = mxDuplicateArray(prhs[4]);  
    backrgb = (unsigned char *)mxGetData( backrgb_in_mxarray );
   
    //associate outputs
    // dummy the size of the mat and try to pre-fill with dummy alpha channel values
    const mwSize dims_B[] = {4,width,height};
    plhs[0] = mxCreateNumericArray(NDIMS_B,dims_B,mxUINT8_CLASS,mxREAL);
    
    image_mat = (unsigned char *)mxGetData( plhs[0] );

    // from random matrix fill
        for( h=0; h<height; h++) {  
        
            for( w=0; w<width; w++) {
                    
                
                image_index = (h * (4 * width)) + (4 * w);
                map_index = (h * width) + w;
                
                // The map uses the same mat index as the output map
                map_value = (int)map[ map_index ];
                
                
                // Fill out first three color entries in relevant matrix index 
                if ( map_value < 1 ) {
     
                    image_mat[ image_index ] = backrgb[0];  //  R
                
                    image_mat[ image_index + 1 ] = backrgb[1];  // G
                 
                    image_mat[ image_index + 2 ] = backrgb[2];  // B
                }
                else {  
                    
                    
                     // lut was a single vector (r1, g1, b1, r2, g2...
                    // The lut_index needs to be expanded so that it can be used properly in the linearized lut table
                    // So if map_index = 1 (for first lookup level above background) then the first lut index is 0.
                    
                    lut_index = (int)(3 * (map_value - 1)); 

                    image_mat[ image_index ] = lut[ lut_index ];  //  R
                    
                    image_mat[ image_index + 1 ] = lut[ lut_index + 1 ];  // G
                 
                    image_mat[ image_index + 2 ] = lut[ lut_index + 2];  // B
                }
                 
                image_mat[ image_index + 3 ] = (unsigned char)0xff;     // Fill the 4th color (alpha) channel with 255 as a filler
                
            }       // height
        }           // width
    

    return;
}

#include "math.h"
#include <matrix.h>
#include <mex.h>
#include "Random.h"
#include "stdio.h"
#include <stdlib.h>
#include <string.h>

#define NDIMS_A 2
#define NDIMS_B 3
#define NDIMS_C 3

SInt16 RandJavaNbit(RandJavaState state, int n_bit) {
    *state = (*state * 0x5DEECE66DLL + 0xBLL) & 0xFFFFFFFFFFFFLL;
    
    if (n_bit == 1){
        return (SInt16) (*state >> 47LL);}
    else if (n_bit == 3){
        return (SInt16) (*state >> 45LL);}
    else{
        return (SInt16) (*state >> 40LL); }
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
    
    mxArray *spikes_array;
    
    mxArray *map_in_mxarray;
    unsigned short *map;
    
    mxArray *backrgb_in_mxarray;
    
    unsigned char *lut, *backrgb, *prefilled_seq, *spikes;
    int w, h, width, height, image_index, map_index, lut_index, noise_type, n_bits, cnt, map_value;
    float probability;
    
    //sta variables
    int nframes, ncells, sta_length;
    int frames, shift, cell;
    
    float *all_sta;
    
    //get seed
    state = (SInt64 *)mxGetData(prhs[0]);
    
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
    
       
    //RGB/BW, Gaussian/Binary, sparse probability flags
    noise_type = mxGetScalar(prhs[6]);
    n_bits = mxGetScalar(prhs[7]);
    probability = mxGetScalar(prhs[8]);
    
    // sta variables
    nframes = mxGetScalar(prhs[9]);
    ncells = mxGetScalar(prhs[10]);
    sta_length = mxGetScalar(prhs[11]);
    
    // spikes
    spikes_array = mxDuplicateArray(prhs[12]);
    spikes = (unsigned char *)mxGetData( spikes_array );
    int spikes_2d[nframes][ncells];
    int i, j, spikes_index;
    for (j=0; j<ncells; j++){
        for (i=0; i<nframes; i++){        
            spikes_index = nframes * j + i;         
            spikes_2d[i][j] = spikes[spikes_index];
//             mexPrintf("%d\t", spikes_2d[i][j]);
        }
//         mexPrintf("%d\n");
    }
    //sta
    float sta[ncells][sta_length][3*width*height];
    for (j=0; j<ncells; j++){
        for (i=0; i<sta_length; i++){
            for (cnt=0;cnt<3*width*height; cnt++){
                sta[j][i][cnt] = 0;
            }
        }
    }

    //associate output
    const mwSize dims_C[]= {3*width*height, sta_length, ncells};
    plhs[0] = mxCreateNumericArray(NDIMS_C,dims_C,mxSINGLE_CLASS,mxREAL);
    all_sta = (float *)mxGetData( plhs[0] );
    
    int map_lined[3][width][height];
    prefilled_seq = (unsigned char *)(map_lined);
    
    int if_fill = 1;
    
    // create movie
    for (frames = 0; frames<nframes-sta_length; frames++){  
        
        //create one frame
        for( h=0; h < height; h++) {
            
            image_index = 3 * h * width;
            
            for( w=width; w != 0 ; w--) {
                
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
            }
        }
        
        // fill sta if spikes occured 
        if (frames<(nframes-sta_length)){
            for (cell=0; cell<ncells; cell++){ // goes over cells
                for (shift=0;shift<sta_length; shift++){ // goes over sta length
                    for (i=0;i<spikes_2d[frames + shift][cell];i++){
                        for (cnt = 0; cnt<3*width*height; cnt++){
                            sta[cell][shift][cnt] += (((float) prefilled_seq[cnt])-128)/256;
                        }
                    }
                }
            }
        }
    }
    

    memcpy(all_sta, sta, sizeof(float)*ncells*sta_length*3*width*height);

    //    mexPrintf("\n%d\t %d ", cnt_spikes, spikes[cnt_spikes] );

               
    return;
} 
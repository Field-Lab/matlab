#include "io64.h"
#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define NDIMS_A 2
#define NDIMS_B 3



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{

	FILE *fp;
	char *filename;
	size_t elements_read; 
    unsigned char *read_vector;
    int width, height, header_len, frame_num, buflen, status;

    int64_T offset;
	unsigned char *image_mat; // unsigned char is uint8
	int x, w, h, read_index, fill_index, flip_flag;

    
    // input protection
	if( nrhs != 6 ) {
		mexErrMsgTxt("Need exactly six inputs ... (filename, width, height, header_len, current_frame_num, flip_flag)");
	}
	if( nlhs > 2 ) {
		mexErrMsgTxt("Too many outputs.");
	}
	if( !mxIsChar(prhs[0]) ) {
		mexErrMsgTxt("Argument must be base filename string.");
	}

    /* Get the length of the input string. */
    buflen = (mxGetNumberOfElements(prhs[0]) + 1);

    /* Allocate memory for input and output strings. */
    filename = mxCalloc(buflen, sizeof(char));

    
    //associate pointers
   // sp = mxGetData(prhs[0]);
    status = mxGetString(prhs[0], filename, buflen);
    
    // Retrieve input scalar size values
    width = mxGetScalar(prhs[1]);
    height = mxGetScalar(prhs[2]);
    header_len = mxGetScalar(prhs[3]);
    frame_num = mxGetScalar(prhs[4]);
    flip_flag =  mxGetScalar(prhs[5]);
    
    const mwSize arraysize[]={(width*height*3), 1};
    
	fp = fopen(filename, "r");    // rb since we know it isnt a text file. read binary
	
    if( !fp ) {
		mexErrMsgTxt("File will not open.");
	}
    
    
    // fseek goes here
    //fseek(fp, (header_len * sizeof( char )), 0);
    //fseek(fp, (header_len * sizeof( unsigned char )), 0);
    //fseek(fp, (((frame_num * arraysize[0]) + header_len) * sizeof(unsigned char)), 0);
	    
    // http://www.gridlabd.org/documents/doxygen/1.1/mexfio64_8c-source.html
    offset = (int64_T)((((int64_T)frame_num * (int64_T)arraysize[0]) + (int64_T)header_len) * (int64_T)sizeof(unsigned char));
	
    x = setFilePos( fp, (fpos_T*) &offset);
    
    // Create a space to put the date
    
    plhs[0] = mxCreateNumericArray( NDIMS_A, arraysize, mxUINT8_CLASS, mxREAL );  // will need to change format
    
    
    // Associate the pointer
	read_vector = (unsigned char *)mxGetData( plhs[0] );
    
    elements_read = fread( read_vector, sizeof(unsigned char), arraysize[0], fp );
    
    fclose(fp);
    
	const mwSize dims_B[]= {4,width,height};
    plhs[1] = mxCreateNumericArray(NDIMS_B,dims_B,mxUINT8_CLASS,mxREAL);
    
    image_mat = (unsigned char *)mxGetData( plhs[1] );

	
	// Now we have the goods, time to fill the matrix
	// we are emulating the reshape command of matlab
	// movie_frame(1:3,:,:) = reshape(temp, 3, stim_obj.span_width, stim_obj.span_height); 
	//  - this means the movie pixels are read off from the vector as little columns of threes
	// 
	// from random matrix fill
    
    for( h=0; h<height; h++) {  
        
        for( w=0; w<width; w++) {
                    
                    fill_index = (h * (4 * width)) + (4 * w);
                    
                    switch (flip_flag) {
                        case 1:         // Normal: no flips
                            read_index = (h * (3 * width)) + (3 * w);
                            break;
                            
                        case 2:         // Flip vertical
                            read_index = ((height - h) * (3 * width)) + (3 * w);
                            break;
                            
                        case 3:         // Flip left-right
                            read_index = (h * (3 * width)) + (3 * (width - w));
                            break;
                            
                        case 4:         // Flip left-right AND vertical
                            read_index = ((height - h) * (3 * width)) + (3 * (width - w));
                            break;
                            
                    }  // flip switch
                    
                    // Fill out matrix
                    image_mat[ fill_index ] = read_vector[ read_index ];  //  R
                
                    image_mat[ fill_index+1 ] = read_vector[ read_index+1 ];  // G
                 
                    image_mat[ fill_index+2 ] = read_vector[ read_index + 2];  // B
                 
                    image_mat[ fill_index+3 ] = (unsigned char)0xff;     // Fill the 4th color (alpha) channel with 255 as a filler
                  
        }       // height
    }           // width

    return;
    
}
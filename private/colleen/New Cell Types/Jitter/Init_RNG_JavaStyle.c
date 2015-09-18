#include <matrix.h>
#include <mex.h>   
#include <Random.h>

#define NDIMS 2


// state is a pointer to a 1-element SInt64 location in memory
// The basic operation is to set the value stored in the the state-address to the seed mumbo-jumbo


void RandJavaSetSeed(SInt32 seed, RandJavaState state) 
{
	*state = (seed ^ 0x5DEECE66DLL) & 0xFFFFFFFFFFFFLL;
	return;
}

// NB on the bizzaro strings: 0x is merely the historical prefix used by C to tell the parser that a hex number is what follows
// The LL suffix stands for the fact that these are _literals_ of the size Long-long int

// hex: 5DEECE66D = decimal 25214903917 = 010111011110111011001110011001101101 bin.         [9 hex characters]
// hex: FFFFFFFFFFFF = decimal 281474976710655 = 111111111111111111111111111111111111111111111111 binary [12 hex characters]

// So: we take bitwise XOR with seed. Then do a bitwise AND with the FFFFF-ect. 
// This is what gets put into the location whose pointer is state. 
     
/*  perform a bitwise XOR operation ("^")
  with the number: 25214903917 (which is given in hexadecimal notation by: 5DEECE66D)
  whose binary representation is: 010111011110111011001110011001101101

  (It turns out the 0x prefix merely tells the computer that the following 5DEECE66D is a number in hexidecimal notation. The LL suffix tells the computer that his is a "literal" number. Of the type "long-long integer".)

  then take that mess and perform a bitwise AND operation ("&")
  with the number 281474976710655 (the hex number FFFFFFFFFFFF)
  whose binary representation is: 111111111111111111111111111111111111111111111111

  then store that whole mess in the memory location that is specified by the variable state. 
*/


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])    
{
    /* We take in a seed and return the state variable. Which is the pointer to the memory location where we put our random numbers */
    
    
    SInt32 seed;  
    SInt64 *state;
    const mwSize dims[]={1,1};
  

    //associate outputs
    
    plhs[0] = mxCreateNumericArray(NDIMS,dims,mxINT64_CLASS,mxREAL);
    
    state = (SInt64 *)mxGetData( plhs[0] );
    
    //associate pointers
    
    seed = mxGetScalar(prhs[0]);

    
    //Actually do something
    RandJavaSetSeed(seed, state); 

    return;
}
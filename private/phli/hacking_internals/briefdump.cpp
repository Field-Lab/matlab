#include "mex.h"
#include "mxinternals.h"


/* Print one-liner describing mxArray, including any crosslinks */
void briefdump(mxArray *in)
{
  int *tmp;
  int i;

  if (in->name && *in->name)
    printf("[%.32s]", in->name);
  printf("(%8s)  ", mxGetClassName(in));
  printf("Address: %p", in);

  if (in->crosslink)
    printf("  <linked to %p>", in->crosslink);

  if (in->refcount)
    printf("  <refcount %d>", in->refcount);

  if (in->number_of_dims == 2) {
    printf("   [%i %i]", in->rowdim, in->coldim);
    if (in->rowdim*in->coldim == 1 && mxGetClassID(in)==mxDOUBLE_CLASS)
      printf(" scalar value: %f", *((double*)(in->data.number_array.pdata)));
    printf("\n");
  } else {
    printf("   [");
    tmp = (int *)(in->rowdim);
    for (i=0; i<in->number_of_dims; i++) {
      printf("%i ", tmp[i]);
    }
    printf("]");
  }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  if(nrhs < 1) 
    mexErrMsgTxt("One input required.");  

  briefdump((mxArray*) prhs[0]);
}
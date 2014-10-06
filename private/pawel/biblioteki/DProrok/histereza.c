#include "mex.h"
#include "math.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
const mxArray *prhs[])
{
int m,n,i;
double *A, *B;
/* check for arguments here, omitted */
m = mxGetM(prhs[0]);
n = mxGetN(prhs[0]);
A=mxGetPr(prhs[0]);
/* allocate the answer */
plhs[0] = mxCreateDoubleMatrix(n, m, mxREAL);
B = mxGetPr(plhs[0]);
for(i=0;i<m*n;i++)
{
    if(*(A+i)==0){*(B+i)=0;}
    else if(*(A+i)==1)
        {
        *(B+i)=1;
        i++;
        for(i;i<m*n;i++)
        {
            if(*(A+i)!=2){*(B+i)=0;}
            else
                break;
        }
        }
    else *(B+i)=0;












}
}

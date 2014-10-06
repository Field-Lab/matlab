#include <mex.h>
#include <cmath>

/*
 * Was segfaulting when test created local array for results; probably was
 * killing the stack or else stack was going out of scope while other
 * threads still accessing it; could test adding a barrier/join?
 */

void filltable(double *table, int size) { 
    #pragma omp parallel for
    for (int n = 0; n < size; ++n)
        table[n] = std::sin(2 * M_PI * n / size);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    int size = 1000000;
    plhs[0] = mxCreateDoubleMatrix(size, 1, mxREAL);
	double* table = mxGetPr(plhs[0]);
    filltable(table, size);
}
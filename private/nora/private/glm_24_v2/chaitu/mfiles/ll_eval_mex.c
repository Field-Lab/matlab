/*=================================================================
 *
 * ll_eval_mex.C    .MEX file for calculating the cifs for given spike trains
 * input2 = ll_eval_mex(D,input,PS,CP)
 *
 * Inputs:  D - list of spikes - first col is time index, second col is neuron
 *          input - filtered stimulus + base rate for each neuron (T X Nneuron)
 *          PS - postspike filters per neuron (Mhist x Nneuron)
 *          CP - coupling filters per neuron (Mcoup x Nneuron*(Nneurons_eff-1))
 *          baseidx - indices of spike trains (columns of D) that corresp. to cols of input
 *
 *Outputs: input2 - the log of the resulting CIF for these spike trains (T x Nneuron)
 *=================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  double *D, *input1, *ps, *cp, *input2, *filt_pr, *baseidx;
  int maxt,nneurons,nneurons_eff,mh, mc, nsp, ub_ps,ub_cp;
  int i,j,k,sp_time, sp_neuron, neuron_i_idx;

  if (nrhs < 3) { 
    mexErrMsgTxt("need at least 3 input args:  ll_eval_mex(D,input1,PS)");
 }

  /*  Read inputs into local variables */
  nsp = mxGetM(prhs[0]);
  D = mxGetPr(prhs[0]);
  
  input1 = mxGetPr(prhs[1]);
  maxt = mxGetM(prhs[1]);  
  nneurons = mxGetN(prhs[1]);
  nneurons_eff = (int)(D[nsp*2-1]); /* mxGetN(prhs[0]); */
  /* printf("nneurons_eff is %d\n",nneurons_eff); */
  
  mh = mxGetM(prhs[2]);
  ps = mxGetPr(prhs[2]);

  if (mxGetN(prhs[0]) != 2) {
      mxErrMsgTxt("invalid format for spike list!\n");
  }
  
  if (mxGetN(prhs[2]) != nneurons) {
      mxErrMsgTxt("invalid dimensions for ps filters!\n");
  }
  
  if (nneurons_eff > 1) {
      if (nrhs < 3) {
          mxErrMsgTxt("need a coupling filter!\n");
      }
      if (mxGetN(prhs[3]) != nneurons*(nneurons_eff-1)) {
          mxErrMsgTxt("invalid dimensions of coupling filters. mdim = %d, ndim = %d, %d, %d!\n",mxGetM(prhs[3]),mxGetN(prhs[3]),mh,nneurons*(nneurons_eff-1));
      }
      mc = mxGetM(prhs[3]);
      cp = mxGetPr(prhs[3]);
      
      if (nrhs < 4) {
          mxErrMsgTxt("need to specify base neuron indices!\n");
      }
      if (mxGetM(prhs[4]) != nneurons) {
          mxErrMsgTxt("invalid number of base indices! must be same as number of cols in input matrix!\n");
      }
      baseidx = mxGetPr(prhs[4]);
      
      /* ***
      printf("The neuron indices are: ");
      for (i=0;i<nneurons;i++) {
          printf("\t %d",(int)(baseidx[i]));
      }
      printf("\n");
      ***/

  }
  
  /* Create output variable */
  plhs[0] = mxCreateDoubleMatrix(maxt,nneurons,mxREAL);
  input2 = mxGetPr(plhs[0]);
  /* Copy the input */
  memcpy((void *)input2,(void *)input1,sizeof(double)*maxt*nneurons);

 for (k=0;k<nsp;k++) { /* Iterate through spikes */
    sp_time = (int)(D[k]-1); /* time of spike */
    sp_neuron = (int)(D[k+nsp]-1); /* neuron that spiked (index in spike list) */

     /* printf("Found spike by neuron %d at time %d\n",sp_neuron,sp_time); */
    
    /* Determine the timerange to be updated */
    ub_ps = maxt;
    if (sp_time+mh+1 < maxt) {
        ub_ps = sp_time+mh+1;
    }
    
    ub_cp = maxt;
    if (sp_time+mc+1 < maxt) {
        ub_cp = sp_time+mc+1;
    }
    
    /* Iterate through neurons (up to nneurons) affected by this spike */
    for (i=0; i < nneurons; i++) {
        
        if (nneurons_eff > 1) {
            neuron_i_idx = (int)baseidx[i]-1; /* index of neuron i in the spike list */
        } else {
            neuron_i_idx = 0;
        }
                
        /* printf("neuron_i_idx is %d\n",neuron_i_idx); */
        
        if (neuron_i_idx == sp_neuron) { /* neuron baseidx[i] spiked */
            filt_pr = ps + i*mh;

            for (j=sp_time+1;j<ub_ps;j++) {
               input2[i*maxt+j] += filt_pr[j-sp_time-1];
            }
            
            /* printf("Updating postspike term of neuron %d\n",neuron_i_idx+1); */
        } else {
            if (sp_neuron < neuron_i_idx) { /* a neuron before i spiked */
                filt_pr = cp + mc*(i*(nneurons_eff-1)+sp_neuron);
            } else { /* a neuron after i spiked */
                filt_pr = cp + mc*(i*(nneurons_eff-1) + sp_neuron-1);
            }

            for (j=sp_time+1;j<ub_cp;j++) {
                input2[i*maxt+j] += filt_pr[j-sp_time-1];
            }
            
           /* printf("Updating coupling term of neuron %d to %d\n",sp_neuron,neuron_i_idx+1); */
        }
        
    }
    
 }
  
 


}




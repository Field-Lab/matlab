% NEW NAMING OF GLM_WRAP AND GLM_EXECUTE

% 2015-07-14
% Both moved to main glm folder because they are the sites of major edits.
% glm_wrap should only call glm/wrap_bookkeep
% glm_execute should only call glm/glm_core



% Version 3: 
% glm_wrap:
%   - add load input argument to the runoptions
%   - made debug version actually run on shorter input (just 5 blocks)
%   - fixed corresponding test call
%   - add a 50% fit option
%    Need to work in outputNL fitting!
% glm_exectue:
%   - none
% glm_execute_OutputNL is now up and running


% Version 2: 
% glm_execute: started 2015-07-14
%   - enable optional_arg input,(for GLMPars for now)
%   - got rid of troubleshoot input
%   - register input nonlinearities directly to fittedGLM
%   - optional_arg as means of better initialization
%   - enable one line plot note to enter into printglmfit
%   - change direcotories to print
% glm_wrap:
%   - added separate call to glm_execute_InputNL_IteratedOpt

% Version 1: started 2015-07-13  adopted 2015-07-14
% glm_wrap: 
%   -added a special_arg input designed to handle concepts such as
%   constrained searching, and fitting non-linearities
%   -took out unused call to doubleopt
% glm_execute:
%   -enables for constraining PS filter using change of basis, fmincon
%   -cleaner calling of fitting algorithm
%   -PS constrain should work robustly for all cases

% Version 0: up to 2015-07-14
% glm_execute: Coupling, xval measures, printing, rk2, rk1 all integrated
% glm_wrap: loads cells, sets directories, loads and process spikes and stim




%%% BELOW IS ALL BEFORE 2015-04-01
%{
GLM starts off at GLM_AH Version 6_2
GLM started on 2015-03-30

Goals:
Boot out the Analysis Directory
Clean up datarun nomenclature
get working for rk1 and rk2
then look at input pt non-linearities
combine loadmoviemat with StimulusPars (Directories_Params_v24)
make the wrap just two files


to do: working copy of rk1 rk2
get rid of UOP
add mfilestamp to print function


clear; clc;
restoredefaultpath; 
glmpath_Alligator
exps = [1]; 
stimtypes = [1];
celltypes = [2];
cell_subset = 'debug'; 
glm_settings{1}.type = 'debug'
glm_settings{1}.name = 'true';
runoptions.replace_existing = true;
glm_wrap(exps,stimtypes,celltypes,cell_subset,glm_settings,runoptions)
                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            49.1197                           560                
     1            49.1197             10            560           4
     2           -987.807            2.5             68           0
     3           -1054.74         2.3853           51.9           3

Version 5c All traces of conductance based eliminated.
Version 5b glmwrap cleaned up, commented, some name changes 
Version 5a cut out all noncalled scripts and functions
Version 5  cut out non-essential code!


Version 4 works for rk1 rk2 and STA
All OLD CODE IS STILL THERE!!!!


Version 3_e
Further cleaned Directories_Params_v24


Version 3_d.
Moved Datarun to "more appropriate locations" within Blocked Spikes dir
All four datasets should work


Version 3_c (no longer work)
Got rid of Java path  (glmpath)
Introduce manual loading of dataruns (this is better!!)
Directories_Params_v24 now just does Dirs, params,  no datarun
Works for White noise 2012-08-09-3 .. need to get better about which 
directories to use.

Version 3_b Working Copy!
XVAL is now performed in glm_execute
works for STA

Version 3_a not a working copy!
added computations to glmwrap_func for testmovies, testspikes
made variable naming in terms of test fit cleaner in glmwrap_func
glm_execute now doesn't work

Version 3: Bring XVAL CODE into glm_execute

Version 2_a
glm_core isolated to glm_execute and its dependencies
path eliminates sub-directories of glm_core / chaitu
- STA works up to xval, need to get rk2 up as well

Version 2: Work on glm_execute (all inside glm_core) bring

%%%%%%%%%%%%%
Version 1_d
Moved glm_code computations to glm_core

Version 1_c
Most other subfolder moved into oldcode_dontdeleteyet
glmwrap_func completely independent (asides from wrap_bookkeep)

Version 1_b
Formalized path changes glmpath_Alligator
finished sealing glmwrap_func
wrap_bookkeep  DONE
bookkeep eliminated

Version 1_a
Brought generalcomputation_AKH and External_Code into the glm folder
Subroutined some of Directories_Params_v23

Version 1: Resort Code: wrap_computations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Version 0
Works! 2015-03-30
new folder name  Matlab_code/glm
output:  NSEM_Home/GLMOutput_Raw

clear; clc;
restoredefaultpath; 
glmpath_Alligator
exps = [1]; 
stimtypes = [1];
celltypes = [2];
cell_subset = 'debug'; 
glm_settings{1}.type = 'debug'
glm_settings{1}.name = 'true';
runoptions.replace_existing = true;
glm_wrap(exps,stimtypes,celltypes,cell_subset,glm_settings,runoptions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Output: exp 1 stimtyp1 celltypes 2  debug  cid 1471
### running: WN expA OFFPar_1471: debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams ###
### fit type: debug_fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams ###
                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            49.1197                           560                
     1            49.1197             10            560           4
     2           -987.807            2.5             68           0
     3           -1054.74         2.3853           51.9           3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


exps = [3]; 
stimtypes = [2];
celltypes = [1];
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'rk2';
cell_subset = 'shortlist';
runoptions.replace_existing = true;
glmwrap_func(exps,stimtypes,celltypes,cell_subset,changes_cell,runoptions)





clear
exps = [1]; 
stimtypes = [1];
celltypes = [2];
changes_cell{1}.type = 'filter_mode';
changes_cell{1}.name = 'rk1';
changes_cell{2}.type = 'debug'
changes_cell{2}.name = 'true';
cell_subset = 'debug';
runoptions.replace_existing = true;
glm_wrap(exps,stimtypes,celltypes,cell_subset,changes_cell,runoptions)
### running: WN expA OFFPar_1471: debug_rk1_MU_PS_noCP_p8IDp8/standardparams ###
### fit type: debug_rk1_MU_PS_noCP_p8IDp8/standardparams ###
                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            50.6303                           560                
     1            50.6303        9.23728            560           3
     2           -836.639        2.30932           9.57           0
     3           -836.639        4.61864           9.57           2
     4           -880.634        1.15466           60.2           0




MORE OUTPUT THAT IS NOT DEBUGGED
MORE OUTPUT THAT IS NOT DEBUGGED
MORE OUTPUT THAT IS NOT DEBUGGED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

### running: WN expA OFFPar_1471: fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams ###
### fit type: fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams ###

                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            1379.38                       1.8e+04                
     1            1379.38             10        1.8e+04           4
     2           -33171.3            2.5       1.42e+03           0
     3           -35078.5        2.39032       1.33e+03           3
     4           -35272.2        1.23018            126           3
     5           -35310.3         1.8414           26.9           7
     6           -35321.4        2.60222           14.9          13
     7           -35325.2        3.34606           4.44          25
     8             -35327        2.01323           4.82          25
     9           -35327.8        1.19082            1.8          25
    10           -35328.2        1.07951           3.06          25
    11           -35328.4       0.447972          0.988          25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### running: WN expA OFFPar_1471: rk1_MU_PS_noCP_p8IDp8/standardparams ###
### fit type: rk1_MU_PS_noCP_p8IDp8/standardparams ###
                                Norm of      First-order 
 Iteration        f(x)          step          optimality   CG-iterations
     0            1395.72                       1.8e+04                
     1            1395.72        10.4007        1.8e+04           3
     2           -28688.1            2.5       1.74e+03           0
     3           -28688.1              5       1.74e+03           2
     4           -30102.1           1.25       1.77e+03           0
     5           -32232.3            2.5       9.74e+03           1
     6           -34766.5          0.625       2.74e+03           3
     7           -35293.8        1.05612            574           4
     8           -35343.2           1.25           46.6           8
     9           -35358.1        2.77749           12.5          17
    10           -35364.3              5           23.3          46
    11           -35365.4        2.60641          0.704          45
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%}


%%%% GLM_AH VERSTION TRACKER %%%%
%{

clear
exps = [1]; 
stimtypes = [1];
celltypes = [2];
changes_cell{1}.type = 'debug'
changes_cell{1}.name = 'true';
cell_subset = 'shortlist';
runoptions.replace_existing = true;
glmwrap_func(exps,stimtypes,celltypes,cell_subset,changes_cell,runoptions)

%%% GLM_AH  Version Tracker %%%%
Version 6_2  output goes to glm_AH_

Version 6_2_d change the Bertha directories (Volumes/Labs/Users/akheitman)
Version 6_2_c worked as of 2015-03-29
Get input non-linearities to work
Option higher tolx values for rk2  (not -4 but rather -6)
Turns out rk2 tolfun -4 is SUFFICIENT!
add runotpion to glmwrap_func
add special change competency to the wrap
got rid of display "NewHessCorrection"


Version 6_1
Up and working .. Running renewed rk2 fits here! (tolx of 4)
Finished up rast_bps_findPS

Version 6_0
Change  glm_wrap_func calling sequence
Change "generate_optimalcondtioned folder" to "raster_performance"  


Version 5_3
Everything works for rastbps_wrap  (still need to do findPS filter, almost
done)

Version 5_0 
Added generate_optimalconditioned folder

Versin 4_2
just more plots and better glm_convergence
2015-02-09

Version 4_1
New folder for glm_convergence
More work on plots
Saved 2015-02-05

Version 4_0
Set up PrototypePlots folder in the NSEM_Home as input/output of Plots
Now works on Bertha with new NSEM_Home and Code Home!!

Version 3_3
continue work in Plots_2015


Version 3_2
New Hacks to eval_xvalperformance_NEW version 1 used... rather than version 0
Once again revisit making plots of comparison
New Folder Plots_2015 (shouldn't be in the path)

Version 3_1_A
No actual alterations
Commented
- glm_convex_optimizationfunction_withNL
- glm_opt_NonConvex_function
- maxiter set to 300;

Version 3_1 
Should only affect rk2-ConductanceBased
Big Change in  glm_opt_NonConvexfunction_withNL (evaluation of lcif_derivative_NL)
Seems to be working ... racking up big counts of CG-iterations
maxiter at 100;

Version 3_0 working
Get rk2-Conductance Based up and running!
NEW m-files
glm_opt_NonConvexfunction_withNL
maxiter at 100;
rk2-Conductance Based  not converging for NSEM
rk2-ConductanceBased WN converges ~30 but poor PS Filter

Version 2_3
Hess Correct back on .. got rid of .5
Max_Iter reset to 100
Big changes to glm_nonconvex for correect gradient and hessian evluation
Everything back under control!! Fast Convergence rk1 NSEM in 20-40 iterations!
HUGE SUCCESS: fast convergence of NSEM rk1/rk2 (~30 iterations)
NEW m-files
glm_opt_NonConvex_function
glm_nonconvex_HessianCorrection


Version 2_2
Restore conditions of glm_AH_24 to see if NSEM will converge with rk1
tolfun universally set to 5
Hessian correction turned off
init conditions set to .01 for each param
save fminunc output so we can track iteration numbers 
    inside glm_execute: fittedGLM.fminunc_output = output; 


Version 2_1
Start commenting out the glm_nonconvex_optimizationfunction
New Hess Correction
Performance:
- NSEM rk1,rk2 fits aren't finishing in under 250 iterations
- NSEM fits need a better initial condition! + iteration tracker!
- all WN fits are just fine!
- NSEM fixed SP fits are still good

Version 2_0
Non-zero initial conditions.. so non-convex filters don't get stuck! WORKS!
Commented glm_convex_optimizationfunction
Old Hess_Correction

Versions 2-3 largely working within
glm_nonconvex_optimizationfunction
and its corresponding glory

Version 2
Re-examine Hessian Correction

Version 1_2
Turn Hess correct back on
Zeros for initial condition .. non-convex filters get stuck!

Versions 1_1:
Made p8IDI8p the default stimulus
Hess Correct Off

Version 1_0:
Got it up and running
%edit inside glm_wrap_func
GLM_Setting
%edit inside glm_exectue%
glm_execute
prep_paramind
prep_stimcelldependentGPXV
glm_convex_optimizationfunction_withNL (unchagned..used ConductanceBased-HardRect)
Hess Correct Off

Version 1: 
Work on fixedSP-ConductanceBased
Will be a convex version to explore Conductance Based for NSEM
Hess Correct Off

Version 0_4
glmwrap deleted
p_init set to a vector of zeros

Version 0_3
No Hess Correction
Corrected Calling Sequence for Conductance Based Model
Testing on Conductance Based
Done 2-15-01-13

Version 0_2 
Confirmed to work with rk1 and rk2.
Does not work with Conductance Based Model 
Has no Hess Correction.
HESS CORRECTION CONFIRMED TO CONFUSE RATHER THAN HELP FOR WN
NOT DONE CORRECTLY YET
Done 2015-01-13

Version 0_1
Works.  
Also runs with Conductance Based but not with same precision as
before winter break.
Added glmwrap_func for easier calling from Bertha.
Done 2015-01-11


Version 0 
Works. Adopt Rokhlin style Versioning system.
Done 2015-01-06

%}

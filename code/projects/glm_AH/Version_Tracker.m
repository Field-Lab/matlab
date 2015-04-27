%{

Version 6_2
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

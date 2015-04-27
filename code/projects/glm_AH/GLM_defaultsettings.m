

GLMType.fit_type = 'WN'; GLMType.map_type = 'mapPRJ'; 
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
end

if i_run == 2
GLMType.fit_type = 'NSEM'; GLMType.map_type = 'mapPRJ'; 
%GLMType.cone_model = 'linear_timekernel_shift0';  GLMType.cone_sname ='timekernel';
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';
GLMType.cone_model = 'DimFlash_092413Fc12_shift0'; GLMType.cone_sname = 'timekernelCONEMODEL';
end



%GLMType.cone_model = 'linear_timekernel_shift0';  GLMType.cone_sname ='timekernel';
%GLMType.cone_model = '8pix_Identity_8pix'; GLMType.cone_sname='p8IDp8';
%GLMType.cone_model = '8pix_Model1_1e4_8pix'; GLMType.cone_sname = 'p8Mod1Max1e4p8';
%GLMType.cone_model = 'linear_timekernel_shift1';  GLMType.cone_sname = 'timekernel';


%GLMType.cone_model = 'TimeConvolve_DimFlash_shift2'; GLMType.cone_sname = 'timesmoothed';
%GLMType.k_filtermode = 'OnOff_hardrect_fixedSP_STA'; GLMType.fixedSPlength = 13;  GLMType.fixedSP_nullpoint = 'mean'; 
GLMType.nullpoint = 'mean'; 
%GLMType.stimfilter_mode = 'rk2-ConductanceBased';
GLMType.stimfilter_mode = 'rk2';
%GLMType.stimfilter_mode = 'fixedSP_rk1_linear';
%GLMType.input_pt_nonlinearity      = false;
%GLMType.input_pt_nonlinearity_type = 'piece_linear_aboutmean';
%GLMType.input_pt_nonlinearity_type = 'piece_linear_shiftmean';
%GLMType.input_pt_nonlinearity_type = 'polynomial_order5_part4';
%GLMType.input_pt_nonlinearity_type = 'piecelinear_fourpiece_eightlevels';
%GLMType.input_pt_nonlinearity_type =  'oddfunc_powerraise_aboutmean';
%GLMType.input_pt_nonlinearity_type =  'log';
%GLMType.input_pt_nonlinearity_type =  'exp';
%GLMType.input_pt_nonlinearity_type = 'polynomial_androot_order2_search2';  % order plus minus 2
%GLMType.input_pt_nonlinearity_type = 'polynomial_androot_order2_search3';

%GLMType.postfilter_nonlinearity      =  true;
%GLMType.postfilter_nonlinearity_type = 'ConductanceBased_HardRect';
%GLMType.postfilter_nonlinearity_type =  'oddfunc_powerraise_aboutmean';
%GLMType.postfilter_nonlinearity_type =  'piece_linear_aboutmean';
%GLMType.postfilter_nonlinearity_type = 'raisepower_meanafter';
%Type = 'Stim_Nonlinearity'; modification = 'raisepower_meanafter'
GLMType.CONVEX = false; % with relation to the filters .. are parameters used linearly in the GLM. 
%GLMType.DoubleOpt = true;
%GLMType.DoubleOpt_Manual = true;
%GLMType.stimfilter_mode = 'rk1';
GLMType.specialchange = false;
%GLMType.specialchange_name = 'Conductance_Based';
%}

GLMType.TonicDrive = true;
GLMType.StimFilter = true;
GLMType.PostSpikeFilter = true;
GLMType.CouplingFilters = false;
GLMType.fixed_spatialfilter = true;
GLMType.func_sname = 'glmwrap_24';
GLMType.fullmfilename =mfilename('fullpath'); 
i_exp = 1; i_cell = 1;

GLMType.fitname  = GLM_fitname(GLMType);   
troubleshoot.doit    = true;
troubleshoot.plotdir = '/Users/akheitman/Matlab_code/troubleshooting_plots'
troubleshoot.name    = 'singleopt';
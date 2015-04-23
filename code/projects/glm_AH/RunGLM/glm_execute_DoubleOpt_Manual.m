% The DoubleOpt analogues of glm_execute
% Designed to implement a second optimization (usually nonlinearity to the
% stimulus)
% 
% Hack Code .. just to make it work 2014-05-24
% Updated 2014-06-10.   Works smoothly now
%  for 1 parameter searches uses fminbnd rather than fminunc
%  fminbnd is much faster here. when provided with reasonable bounds

% BND uses a different optimizing routine which is bounded


% CALLS which use GLMType:
%  prep_paramindGP
%  prep_stimcelldependentGPXV
%  glm_execute_DoubleOpt_InnerLoop

function [fittedGLM, manual_search] = glm_execute_DoubleOpt_Manual(GLMType, spikes, fitmovie, inputstats, glm_cellinfo,troubleshoot)
%% Get rid of all time,  put all inputs into bins

fittedGLM.cellinfo = glm_cellinfo;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load up GLMParams compute some universal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GLMPars           = GLMParams;


if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end

if GLMType.DoubleOpt, GLMPars.optimization.tolfun = 4; end
if GLMType.debug, GLMPars.optimization.tolfun = 1; end
fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;

frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;


% Perhaps we should combine this! With convolving with spikes !
bin_size      = t_bin;
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
if GLMType.CouplingFilters
    basis_params  = GLMPars.spikefilters.cp;
    cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
end
clear bin_size basis_params

% Convolve Spike Times with appropriate basis
% Think about flushing dt out to the wrapper
% Take care of all timing in glm_execute or in glmwrap.
t_bin        = t_bin;
home_sptimes = spikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    fixed_covariates.PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.CouplingFilters;
    basis = cp_basis';
    display('figure out coupling here!  CP_bin');
end

if GLMType.TonicDrive
    fixed_covariates.MU_bin = ones(1,bins);
end

[paramind] =  prep_paramindGP(GLMType, GLMPars); 
%p_init     =  .01* ones(paramind.paramcount,1);
%p_init     = zeros(paramind.paramcount,1);
p_init     = .01* ones(paramind.paramcount,1);



% Find the correct stimulus related input term


%%
GLMType_start = GLMType;
GLMType_start.input_pt_nonlinearity = false;
[fstar, pstar] = glm_execute_DoubleOpt_InnerLoop([], '',glm_cellinfo, GLMType_start, GLMPars,fixed_covariates,home_spbins,fitmovie,inputstats, t_bin,p_init);

linear_pstar = pstar; linear_fstar= fstar;
p_guess = pstar;
clear fstar pstar
%%


if isfield (GLMType,'input_pt_nonlinearity')  && GLMType.input_pt_nonlinearity
    if strcmp(GLMType.input_pt_nonlinearity_type , 'piece_linear_aboutmean')
        nonlinpar = [.2:.2:2]';
        lin_index = find(nonlinpar == 1);
    end
     if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
        nonlinpar_ratio = [.1:.1:2]';
        nonlinpar_shift = [-.25:.05:.25]';
        
        manual_search.nonlinpar_ratio = nonlinpar_ratio;
        manual_search.nonlinpar_shift = nonlinpar_shift;
        
        
        nonlinpar0 = .133*ones(length(nonlinpar_ratio)*length(nonlinpar_shift),2);
        nonlinpar0(:,1) = repmat(nonlinpar_ratio,[length(nonlinpar_shift), 1]);
        
        col2 = repmat(nonlinpar_shift, [1, length(nonlinpar_ratio)]);
        col2 = reshape(col2', [length(nonlinpar_ratio)*length(nonlinpar_shift),1]);
        nonlinpar0(:,2) = col2;
        
        nonlinpar = nonlinpar0;
        lin_index = intersect(find(nonlinpar(:,1) == 1) , find(nonlinpar(:,2) == 0) );
     end
    
     if strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_order3_part5')
        % Call External Code ext_partition
        % each row is different parameter set
        % each column is weight attributed
        % Col 1 Linear, Col 2 Quadratic, Col 3 SqRoot
        % Col 4 Third Order, Col 5 3rd root
        nonlinpar = (1/5) * ext_partitions(5,[1,1,1,1,1]);
        lin_index = find(nonlinpar(:,1) == 1);

        manual_search.nonlinpar_polyweights = nonlinpar;
        manual_search.nonlinpar_columns = '[Linear, Quad, Sqroot, Third, ThirdRoot]'; 
     end
    
     if strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_order5_part4')
        % Call External Code ext_partition
        % each row is different parameter set
        % each column is weight attributed
        % Col 1 Linear, Col 2 Quadratic, Col 3 SqRoot
        % Col 4 Third Order, Col 5 3rd root
        nonlinpar = (1/4) * ext_partitions(4,[1,1,1,1,1,1,1,1,1]);
        lin_index = find(nonlinpar(:,1) == 1);

        manual_search.nonlinpar_polyweights = nonlinpar;
        manual_search.nonlinpar_columns = '[Linear, Quad, Sqroot, Third, ThirdRoot]'; 
     end
    
     
    if strcmp(GLMType.input_pt_nonlinearity_type,'piecelinear_fourpiece_eightlevels')
        
        nonlinpar = (4/8) * ext_partitions(8, [1 1 1 1]);
        lin_index = intersect( intersect(find(nonlinpar(:,1) == 1) , find(nonlinpar(:,2) == 1) ) , ...
                    intersect(find(nonlinpar(:,3) == 1) , find(nonlinpar(:,4) == 1) ) ) ;
        
        manual_search.nonlinpar_polyweights = nonlinpar;
        manual_search.nonlinpar_columns = '[slopeq1,slopeq2,slopeq3,slopeq4]'; 
    end
     
     
     
    
end

objective_value = zeros( size(nonlinpar,1) , 1);
opt_pval = zeros(size(nonlinpar,1),length(p_guess));




%% Run through
if isfield(GLMType,'input_pt_nonlinearity')  && GLMType.input_pt_nonlinearity 
    for i_try = 1:size(nonlinpar,1)
        [fstar, pstar] = glm_execute_DoubleOpt_InnerLoop(nonlinpar(i_try,:), 'input_pt_nonlinearity',...
            glm_cellinfo, GLMType, GLMPars,fixed_covariates,home_spbins,fitmovie,inputstats,t_bin,p_guess);
        
        if rem(i_try,20) == 1        
            display(sprintf('Working on Parm index %d out of %d: opt_value of %d', i_try,size(nonlinpar,1), fstar));
        end
        objective_value(i_try) = fstar;
        opt_pval(i_try,:) = pstar';
    end
end

manual_search.nonlinpar = nonlinpar;
manual_search.dimension = size(nonlinpar,2);
manual_search.lin_index = lin_index;
manual_search.opt_pval  = opt_pval;
manual_search.objective_value = objective_value;
if isfield (GLMType,'input_pt_nonlinearity')  && GLMType.input_pt_nonlinearity
    manual_search.NL = GLMType.input_pt_nonlinearity_type;
end




fstar = min(objective_value);
param_index = find(objective_value == fstar);
nonlinstar = nonlinpar(param_index,:);
nonlin_pstar = (opt_pval(param_index,:))';


%% Same Procedure as glm_execute_Double_Opt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plug optimized parameter back into GLMType so we hold onto them%%%%
if isfield(GLMType,'input_pt_nonlinearity')  && GLMType.input_pt_nonlinearity
	if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        fittedGLM.GLMPars.others.point_nonlinearity.increment_to_decrement = nonlinstar;
    elseif strcmp(GLMType.input_pt_nonlinearity_type , 'raisepower_meanafter')
        fittedGLM.GLMPars.others.point_nonlinearity.scalar_raisedpower = nonlinstar;
    
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
        fittedGLM.GLMPars.others.point_nonlinearity.increment_to_decrement  = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.shiftmean               = nonlinstar(2);
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
        fittedGLM.GLMPars.others.point_nonlinearity.scalar_raisedpower_aboutmean = nonlinstar;
    elseif strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.linear    = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.quadratic = nonlinstar(2);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.sqroot    = nonlinstar(3);
	elseif strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2_search2')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.linear         = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.quadoversqroot = nonlinstar(2);
	elseif strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2_search3')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.linear    = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.quadratic = nonlinstar(2);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.sqroot    = nonlinstar(3);
    
    elseif strcmp(GLMType.input_pt_nonlinearity_type,'polynomial_order3_part5')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.linear    = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.quadratic = nonlinstar(2);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.sqroot    = nonlinstar(3);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.thirdpower= nonlinstar(4);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.thirdroot = nonlinstar(5);
        
    elseif strcmp(GLMType.input_pt_nonlinearity_type,'polynomial_order5_part4')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.linear    = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.quadratic = nonlinstar(2);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.sqroot    = nonlinstar(3);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.thirdpower= nonlinstar(4);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.thirdroot = nonlinstar(5);       
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.fourthpower= nonlinstar(6);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.fourthroot = nonlinstar(7);       
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.fifthpower= nonlinstar(8);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.fifthroot = nonlinstar(9);
        
    elseif strcmp(GLMType.input_pt_nonlinearity_type,'piecelinear_fourpiece_eightlevels')
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.slope_quartile_1 = nonlinstar(1);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.slope_quartile_2 = nonlinstar(2);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.slope_quartile_3 = nonlinstar(3);
        fittedGLM.GLMPars.others.point_nonlinearity.coefficients.slope_quartile_4 = nonlinstar(4);
        
    else
        error('youmessed up assigning the optimized parameters back into your data structure')
    end

    fittedGLM.input_pt_nonlinearity_type   = GLMType.input_pt_nonlinearity_type;
    fittedGLM.pt_nonlinearity_param        = nonlinstar;
end

if isfield(GLMType, 'postfilter_nonlinearity') && GLMType.postfilter_nonlinearity
	lcif_nonlinearity.type = GLMType.postfilter_nonlinearity_type;
	if strcmp(GLMType.postfilter_nonlinearity_type ,  'piece_linear_aboutmean')  
            lcif_nonlinearity.increment_to_decrement = nonlinstar;
            lcif_nonlinearity.filter_index = paramind.X;
            allparams = 1:paramind.paramcount;
            lcif_nonlinearity.linear_index = setdiff(allparams, lcif_nonlinearity.filter_index);

    end
    if strcmp(GLMType.postfilter_nonlinearity_type ,  'oddfunc_powerraise_aboutmean')   
            lcif_nonlinearity.scalar_raisedpower = nonlinstar;
            lcif_nonlinearity.filter_index = paramind.X;
            allparams = 1:paramind.paramcount;
            lcif_nonlinearity.linear_index = setdiff(allparams, lcif_nonlinearity.filter_index);

	end
    fittedGLM.GLMType.lcif_nonlinearity = lcif_nonlinearity;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Unpack the output into filters  %%%%%
% Do this so we don't have to reinterpret the parameters, just have a final
% filter which we can take home and interpret
% this is admittedly a bit ugly .. 


for i_fit = 1:2
    clear linearfilters rawfit
    if i_fit == 1, pstar = linear_pstar; end
    if i_fit == 2, pstar = nonlin_pstar; end
    
    bpf = GLMPars.bins_per_frame;
    rawfit.opt_params        = pstar;
    rawfit.paramind          = paramind;
    rawfit.objective_val     = fstar;


    clear linearfilters
    linearfilters.note = 'These terms, convolved with the covariates, make the log of the conditional intensity function';
    if isfield(paramind, 'MU')
            linearfilters.TonicDrive.Filter      = pstar(paramind.MU);
            linearfilters.TonicDrive.note        ='no convolution necessary, this term is part of the lcif for every bin';
    end
    if isfield(paramind, 'PS')
            rawfit.ps_basis = ps_basis;
            linearfilters.PostSpike.Filter     = ps_basis * pstar(paramind.PS);
            linearfilters.PostSpike.startbin   = 1;  
            linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
            linearfilters.PostSpike.note0       = 'Filter starts at "startbin" bins after the spikebin';
    end
    if isfield(paramind, 'CP')
            rawfit.cp_basis = cp_basis;
            error('need to fill in coupling..  Nice way to handle it')
    end




    center_coord    = glm_cellinfo.slave_centercoord;
    ROI_length      = GLMPars.stimfilter.ROI_length;
    stimsize.width  = size(fitmovie,1);
    stimsize.height = size(fitmovie,2); 
    ROIcoord        = ROI_coord(ROI_length, center_coord, stimsize);
    rawfit.ROIcoord = ROIcoord;
    clear stimsize center_coord;
    WN_STA           = double(glm_cellinfo.WN_STA); 
    [STA_sp,STA_time]= spatialfilterfromSTA(WN_STA,ROIcoord.xvals,ROIcoord.yvals);
    if GLMType.CONVEX
        if strcmp(GLMType.stimfilter_mode, 'fixedSP_rk1_linear')    
            timefilter           = pstar(paramind.X);
            stimfilter           = STA_sp * (timefilter');
            stimfilter           = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.X)]);
            rawfit.spatialfilter = STA_sp;
            linearfilters.Stimulus.Filter             = stimfilter;
            linearfilters.Stimulus.Filter_rank        = 1;
            linearfilters.Stimulus.space_rk1          = reshape(STA_sp, [ROI_length,ROI_length]);
            linearfilters.Stimulus.time_rk1           = pstar(paramind.X);
            %linearfilters.Stimulus.WN_note            = 'use WN STA as a reference to compare to fitted filters'
            %linearfilters.Stimulus.WN_STA             = WN_STA;
            %linearfilters.Stimulus.WN_STA_space_rk1   = reshape(STA_sp, [ROI_length,ROI_length]);
            %linearfilters.Stimulus.WN_STA_time_rk1    = STA_time;
            linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
            linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
            linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
            linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];
            linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
            linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
            linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
        end
    end
    if ~GLMType.CONVEX    
        if strcmp(GLMType.stimfilter_mode, 'rk1') || strcmp(GLMType.stimfilter_mode, 'rk2')       
            timefilter1  = pstar(paramind.time1);
            spacefilter1 = pstar(paramind.space1);
            stimfilter   = spacefilter1 * timefilter1';

            if strcmp(GLMType.stimfilter_mode, 'rk2')
                timefilter2  = pstar(paramind.time2);
                spacefilter2 = pstar(paramind.space2);
                stimfilter  = spacefilter1 * timefilter1' + spacefilter2 * timefilter2';
            end

            stimfilter = reshape(stimfilter, [ROI_length,ROI_length,length(paramind.time1)]);
            linearfilters.Stimulus.Filter             = stimfilter;
            linearfilters.Stimulus.Filter_rank        = 1;
            linearfilters.Stimulus.time_rk1           = timefilter1;
            linearfilters.Stimulus.space_rk1          = reshape(spacefilter1,[ROI_length,ROI_length]);
            linearfilters.Stimulus.x_coord            = ROIcoord.xvals;
            linearfilters.Stimulus.y_coord            = ROIcoord.yvals;
            linearfilters.Stimulus.frame_shifts       = [0:1:(GLMPars.stimfilter.frames-1)];
            linearfilters.Stimulus.bin_shifts         = [0:bpf:(GLMPars.stimfilter.frames-1)*bpf];
            linearfilters.Stimulus.note1              = 'Filter is in [x,y,"frames before current bin"]';
            linearfilters.Stimulus.note2              = 'Recall each bin is housed in a frame (multiple bins per frame';
            linearfilters.Stimulus.note3              = 'frame_shifts describes the transfrom from time index to frames ahead of current bin';
        end  
    end
    
    if i_fit == 1, linearfilters_linstim    = linearfilters; rawfit_linstim    = rawfit; end
    if i_fit == 2, linearfilters_nonlinstim = linearfilters; rawfit_nonlinstim = rawfit; end 
    clear rawfit linearfilters
end
fittedGLM.linstim.rawfit               = rawfit_linstim ;
fittedGLM.linstim.linearfilters        = linearfilters_linstim;


fittedGLM.rawfit           = rawfit_nonlinstim ;
fittedGLM.linearfilters    = linearfilters_nonlinstim;



fittedGLM.note = 'in theory, linearfilters and t_bin/ binsperframe is sufficient for xval and simulation'; 
fittedGLM.fit_time = datestr(clock);
fittedGLM.writingcode = mfilename('fullpath');


end
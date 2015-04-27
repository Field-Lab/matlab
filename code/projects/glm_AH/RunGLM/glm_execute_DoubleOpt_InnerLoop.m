% Hacked code 2014-05-24
% Modulate the movie here and do a single optimization run
% we need a function to optimize .. this seems to be the cleanest way to do
% so

function [fstar,pstar] = glm_execute_DoubleOpt_InnerLoop(outer_optimizing_par,optimizing_type,glm_cellinfo, GLMType, GLMPars,fixed_covariates,home_spbins,fitmovie,inputstats,t_bin,p_init,troubleshoot)


[paramind] =  prep_paramindGP(GLMType, GLMPars); 
if strcmp(optimizing_type , 'input_pt_nonlinearity')
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        GLMPars.others.point_nonlinearity.increment_to_decrement = outer_optimizing_par;
    end
    if strcmp(GLMType.input_pt_nonlinearity_type , 'raisepower_meanafter')
        GLMPars.others.point_nonlinearity.scalar_raisedpower = outer_optimizing_par;
    end
    if strcmp(GLMType.input_pt_nonlinearity_type, 'oddfunc_powerraise_aboutmean')
        GLMPars.others.point_nonlinearity.scalar_raisedpower_aboutmean = outer_optimizing_par;
    end
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_shiftmean')
        GLMPars.others.point_nonlinearity.increment_to_decrement = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.shiftmean              = outer_optimizing_par(2);
    end
    if strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2')
        GLMPars.others.point_nonlinearity.coefficients.linear    = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.quadratic = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients.sqroot    = outer_optimizing_par(3);
    end
    if strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2_search2')
        GLMPars.others.point_nonlinearity.coefficients.linear    = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.quadoversqroot = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients;
    end
    if strcmp(GLMType.input_pt_nonlinearity_type, 'polynomial_androot_order2_search3')
        GLMPars.others.point_nonlinearity.coefficients.linear    = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.quadratic = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients.sqroot    = outer_optimizing_par(3);
    end
    if strcmp(GLMType.input_pt_nonlinearity_type,'polynomial_order3_part5')
        GLMPars.others.point_nonlinearity.coefficients.linear    = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.quadratic = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients.sqroot    = outer_optimizing_par(3);
        GLMPars.others.point_nonlinearity.coefficients.thirdpower= outer_optimizing_par(4);
        GLMPars.others.point_nonlinearity.coefficients.thirdroot = outer_optimizing_par(5);
    end
    
    if strcmp(GLMType.input_pt_nonlinearity_type,'polynomial_order5_part4')
        GLMPars.others.point_nonlinearity.coefficients.linear      = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.quadratic   = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients.sqroot      = outer_optimizing_par(3);
        GLMPars.others.point_nonlinearity.coefficients.thirdpower  = outer_optimizing_par(4);
        GLMPars.others.point_nonlinearity.coefficients.thirdroot   = outer_optimizing_par(5);
        GLMPars.others.point_nonlinearity.coefficients.fourthpower = outer_optimizing_par(6);
        GLMPars.others.point_nonlinearity.coefficients.fourthroot  = outer_optimizing_par(7);
        GLMPars.others.point_nonlinearity.coefficients.fifthpower  = outer_optimizing_par(8);
        GLMPars.others.point_nonlinearity.coefficients.fifthroot   = outer_optimizing_par(9);
    end
    
    
    if strcmp(GLMType.input_pt_nonlinearity_type,'piecelinear_fourpiece_eightlevels')
        GLMPars.others.point_nonlinearity.coefficients.slope_quartile_1   = outer_optimizing_par(1);
        GLMPars.others.point_nonlinearity.coefficients.slope_quartile_2   = outer_optimizing_par(2);
        GLMPars.others.point_nonlinearity.coefficients.slope_quartile_3   = outer_optimizing_par(3);
        GLMPars.others.point_nonlinearity.coefficients.slope_quartile_4   = outer_optimizing_par(4);
    end
    
    
    
end
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);


if exist('troubleshoot','var') && troubleshoot.doit
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA,troubleshoot);
else
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
end


% Load up lcif_nonlinearity if we are looking at a postfilter_nonlinearity
if strcmp(optimizing_type , 'postfilter_nonlinearity')
    lcif_nonlinearity.type = GLMType.postfilter_nonlinearity_type;
    if strcmp(GLMType.postfilter_nonlinearity_type ,  'piece_linear_aboutmean') 
        lcif_nonlinearity.increment_to_decrement = outer_optimizing_par;
        lcif_nonlinearity.filter_index = paramind.X;
        
        allparams = 1:paramind.paramcount;
        lcif_nonlinearity.linear_index = setdiff(allparams, lcif_nonlinearity.filter_index);
    end
    if strcmp(GLMType.postfilter_nonlinearity_type ,  'oddfunc_powerraise_aboutmean')   
            lcif_nonlinearity.scalar_raisedpower = outer_optimizing_par;
            lcif_nonlinearity.filter_index = paramind.X;
            allparams = 1:paramind.paramcount;
            lcif_nonlinearity.linear_index = setdiff(allparams, lcif_nonlinearity.filter_index);

	end
   GLMType.lcif_nonlinearity = lcif_nonlinearity;
end






bins = size(X_bin,2);
clear WN_STA center_coord



%%% Normally 'display' is set to  'iter' %%%
%%% but turned it off this time, case double opt %%
optim_struct = optimset(...
   'derivativecheck','off',...
   'diagnostics','off',...  % 
   'display','off',...  %'iter-detailed',.. .   %%   
   'funvalcheck','off',... % don't turn this on for 'raw' condition (edoi).
   'GradObj','on',...
   'largescale','on',...
   'Hessian','on',...
   'MaxIter',GLMPars.optimization.maxiter,... % you may want to change this
   'TolFun',10^(-(GLMPars.optimization.tolfun)),...
   'TolX',10^(-(GLMPars.optimization.tolx)));

if GLMType.CONVEX
        glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
        % Maybe move this inside of the stimulus preparation % 
        bpf         = GLMPars.bins_per_frame;
        shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
        X_bin_shift = prep_timeshift(X_bin,shifts);
        if isfield(paramind, 'MU')
            glm_covariate_vec( paramind.MU , : ) = fixed_covariates.MU_bin;
        end
        if isfield(paramind, 'X')
            glm_covariate_vec( paramind.X , : ) = X_bin_shift;
        end
        if isfield(paramind, 'PS')
            glm_covariate_vec( paramind.PS , : ) = fixed_covariates.PS_bin;
        end
        if isfield(paramind, 'CP')
            glm_covariate_vec( paramind.CP , : ) = fixed_covariates.CP_bin;
        end
        
        if isfield(GLMType, 'lcif_nonlinearity')
            [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction_withNL(p,glm_covariate_vec,home_spbins,t_bin, GLMType.lcif_nonlinearity),p_init,optim_struct);
        else
            [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),p_init,optim_struct);
        end
end
if ~GLMType.CONVEX
        bpf               = GLMPars.bins_per_frame;
        frame_shifts      = 0:1:(GLMPars.stimfilter.frames-1);
        % dnote part that is convex
        convex_cov = NaN(paramind.convParams , bins );  % convex covariates
        if isfield(paramind, 'MU')
            convex_cov( paramind.MU , : ) = fixed_covariates.MU_bin;
        end
        if isfield(paramind, 'PS')
            convex_cov( paramind.PS , : ) = fixed_covariates.PS_bin;
        end
        if isfield(paramind, 'CP')
            convex_cov( paramind.CP , : ) = fixed_covariates.CP_bin;
        end
        filtertype = GLMType.stimfilter_mode;
         %  glm_nonconvex_optimizationfunction(p_init,GLMType.stimfilter_mode,nonconvex_paramind,convex_covariate_vec,X_frame,frame_shifts, bpf, home_spbins,t_bin)
        if isfield(GLMType, 'lcif_nonlinearity')
            [pstar fstar eflag output] = fminunc(@(p) glm_nonconvex_optimizationfunction(p,...
            filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins,t_bin, GLMType.lcif_nonlinearity),p_init,optim_struct);
        else
        [pstar fstar eflag output] = fminunc(@(p) glm_nonconvex_optimizationfunction(p,...
            filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins,t_bin),p_init,optim_struct);
        end
end
end
    
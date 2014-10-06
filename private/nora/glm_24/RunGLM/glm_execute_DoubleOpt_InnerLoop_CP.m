% Hacked code 2014-05-24
% Modulate the movie here and do a single optimization run
% we need a function to optimize .. this seems to be the cleanest way to do
% so

function [fstar,pstar] = glm_execute_DoubleOpt_InnerLoop_CP(outer_optimizing_par,optimizing_type,glm_cellinfo, GLMType, GLMPars,fixed_covariates,home_spbins,fitmovie,t_bin,troubleshoot)

% NBCoupling
[paramind] =  prep_paramindGP_CP(GLMType, GLMPars, length(glm_cellinfo.pairs)); 
p_init     =  .01* ones(paramind.paramcount,1);
if strcmp(optimizing_type , 'input_pt_nonlinearity')
    if strcmp(GLMType.input_pt_nonlinearity_type, 'piece_linear_aboutmean')
        GLMPars.others.point_nonlinearity.increment_to_decrement = outer_optimizing_par;
    end
    if strcmp(GLMType.input_pt_nonlinearity_type , 'raisepower_meanafter')
        GLMPars.others.point_nonlinearity.scalar_raisedpower = outer_optimizing_par;
    end
end
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);


if exist('troubleshoot','var') && troubleshoot.doit
[X_frame,X_bin]    = prep_stimcelldependentGP(GLMType, GLMPars, fitmovie, center_coord, WN_STA,troubleshoot);
else
[X_frame,X_bin]    = prep_stimcelldependentGP(GLMType, GLMPars, fitmovie, center_coord, WN_STA);
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
    % NBCoupling 05-28-14
    if isfield(paramind, 'CP')
        for j_pair=1:length(glm_cellinfo.pairs)
            glm_covariate_vec( paramind.CP{j_pair} , : ) = fixed_covariates.CP_bin{j_pair};
        end
    end
    % end NBCoupling
        [pstar fstar eflag output]     = fminunc(@(p) glm_convex_optimizationfunction(p,glm_covariate_vec,home_spbins,t_bin),p_init,optim_struct);
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
    % NBCoupling 05-28-14
    if isfield(paramind, 'CP')
        for j_pair=1:length(glm_cellinfo.pairs)
           convex_cov( paramind.CP{j_pair} , : ) = fixed_covariates.CP_bin{j_pair};
        end
    end
    % end NBCoupling
        filtertype = GLMType.stimfilter_mode;
         %  glm_nonconvex_optimizationfunction(p_init,GLMType.stimfilter_mode,nonconvex_paramind,convex_covariate_vec,X_frame,frame_shifts, bpf, home_spbins,t_bin)
        [pstar fstar eflag output] = fminunc(@(p) glm_nonconvex_optimizationfunction(p,...
            filtertype,paramind,convex_cov,X_frame,frame_shifts, bpf, home_spbins,t_bin),p_init,optim_struct);
end
end
    
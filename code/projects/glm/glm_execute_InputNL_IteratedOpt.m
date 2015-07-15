function [fittedGLM] = glm_execute_InputNL_IteratedOpt(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,optional_arg)
% Part 1: Fit GLM pre input NL
% Part 2: Iterate through NL fits

glm_cellinfo0 = glm_cellinfo;


GLMPars = GLMParams;
if strcmp(GLMType.input_pt_nonlinearity_type, 'log_powerraise')
    % Run with non-linearity parameter set to identity to get first
    % estimate of filter
    NL_Par_0 = 0;
    
    search_min = -1;
    search_max = 1;
    
    GLMPars.others.point_nonlinearity.log_powerraise = NL_Par_0;
    options{1}.name = 'GLMPars';
    options{1}.GLMPars = GLMPars;
    
    glm_cellinfo.cell_savename = sprintf('%s_preNL', glm_cellinfo0.cell_savename);
    glm_cellinfo.d_save        = sprintf('%s/earlierfits', glm_cellinfo0.d_save);
    if ~exist(glm_cellinfo.dsave,'dir'), mkdir(glm_cellinfo.d_save); end
    
    fittedGLM = glm_execute(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,options);
    
    
    %%% test  this should return same value as fittedGLM.rawfit.objective_val
    %{
    t_bin = fittedGLM.t_bin;
    bins  = fittedGLM.bins_per_frame * size(fitmovie,3);
    pstar = fittedGLM.rawfit.opt_params;
    [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins);
    objval = subR_modinputNL_findobj(NL_Par_0, lcif_nonstim.total, pstar, ...
        GLMType, GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins);
    %}
end


loops = 2;
for i_loop = 1:loops
    t_bin = fittedGLM.t_bin;
    bins  = fittedGLM.bins_per_frame * size(fitmovie,3);
    pstar = fittedGLM.rawfit.opt_params;
    [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins);

    optim_struct = optimset(...
       'derivativecheck','off',...
       'diagnostics','off',...  % 
       'display','iter',...  %'iter-detailed',... 
       'MaxFunEvals',50,... % you may want to change this
       'TolX',.01) ;
    [NL_Par_star fstar eflag] = fminbnd(@(NL_Par) ...
                subR_modinputNL_findobj(NL_Par,lcif_nonstim.total, pstar,...
                GLMType, GLMPars, fitspikes, fitmovie, inputstats,...
                glm_cellinfo, t_bin,bins),...
                search_min,search_max, optim_struct);
    
    % search over smaller subdomains        
    dist_1 = abs(NL_Par_star - search_min);
    dist_2 = abs(NL_Par_star - search_max);
    search_min = NL_Par_star - .5*min(dist_1,dist_2)
    search_max = NL_Par_star + .5*min(dist_1,dist_2)
    
            
    GLMPars.others.point_nonlinearity.log_powerraise = NL_Par_star;
    options{1}.name    = 'GLMPars';
    options{1}.GLMPars = GLMPars;
    options{2}.name    = 'p_init';
    options{2}.p_init  = pstar;
    glm_cellinfo.cell_savename = sprintf('%s_afterround_%d', glm_cellinfo0.cell_savename,i_loop);
    glm_cellinfo.d_save        = sprintf('%s/earlierfits', glm_cellinfo0.d_save);
    if ~exist(glm_cellinfo.dsave,'dir'), mkdir(glm_cellinfo.d_save); end
    
    if i_loop == loops
        glm_cellinfo = glm_cellinfo0;
    end
    fittedGLM = glm_execute(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,options);
end



end


function obj_val = subR_modinputNL_findobj(NL_Params, lcif_external, pstar, GLMType, GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins)
% AKHEITMAN 2015-07-15
% subRoutine which will get optimized for inding input NL


% Part 1: Unpack non-linearity param
% Part 2: Find Stim driven lcif (with input
% Part 3: Add in external lcif and find obj_val 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Unpack non-linearity param
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(GLMType.input_pt_nonlinearity_type, 'log_powerraise')
    GLMPars.others.point_nonlinearity.log_powerraise = NL_Params;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Find Stim driven lcif (with input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord
if GLMType.CONVEX
    glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    pstar = pstar';
    lcif_stim = pstar(paramind.X) *X_bin_shift;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: Add in external lcif and find obj_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );
total_lcif = lcif_stim + lcif_external;
cif        = exp(total_lcif);
obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));

end
function [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins)
% AKHEITMAN 2015-07-15  extract non-stim components of lcif


% Construct PS_Basis
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
    
    if isfield(GLMType, 'PS_constrain')
        if strcmp(GLMType.PS_constrain.type, 'PS_netinhibitory_domainconstrain_COB')
            ps_basis_0 = ps_basis; clear ps_basis
            v        = sum(ps_basis_0,1);
            v        = v / norm(v) ;
            orthog_v = null(v);
            COB      = [v', orthog_v] ;
            ps_basis = (inv(COB) * ps_basis_0')' ;
        end
    end
            
end

% Find PS_bin and MU_bin (pre-param multiply)
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


% Multiply by param for final log-cif values
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
pstar = pstar';
total_lcif = 0;
if isfield(paramind, 'MU')
    lcif_nonstim.components.mu   = pstar(paramind.MU)*MU_bin;
    total_lcif = total_lcif + lcif_nonstim.components.mu;
end
if isfield(paramind, 'PS')
    lcif_nonstim.components.ps   = pstar(paramind.PS)*PS_bin;
    total_lcif = total_lcif +  lcif_nonstim.components.ps;
end
lcif_nonstim.total = total_lcif;



end




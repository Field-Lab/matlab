% AKHEITMAN 2015-07-19  started.. 
% First version up and working
function [fittedGLM] = glm_execute_OutputNL(postfilter_NL,preoutputNL_fittedGLM, GLMType,...
    fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,optional_arg)

% Part 1: Extract lcif components
% Part 2: 

%%%%%%%%%%%%
% Part 1 %%%
%%%%%%%%%%%%
fittedGLM= preoutputNL_fittedGLM;
t_bin = fittedGLM.t_bin;
bins  = fittedGLM.bins_per_frame * size(fitmovie,3);
pstar = fittedGLM.rawfit.opt_params;
[lcif_fit.nonstim] = subR_lcif_nonstim(pstar, fittedGLM.GLMType,fittedGLM.GLMPars,fitspikes,t_bin,bins);
[lcif_fit.stim.preNL objval] = subR_findobj_lcifstim(lcif_fit.nonstim.total, pstar, ...
    fittedGLM.GLMType, fittedGLM.GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins);  

[lcif_test.stim.preNL ] = subR_findobj_lcifstim([], pstar, ...
    fittedGLM.GLMType, fittedGLM.GLMPars, [], testmovie, inputstats, glm_cellinfo, t_bin,bins); 


t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );

display(sprintf('objective value from preNL fitted GLM::  %d',fittedGLM.rawfit.objective_val))
display(sprintf('objective value subR should be ~ equal::  %d',objval));              



%%%%%%%%%%%%
% Part 2 %%%
%%%%%%%%%%%%
NL_Input.input_fit  = lcif_fit.stim.preNL  / std(lcif_fit.stim.preNL);
NL_Input.input_test = lcif_test.stim.preNL / std(lcif_fit.stim.preNL);

if strcmp(postfilter_NL.type, 'Logistic_2Par_fixMU')    
    NL_Input.y_int              = exp( lcif_fit.nonstim.components.mu(1) );
    NL_Input.rawfilteroutput    = lcif_fit.stim.preNL;
    NL_Input.scale_rawtoNLinput = std(lcif_fit.stim.preNL);
    
    if isfield(lcif_fit.nonstim.components, 'ps')
        NL_Input.fit_lcif_additionaldrive = lcif_fit.nonstim.components.ps;
    else
        NL_Input.fit_lcif_additionaldrive = 0;
    end
    NL_Output = subR_LogFixMu_fmincon(NL_Input,home_spbins,t_bin);
    
    
end

end


function NLOutput = subR_LogFixMu_fmincon(NL_Input,home_spbins,t_bin)
% 2015-06-24  AKHeitman

% NL_Input needs: y_int (value with no drive)
%                input_test(corresponding input to the NL for fit data)
%                input_fit(corresponding input to the NL for test data) 

optim_struct = optimset(...
                    'derivativecheck','off','diagnostics','off',...  % 
                    'display','iter','funvalcheck','off',... 
                    'MaxIter',100,'TolFun',10^(-6),'TolX',10^(-9) );
lowerbound = [NL_Input.y_int+1 .1];
upperbound = [1000         100];
LOGI_Params0 = [100,1];
[LOGI_Params_Opt, new_objval, eflag, output] = fmincon(@(LOGI_Params) subR_objval_LOGISTIC...
    (LOGI_Params, NL_Input.y_int,NL_Input.input_fit, NL_Input.fit_lcif_additionaldrive, home_spbins,t_bin),...
    LOGI_Params0,[],[],[],[],lowerbound,upperbound,[],optim_struct);

[~, lcif_LOGI_crossvaltest] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_test, 0, [],t_bin); 

[~, lcif_LOGI_fit] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin); 


NLOutput.note_metric = 'Logistic with Null Rate fixed to tonic drive. Driven by linear filter output,normalized to std 1';
NLOutput.param_string = sprintf('2 Params fit with fmincon,  Optimal Max Rate: %1.1e, Optimal Slope: %1.1e',...
    LOGI_Params_Opt(1), LOGI_Params_Opt(2) );
NLOutput.maxrate = LOGI_Params_Opt(1);
NLOutput.slope   = LOGI_Params_Opt(2);
NLOutput.fmincon_output = output;
NLOutput.eflag = eflag;
NLOutput.crossvaltest_finalrate = exp(lcif_LOGI_crossvaltest);
NLOutput.new_objval = new_objval;
NLOutput.fit_rate = exp(lcif_LOGI_fit);
NLOutput.codename = mfilename('fullpath');
clear lcif_LOGI_crossvaltest dummy lowerbound upperbound LOGI_Params0
% Check Code
%{
lcif   = baseGLM.lcif_fit.stim + baseGLM.lcif_fit.mu; 
cif    = exp(lcif);
objval = -( sum( lcif(home_spbins) ) - t_bin * sum(cif) );
pars = [10,1];
objval = objval_LOGISTIC(pars,NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin)
%}
end
function [objval lcif_LOGI] = subR_objval_LOGISTIC(LOGI_PARAMS, Y_INT, lcif_intoLOGI, lcif_ext, spikebins,t_bin)
MAX     = LOGI_PARAMS(1);
RATE    = LOGI_PARAMS(2);


OFFSET  = log( (MAX/Y_INT) - 1  ) / RATE;                    

lcif_LOGI = log(MAX ./ (1 + exp(-RATE * (lcif_intoLOGI- OFFSET) )));

lcif   = lcif_LOGI + lcif_ext;
cif    = exp(lcif);
objval = -( sum( lcif(spikebins) ) - t_bin * sum(cif) );
end
function [lcif_stim, obj_val] = subR_findobj_lcifstim(lcif_external, pstar, GLMType, GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins)
% AKHEITMAN 2015-07-15
% subRoutine which will get optimized for inding input NL

% Part 1: Find Stim driven lcif (with input
% Part 3: Add in external lcif and find obj_val 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Find Stim driven lcif (with input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[paramind] =  prep_paramindGP(GLMType, GLMPars); 
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord



if GLMType.CONVEX
   
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    pstar = pstar';
    lcif_stim = pstar(paramind.X) *X_bin_shift;
end


if ~GLMType.CONVEX
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    pstar = pstar;
    
    
    frame_shifts = 0:1:(GLMPars.stimfilter.frames-1);
    if isfield(GLMPars.stimfilter,'frames_negative')
        frame_shifts = -(GLMPars.stimfilter.frames_negative):1:(GLMPars.stimfilter.frames-1)*1;
    end
    pixels  = size(X_frame,1);
    frames  = size(X_frame,2);
    
    
    % PARAMS TO GET INCORPORTED IN COVARIATE VEC 
    TimeFilter  = pstar(paramind.time1);
    SpaceFilter = pstar(paramind.space1);
    
    % FIND SPATIAL FILTER COVARIATE VEC (USING TIMEFILTER)
    if min(frame_shifts) == 0;
        convolvingFilter = (TimeFilter);        
    elseif min(frame_shifts) > 0
        padzeros         = zeros(min(frame_shifts),1);
        convolvingFilter = [padzeros ; (TimeFilter)]; 
    else
        error('frame_shifts should be >=0')
    end
    timeconvStim = zeros( pixels , (frames+length(convolvingFilter)-1) );
    for i_row = 1:size(X_frame,1)
        timeconvStim(i_row,:) = conv(X_frame(i_row,:) , convolvingFilter);
    end
    timeconvStim = timeconvStim(:,1:frames);
    bins   = bpf * frames;
    A      = repmat(timeconvStim, [ bpf,1]); 
    spatial_covariatevec  = reshape(A, [pixels, bins]);
    
    % FIND TEMPORAL FILTER COVARIATE VEC (USING SPATIAL FILTER)
    spaceconvStim           = SpaceFilter' * X_frame;
    B                       = repmat(spaceconvStim, [ bpf,1]); 
    spaceconvStim_bin       = reshape(B, [1, bins]);
    
    bin_shifts              = bpf *frame_shifts;
    temporal_covariatevec   = prep_timeshift(spaceconvStim_bin,bin_shifts);
    
    
    % STIMULUS ADDED TWICE / DIVIDE EACH COVARIATE VEC BY TWO
    stim_covariate_vec =  [.5*spatial_covariatevec; .5*temporal_covariatevec];    
    lcif_stim = (pstar(paramind.X))'  * stim_covariate_vec;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Add in external lcif and find obj_val
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout > 1
    t_bin        = t_bin;
    home_sptimes = fitspikes.home';
    home_spbins  = ceil(home_sptimes / t_bin);
    home_spbins  = home_spbins(find(home_spbins < bins) );
    total_lcif = lcif_stim + lcif_external;
    cif        = exp(total_lcif);
    obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));
end

end
function [lcif_nonstim] = subR_lcif_nonstim(pstar, GLMType,GLMPars,fitspikes,t_bin,bins)
% AKHEITMAN 2015-07-15  extract non-stim components of lcif


% Construct PS_Basis
if GLMType.PostSpikeFilter
    basis_params  = GLMPars.spikefilters.ps;
    ps_basis      = prep_spikefilterbasisGP(basis_params,t_bin);
    

    if isfield(GLMType, 'special_arg') && isfield(GLMType.special_arg,'PS_Constrain')
            ps_basis_0 = ps_basis; clear ps_basis
            v        = sum(ps_basis_0,1);
            v        = v / norm(v) ;
            orthog_v = null(v);
            COB      = [v', orthog_v] ;
            ps_basis = (inv(COB) * ps_basis_0')' ;
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




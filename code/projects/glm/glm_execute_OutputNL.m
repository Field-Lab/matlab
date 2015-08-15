% AKHEITMAN 2015-07-19  started.. 
% First version up and working
function [fittedGLM] = glm_execute_OutputNL(postfilter_NL,fittedGLM_preoutputNL, GLMType,...
    fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo,neighborspikes,optional_arg)

% Part 1: extract original lcif components
% Part 2: optimize postfilter-NL
% Part 3: 


%%%%%%%%%%%%
% Part 1: extract original lcif components
%%%%%%%%%%%%
fittedGLM= fittedGLM_preoutputNL;
t_bin = fittedGLM.t_bin;
bins  = fittedGLM.bins_per_frame * size(fitmovie,3);
pstar = fittedGLM.rawfit.opt_params;
[lcif_fit.nonstim] = subR_lcif_nonstim(pstar, fittedGLM.GLMType,fittedGLM.GLMPars,fitspikes,t_bin,bins);
[lcif_fit.stim.preNL objval] = subR_findobj_lcifstim(lcif_fit.nonstim.total, pstar, ...
    fittedGLM.GLMType, fittedGLM.GLMPars, fitspikes, fitmovie, inputstats, glm_cellinfo, t_bin,bins);  
[lcif_test.stim.preNL ] = subR_findobj_lcifstim([], pstar, ...
    fittedGLM.GLMType, fittedGLM.GLMPars, [], testmovie, inputstats, glm_cellinfo, t_bin,bins); 
display(sprintf('objective value from preNL fitted GLM::  %d',fittedGLM.rawfit.objective_val))
display(sprintf('objective value subR should be ~ equal::  %d',objval));              
clear fittedGLM


%%%%%%%%%%%%
% Part 2: optimize postfilter-NL 
%%%%%%%%%%%%
t_bin        = t_bin;
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );

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


%%%%%%%%%%%%
% Part 3: cross validated performance
%%%%%%%%%%%%
recorded_raster = fittedGLM_preoutputNL.xvalperformance.rasters.recorded;
if GLMType.PostSpikeFilter
    PS = fittedGLM_preoutputNL.linearfilters.PostSpike.Filter;
    NL_xvalperformance = subR_xvalperformance_givenrate(NL_Output.crossvaltest_finalrate, recorded_raster, t_bin,PS); 
else
    NL_xvalperformance = subR_xvalperformance_givenrate(NL_Output.crossvaltest_finalrate, recorded_raster, t_bin); 
end


%%%%%%%%%%%%
% Part 4: save and print
%%%%%%%%%%%%
% Hack to Spit Out                
fittedGLM.cell_savename        = fittedGLM_preoutputNL.cell_savename;
fittedGLM.cellinfo             = fittedGLM_preoutputNL.cellinfo;
fittedGLM.postfilter_nonlinearity_type         = postfilter_NL.type;
fittedGLM.postfilter_nonlinearity              = NL_Output;
fittedGLM.xvalperformance                      = NL_xvalperformance;
fittedGLM.t_bin                = fittedGLM_preoutputNL.t_bin;
fittedGLM.bins_per_frame       = fittedGLM_preoutputNL.bins_per_frame;
fittedGLM.rawfit.objective_val = NL_Output.new_objval;
fittedGLM.fit_time    = datestr(clock);
fittedGLM.writingcode =  mfilename('fullpath');
fittedGLM.linearfilters.Stimulus         = fittedGLM_preoutputNL.linearfilters.Stimulus;
fittedGLM.linearfilters.Stimulus_rescale = NL_Input.scale_rawtoNLinput;
fittedGLM.linearfilters.Stimulus_rescalenote = ...
    'divide stimfilter by rescaler (to z-score to 1,convolve with stim, to  retrieve pre-Nonlinearity rate';
fittedGLM.linearfilters.Stimulus_note = 'rescaled to have z score of 1, preceeding use of non-linearity'
if GLMType.PostSpikeFilter
    fittedGLM.linearfilters.PostSpike = fittedGLM_preoutputNL.linearfilters.PostSpike;
end
fittedGLM.d_save = glm_cellinfo.d_save;
savename = sprintf('%s/%s',glm_cellinfo.d_save, glm_cellinfo.cell_savename);
eval(sprintf('save %s.mat fittedGLM', savename)); 


mu = lcif_fit.nonstim.components.mu(1);
stimtransform.normalized_filteroutput_fit= NL_Input.input_fit;
stimtransform.cif_withNL_fit             = NL_Output.fit_rate;
stimtransform.cif_rawGLM_fit             = exp( std(lcif_fit.stim.preNL) * NL_Input.input_fit  + mu );
stimtransform.cif_withNL_test            = NL_Output.crossvaltest_finalrate;
stimtransform.cif_rawGLM_test            = exp( std(lcif_fit.stim.preNL) * NL_Input.input_test + mu );
stimtransform.inputNL                    = NL_Input;


subR_plotfittedNL(fittedGLM, fittedGLM_preoutputNL, stimtransform,glm_cellinfo.d_save)
end
function subR_plotfittedNL(fittedGLM, fittedGLM_preoutputNL, stimtransform, savedir)

% Cleaned up AKHeitman 2015-06-24
homedir = pwd;
clf;  
printname = sprintf('DiagNLPlot_%s',fittedGLM.cell_savename);
info    = fittedGLM.cellinfo;
GLMType = fittedGLM_preoutputNL.GLMType;

% text box
subplot(5,1,1)
axis off
set(gca, 'fontsize', 12)
obj_NEW = fittedGLM.rawfit.objective_val;
obj_OLD = fittedGLM_preoutputNL.rawfit.objective_val;
optNL_describe  = fittedGLM.postfilter_nonlinearity.note_metric;
optNL_string    = fittedGLM.postfilter_nonlinearity.param_string;

c = 0; offset = 0; delta_c = 1.1;
text(-offset, 1-0.1*c,sprintf('%s: %s %d: %s-Fit (red): POSTNL refit with %s',...
    info.exp_nm, info.celltype,info.cid, GLMType.fit_type, fittedGLM.postfilter_nonlinearity_type), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Red is original GLM, Blue is Postfilter Nonlinearity'))
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Objval PctChange %d: from %1.2e to %1.2e',...
   round(100*(obj_NEW-obj_OLD)/obj_OLD),obj_OLD,obj_NEW), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_describe), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('%s', optNL_string), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Original Fit: %s',GLMType.fitname), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('NL fit Computated at %s',datestr(clock)), 'interpreter','none')
c = c + delta_c; 
text(-offset, 1-0.1*c,sprintf('Mfile: %s', mfilename('fullpath')), 'interpreter','none' );



% plotting non-linearities
x1 = sort(stimtransform.normalized_filteroutput_fit); 
y1 = sort(stimtransform.cif_rawGLM_fit);
y2 = sort(stimtransform.cif_withNL_fit);

LW = 2;
subplot(5,3,[4 7]); set(gca, 'fontsize', 8); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([-4,4]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity - Central Portion')

subplot(5,3,[5 8]); set(gca, 'fontsize', 8); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([ max(min(x1),-10),0]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Inhibitory Portion')

subplot(5,3,[6 9]); set(gca, 'fontsize', 8); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([0,min(10,max(x1))]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Excitatory Portion')


% plot rasters
subplot(5,1,[4 5]); set (gca, 'fontsize',10)
secs     = 6;
dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
NL_rast  = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([0 , 2*trials]); hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    NL1  = time(find(NL_rast(i_trial,:)));
    % Plot the raster
    plot(rec1, i_trial, 'k.')    
    yshift = i_trial;
    if length(NL1) < 4*length(rec1) 
        if length(NL1) > 0
            plot(NL1, yshift + trials, 'b.')
        end
    end
end
xlabel('seconds'); ylabel('trials')
cd(savedir)
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
cd(homedir)
end

function NL_xvalperformance     = subR_xvalperformance_givenrate( stimdrivenrate, logicalspike, t_bin,PS)
% AKHEITMAN 2015-06-24  it works!
% AKH 2015-07-19 integrated into glm main code
% New version includes PS filter

params.bindur     = t_bin;
params.bins       = length(stimdrivenrate);
params.trials     = size(logicalspike,1);
params.testdur_seconds = params.bindur * params.bins ;   

% Set log-conditional as stim driven only
lcif_teststim = log(stimdrivenrate);
lcif = repmat(lcif_teststim, params.trials,1);  


if exist('PS','var')
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif + lcif_ps;
end


glm_ratepersec  = exp(lcif);
glm_rateperbin  = params.bindur * glm_ratepersec;

spikerate_bin    = size(find(logicalspike(:))) /  size(logicalspike(:));      
model_null0      = spikerate_bin * ones(1, params.bins);
model_null       = repmat(model_null0, params.trials, 1);
null_logprob     = sum(eval_rasterlogprob(logicalspike, model_null, 'binary', 'conditioned'));
[raster_logprob_bin] = eval_rasterlogprob( logicalspike, glm_rateperbin,  'binary', 'conditioned') ;
glm_logprob       = sum(raster_logprob_bin);
glm_bits          = glm_logprob - null_logprob;
glm_bits_perspike = glm_bits / (sum(model_null0));
glm_bits_perbin   = glm_bits / params.bins;
glm_bits_persecond   = glm_bits / params.testdur_seconds;

NL_xvalperformance.note = 'Scores include optimized Non-Linearity';
NL_xvalperformance.logprob_null_raw            = null_logprob;
NL_xvalperformance.logprob_glm_raw      =  glm_logprob;
NL_xvalperformance.logprob_glm_bpspike  =  glm_bits_perspike;
NL_xvalperformance.logprob_glm_bpsec    =  glm_bits_persecond;

lcif_const  = lcif(1,:);
logical_sim = zeros(params.trials, params.bins);
% PS Filter 
if exist('PS','var')
    cif_psgain = exp(PS);
    ps_bins     = length(cif_psgain);
    for i_trial = 1 : size(logicalspike,1)
        cif0         = exp(lcif_const);         
        cif_ps       = cif0;
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins- ps_bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif_ps(i));
                cif_ps(i+1: i + ps_bins) =  cif_ps(i+1: i + ps_bins) .* (cif_psgain');
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
else
    for i_trial = 1 : size(logicalspike,1)
        cif         = exp(lcif_const);         
        binary_simulation = zeros(1,params.bins);
        for i = 1 : params.bins;
            roll = rand(1);
            if roll >  exp(-params.bindur*cif(i));
                binary_simulation(i)= 1;
            end
        end
        logical_sim(i_trial,:) = binary_simulation ;
    end
end
NL_xvalperformance.rasters.note           = 'glmsim includes altered non-linearity';
NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_sim        = logical_sim;
NL_xvalperformance.rasters.bintime        = params.bindur;
end       
function NLOutput = subR_LogFixMu_fmincon(NL_Input,home_spbins,t_bin)
% 2015-06-24  AKHeitman
% integrated into main GLM 2015-07-19

% Calls another subRoutine: subR_objval_LOGISTIC

% NL_Input needs: y_int (value with no drive)
%   input_test(corresponding input to the NL for fit data)
%   input_fit(corresponding input to the NL for test data) 
%   fit_lcif_additionaldrive (all terms outside of the non-linearity)
%       (ie. post-spike and coupling filters, eye-movements etc.)


% First Param: Max Rate  (top of the logistic curve)
% Second Param: 

% Step 1: Initialize Search
% Step 2: Search
% Step 3: Find Corresponding 


%%% Step 1: Initialize Search %%%%%
lowerbound = [NL_Input.y_int+1 .1];
upperbound = [1000         100];
LOGI_Params0 = [100,1];
optim_struct = optimset(...
                    'derivativecheck','off','diagnostics','off',...  % 
                    'display','iter','funvalcheck','off',... 
                    'MaxIter',100,'TolFun',10^(-6),'TolX',10^(-9) );
                
                
%%% Step 2: Search %%%
[LOGI_Params_Opt, new_objval, eflag, output] = fmincon(@(LOGI_Params) subR_objval_LOGISTIC...
    (LOGI_Params, NL_Input.y_int,NL_Input.input_fit, NL_Input.fit_lcif_additionaldrive, home_spbins,t_bin),...
    LOGI_Params0,[],[],[],[],lowerbound,upperbound,[],optim_struct);
[~, lcif_LOGI_crossvaltest] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_test, 0, [],t_bin); 
[~, lcif_LOGI_fit] =  subR_objval_LOGISTIC(LOGI_Params_Opt,...
    NL_Input.y_int,NL_Input.input_fit, 0, home_spbins,t_bin); 


NLOutput.note_metric = 'Logistic with Null Rate fixed to tonic drive. Driven by linear filter output,normalized to std 1';
NLOutput.param_string = sprintf('Nullresponse set at %1.1e, Fmincon fits: Optimal Max Rate: %1.1e, Optimal Slope: %1.1e',...
    NL_Input.y_int, LOGI_Params_Opt(1), LOGI_Params_Opt(2) );
NLOutput.maxrate       = LOGI_Params_Opt(1);
NLOutput.slope         = LOGI_Params_Opt(2);
NLOutput.null_response = NL_Input.y_int;
NLOutput.fmincon_output = output;
NLOutput.eflag = eflag;
NLOutput.crossvaltest_finalrate = exp(lcif_LOGI_crossvaltest);
NLOutput.new_objval = new_objval;
NLOutput.fit_rate = exp(lcif_LOGI_fit);
NLOutput.codename = mfilename('fullpath');
clear lcif_LOGI_crossvaltest dummy lowerbound upperbound LOGI_Params0








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




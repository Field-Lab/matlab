%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AKHeitman 2015-04-02

% Creates structure which dictates GLMType
% Loads cells 
% Loads stimuli / basic stimuli processing
% Loads spike trains / basic spike train processing
% Requires the organizedspikes structure with spike times relative
%    to start of each block of stimulus
% No direct GLM Paramater usage
% Feeds into glm_execute which is located in glm_core directory
% glm_execute along with glm_core 
%    which has no additional code dependencies, no loading of matfiles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wrap_bookeeping Calls %
%  NSEM_BaseDirectories
%  GLM_Settings
%  GLM_fitname
%  NSEM_secondaryDirectories
%  loadmoviematfiles
%  StimulusParams

% Main Call %
%   glm_execute  

% Subroutines at bottom of function
%  subR_concat_fitspikes_fromorganizedspikes
%  subR_createraster
%  subR_concat_fitmovie_fromblockedcell
%  subR_visionSTA_to_xymviCoord
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Unique Subroutine for this wrapper
%   subR_lcifdecomp_fittedGLM
%   subR_lcifstim_fittedGLM  (only stimulus portion .. good for testmovie)


% Version 4 Finds optimal LN (throws away all other terms)
% Version 3 Pulled out minimization functions.  Cleaned, but same saving.
%           2015-06-23
% Version 2 is last hack to work  (closed ~2015-06-16)


% Calling Sequences
%{

clear; clc
exps = [1 2 3 4]; stimtypes = [1 2]; celltypes = [1 2]; 
cell_subset = 'all'; postfilterNL.debug = false;
baseGLM.settings = {};
%baseGLM.settings{1}.type = 'PostSpikeFilter';
%baseGLM.settings{1}.name =  'OFF';
postfilterNL.type        = 'Logistic_fixMU_noPS';
%postfilterNL.type        = 'PieceLinear_FourPieces_2STD';
runoptions.print         = true;
glm_postfilterNL_LNonly(exps,stimtypes,celltypes,cell_subset,baseGLM.settings,postfilterNL,runoptions)

clear
exps = [1 2]; 
stimtypes = [1 2];
celltypes = [1 2];
cell_subset = 'all';
%cell_subset = 'glmconv_4pct';
postfilterNL.debug = false;
postfilterNL.type        = 'Logistic_fixMU_noPS';
baseGLM.settings{1}.type = 'cone_model';
baseGLM.settings{1}.name = 'rieke_linear'
baseGLM.settings{2}.type = 'input_pt_nonlinearity';
baseGLM.settings{2}.name= 'piecelinear_fourpiece_eightlevels';
runoptions.print         = true;
glm_postfilterNL_LNonly(exps,stimtypes,celltypes,cell_subset,baseGLM.settings,postfilterNL,runoptions)


clear
exps = [3]; 
stimtypes = [2 1];
celltypes = [2 1];
%cell_subset = 'all';
cell_subset = 'glmconv_4pct';
postfilterNL.debug = false;
postfilterNL.type        = 'Logistic_fixMU_noPS';
baseGLM.settings{1}.type = 'cone_model';
baseGLM.settings{1}.name = 'rieke_linear'
baseGLM.settings{2}.type = 'input_pt_nonlinearity';
baseGLM.settings{2}.name= 'piecelinear_fourpiece_eightlevels';
runoptions.print         = true;
glm_postfilterNL_LNonly(exps,stimtypes,celltypes,cell_subset,baseGLM.settings,postfilterNL,runoptions)

%}

function glm_postfilterNL_LNonly(exps,stimtypes,celltypes,cell_subset,baseGLM_settings,postfilterNL,runoptions)

% Load core directories and all eligible cells
BD = NSEM_BaseDirectories;
eval(sprintf('load %s/allcells.mat', BD.Cell_Selection));
% Define structure which uniquely defines GLM to be used 
if exist('baseGLM_settings', 'var')
    baseGLM.settings = baseGLM_settings; clear baseGLM_settings;
    baseGLM.Type = GLM_settings('default',baseGLM.settings);
else
    baseGLM.Type = GLM_settings('default');
end
baseGLM.Type.fitname    = GLM_fitname(baseGLM.Type); 
currentdir = pwd;

for i_exp = exps    
    for i_stimtype = stimtypes
        % Load master datarun, bookkeep
        exp_nm  = allcells{i_exp}.exp_nm;
        expname = allcells{i_exp}.expname;
        eval(sprintf('load %s/%s/datarun_master.mat', BD.BlockedSpikes,exp_nm));
        if i_stimtype == 1, stimtype = 'WN';   end
        if i_stimtype == 2, stimtype = 'NSEM'; end
        baseGLM.Type.fit_type = stimtype;        
        % Load and process stimulus
        [StimulusPars, exp_info] = StimulusParams(exp_nm, stimtype, baseGLM.Type.map_type);
        [blockedmoviecell, inputstats, origmatfile] = loadmoviematfile(exp_nm , stimtype, baseGLM.Type.cone_model,'fitmovie');
        [testmovie0]          = loadmoviematfile(exp_nm , stimtype, baseGLM.Type.cone_model,'testmovie');
        testmovie             = testmovie0{1}.matrix(:,:,StimulusPars.slv.testframes);
        baseGLM.Type.fitmoviefile  = origmatfile;
        if postfilterNL.debug
            display('shorten stimulus for post filter debugging mode')
            StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:2);
        end
        fitmovie_concat       = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , StimulusPars.slv); 
        % Directories  
        secondDir.exp_nm        = exp_nm; 
        secondDir.map_type      = baseGLM.Type.map_type; 
        secondDir.stim_type     = stimtype;
        secondDir.fitname       = baseGLM.Type.fitname;
        
        Dirs.WN_STAdir          = NSEM_secondaryDirectories('WN_STA', secondDir); 
        Dirs.organizedspikesdir = NSEM_secondaryDirectories('organizedspikes_dir', secondDir); 
        Dirs.baseglm            = NSEM_secondaryDirectories('loaddir_GLMfit', secondDir);
        
        
        % Hack to get the correct save directory  
        BD_hack = BD;
        if isfield(postfilterNL,'debug') && postfilterNL.debug
            BD_hack.GLM_output_raw = sprintf('%s/PostFilterNL/dbug_%s', BD.GLM_develop_output_raw,postfilterNL.type)
        else 
            BD_hack.GLM_output_raw = sprintf('%s/PostFilterNL/%s', BD.GLM_develop_output_raw,postfilterNL.type);
        end
        savedir  = NSEM_secondaryDirectories('savedir_GLMfit', secondDir,'',BD_hack)
        if ~exist(savedir,'dir'), mkdir(savedir); end
        
        
        
        % Loop through cells 
        for i_celltype = celltypes            
            if i_celltype == 1; cellgroup = allcells{i_exp}.ONP;  celltype = 'ONPar'; end
            if i_celltype == 2; cellgroup = allcells{i_exp}.OFFP; celltype = 'OFFPar'; end
            if strcmp(cell_subset,'all')
                candidate_cells = [allcells{i_exp}.ONP allcells{i_exp}.OFFP]
            elseif strcmp(cell_subset,'shortlist') || strcmp(cell_subset, 'debug') 
                [~,candidate_cells,~]  = cell_list(i_exp, cell_subset); 
                candidate_cells = cell2mat(candidate_cells) ; 
            elseif strcmp(cell_subset,'glmconv_4pct')
                eval(sprintf('load %s/allcells_glmconv.mat', BD.Cell_Selection));              
                conv_column = 2; 
                conv_index_ON = find(allcells_glmconv{i_exp}.ONP_CONV(:,conv_column));
                conv_index_OFF = find(allcells_glmconv{i_exp}.OFFP_CONV(:,conv_column));
                candidate_cells = [allcells{i_exp}.ONP(conv_index_ON) allcells{i_exp}.OFFP(conv_index_OFF)];
            end
            cellgroup = intersect(candidate_cells, cellgroup)
            for i_cell = 1:length(cellgroup)
                %% Actual Computation
                cid = cellgroup(i_cell); 
                cell_savename = sprintf('%s_%d', celltype,cid);
                display(sprintf('working on %s: %s', expname, cell_savename))
                eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', Dirs.organizedspikesdir, cell_savename));
                % Process spikes for glm_execute with proper subroutines
                fitspikes_concat.home  = subR_concat_fitspikes_fromorganizedspikes(organizedspikes.block, StimulusPars.slv);
                testspikes_raster.home = subR_createraster(organizedspikes.block, StimulusPars.slv);
                % load fittedGLM
                eval(sprintf('load %s/%s.mat fittedGLM', Dirs.baseglm, cell_savename));
                glm_cellinfo = fittedGLM.cellinfo;
                
                
                
                % works
                [baseGLM.lcif_fit,baseGLM.objval] =  subR_lcifdecomp_fittedGLM(fittedGLM.rawfit.opt_params,...
                    fittedGLM.GLMType,fittedGLM.GLMPars,fitspikes_concat,fitmovie_concat,inputstats,glm_cellinfo);
                [baseGLM.lcif_test.stim] = subR_lcifstim_fittedGLM(fittedGLM.rawfit.opt_params,...
                    fittedGLM.GLMType,fittedGLM.GLMPars,testmovie,inputstats,glm_cellinfo);
                
                % Necessary for refitting
                t_bin        = fittedGLM.t_bin;
                home_sptimes = fitspikes_concat.home';
                home_spbins  = ceil(home_sptimes / t_bin);
                home_spbins  = home_spbins(find(home_spbins < length(baseGLM.lcif_fit.stim)) );
                
                
                % normalized output of the stimulus filter
                NL_Input.input_fit  = baseGLM.lcif_fit.stim  / std(baseGLM.lcif_fit.stim);
                NL_Input.input_test = baseGLM.lcif_test.stim / std(baseGLM.lcif_fit.stim);
                NL_Input.y_int      = exp( baseGLM.lcif_fit.mu(1) );
                NL_input.rawfilteroutput    = baseGLM.lcif_fit.stim;
                NL_input.scale_rawtoNLinput = std(baseGLM.lcif_fit.stim);
                
                
                % CAREFULLY CONSTRAINED LOCAL SEARCH
                if fittedGLM.GLMType.PostSpikeFilter
                    lcif_stim0 = baseGLM.lcif_fit.stim;                    
                    SCALARS_INIT = [ baseGLM.lcif_fit.mu(1) , 1];
                    lowerbound = (1/10) * SCALARS_INIT;
                    upperbound = 10 * SCALARS_INIT;

                    optim_struct = optimset(...
                    'derivativecheck','off','diagnostics','off',...  % 
                    'display','off','funvalcheck','off',... 
                    'MaxIter',100,'TolFun',10^(-6),'TolX',10^(-9) );
                    
                    display('Running refit of tonic drive without the PS Filter')
                    [SCALARS_OPT, new_objval, eflag, output] = fmincon(@(SCALARS) subR_rescale_stim...
                    (SCALARS, lcif_stim0, home_spbins,fittedGLM.t_bin),...
                    SCALARS_INIT,[],[],[],[],lowerbound,upperbound,[],optim_struct);
                    display(sprintf('Tonic Drive with PS %1.2e hz, Without PS %1.2e hz',...
                        exp( baseGLM.lcif_fit.mu(1) ),   exp(SCALARS_OPT(1))));
                    display(sprintf('Without PS, Stim rescaled a factor %1.2e', SCALARS_OPT(2)))
                    
                    NL_Input.y_int      = exp(SCALARS_OPT(1));
                    clear SCALARS_INIT upperbound lowerbound new_objval eflag output lcif_stim0
                end

                if strcmp(postfilterNL.type, 'Logistic_fixMU_noPS')
                    NL_Output = LogFixMu_fmincon(NL_Input,home_spbins,t_bin);
                end
                
                fittedGLM_preNL = fittedGLM; 
                recorded_raster = fittedGLM_preNL.xvalperformance.rasters.recorded;
                NL_xvalperformance = subR_xvalperformance_LNonly(NL_Output.crossvaltest_finalrate, recorded_raster, fittedGLM_preNL.t_bin);
                clear fittedGLM recorded_raster
                
                % hack to save components in a somewhat consistent manner                
                fittedGLM.cell_savename        = cell_savename;
                fittedGLM.cellinfo             = fittedGLM_preNL.cellinfo;
                fittedGLM.nonlinearity         = postfilterNL.type;
                fittedGLM.NL_Output            = NL_Output;
                fittedGLM.xvalperformance      = NL_xvalperformance;
                
                fittedGLM.stimtransform.normalized_filteroutput_fit= NL_Input.input_fit;
                fittedGLM.stimtransform.cif_withNL_fit             = NL_Output.fit_rate;
                fittedGLM.stimtransform.cif_rawGLM_fit             = exp( std(baseGLM.lcif_fit.stim) * NL_Input.input_fit  + baseGLM.lcif_fit.mu(1) );
                fittedGLM.stimtransform.cif_withNL_test            = NL_Output.crossvaltest_finalrate;
                fittedGLM.stimtransform.cif_rawGLM_test            = exp( std(baseGLM.lcif_fit.stim) * NL_Input.input_test + baseGLM.lcif_fit.mu(1) );
                fittedGLM.stimtransform.inputNL                    = NL_Input;
                
                fittedGLM.t_bin                = fittedGLM_preNL.t_bin;
                fittedGLM.bins_per_frame       = fittedGLM_preNL.bins_per_frame;
                fittedGLM.rawfit.objective_val = NL_Output.new_objval;
                fittedGLM.fit_time    = datestr(clock);
                fittedGLM.writingcode =  mfilename('fullpath');
                fittedGLM.linearfilters.Stimulus         = fittedGLM_preNL.linearfilters.Stimulus ;
                fittedGLM.linearfilters.Stimulus_rescale = NL_input.scale_rawtoNLinput;
                fittedGLM.linearfilters.Stimulus_rescalenote = ...
                    'multiply stimfilter by rescaler,convolve with stim, retrieve pre-Nonlinearity rate';
                fittedGLM.d_save = savedir;
                
                savename = sprintf('%s/%s',savedir, cell_savename);
                eval(sprintf('save %s.mat fittedGLM fittedGLM_preNL', savename));
                
                if runoptions.print
                    subR_plotfittedNL(fittedGLM, fittedGLM_preNL, savedir)
                end
                

            end            
        end
    end
end

end

function subR_plotfittedNL(fittedGLM, fittedGLM_preNL, savedir)

% Cleaned up AKHeitman 2015-06-24
homedir = pwd;
clf;  
printname = sprintf('DiagNLPlot_%s',fittedGLM.cell_savename);
info    = fittedGLM.cellinfo;
GLMType = fittedGLM_preNL.GLMType;

% text box
subplot(3,1,1)
axis off
set(gca, 'fontsize', 12)
obj_NEW = fittedGLM.rawfit.objective_val;
obj_OLD = fittedGLM_preNL.rawfit.objective_val;
optNL_describe  = fittedGLM.NL_Output.note_metric;
optNL_string    = fittedGLM.NL_Output.param_string;

c = 0; offset = 0; delta_c = 1.1;
text(-offset, 1-0.1*c,sprintf('%s: %s %d: %s-Fit (red): POSTNL refit with %s',...
    info.exp_nm, info.celltype,info.cid, GLMType.fit_type, fittedGLM.nonlinearity), 'interpreter','none')
c = c + delta_c;
text(-offset, 1-0.1*c,sprintf('Red is original GLM, Blue includes Postfilter Nonlinearity'))
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
x1 = sort(fittedGLM.stimtransform.normalized_filteroutput_fit); 
y1 = sort(fittedGLM.stimtransform.cif_rawGLM_fit);
y2 = sort(fittedGLM.stimtransform.cif_withNL_fit);

LW = 2;
subplot(3,3,4); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([-4,4]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity - Central Portion')

subplot(3,3,5); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([ max(min(x1),-10),0]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Inhibitory Portion')

subplot(3,3,6); set(gca, 'fontsize', 10); hold on;
plot(x1,y1,'r','linewidth',LW); 
plot(x1,y2,'b','linewidth',LW);
xlim([0,min(10,max(x1))]);
ylabel('Output (Hz)')
xlabel('StimFilter / std(StimFilter)')
title('Nonlinearity: Excitatory Portion')


% plot rasters
subplot(3,1,3); set (gca, 'fontsize',10)
secs     = 6;
dt = fittedGLM.t_bin;
bins     = 120 * 8 * fittedGLM.bins_per_frame;
rec_rast = fittedGLM.xvalperformance.rasters.recorded(:,1:bins);
glm_rast = fittedGLM_preNL.xvalperformance.rasters.glm_sim(:,1:bins);  
NL_rast  = fittedGLM.xvalperformance.rasters.glm_sim(:,1:bins); 
trials   = size(rec_rast,1);
time     = dt*[1:bins];
xlim([0, ceil(time(end))]);
ylim([0 , 3*trials]); hold on
for i_trial = 1:trials
    rec1 = time(find(rec_rast(i_trial,:)));
    glm1 = time(find(glm_rast(i_trial,:)));
    NL1  = time(find(NL_rast(i_trial,:)));
    % Plot the raster
    plot(rec1, i_trial, 'k.')


    yshift = i_trial;
    if length(glm1) < 4*length(rec1) 
        if length(glm1) > 0
            plot(glm1, yshift + trials, 'r.')
        end
    end
    if length(NL1) < 4*length(rec1) 
        if length(NL1) > 0
            plot(NL1, yshift + 2*trials, 'b.')
        end
    end
end
xlabel('seconds'); ylabel('trials')

cd(savedir)
orient landscape
eval(sprintf('print -dpdf %s.pdf',printname))
cd(homedir)
end
function NL_xvalperformance     = subR_xvalperformance_LNonly( stimdrivenrate, logicalspike, t_bin)
% AKHEITMAN 2015-06-24  it works!
params.bindur     = t_bin;
params.bins       = length(stimdrivenrate);
params.trials     = size(logicalspike,1);
params.testdur_seconds = params.bindur * params.bins ;   

% Set log-conditional as stim driven only
lcif_teststim = log(stimdrivenrate);
lcif = repmat(lcif_teststim, params.trials,1);  


% FOR NOW WE IGNORE PS FILTER  LN ONLY!!
%{
if PostSpikeFilter
    PS = otherfilters_NEW.PS; 
    lcif_ps = fastconv(logicalspike , [0; PS]', size(logicalspike,1), size(logicalspike,2) );    
    lcif = lcif + lcif_ps;
end
%}

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

% FOR NOW WE IGNORE PS FILTER  LN ONLY!!
%{
if PostSpikeFilter
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
%}
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
NL_xvalperformance.rasters.note           = 'glmsim includes altered non-linearity';
NL_xvalperformance.rasters.recorded       = logicalspike;
NL_xvalperformance.rasters.glm_sim        = logical_sim;
NL_xvalperformance.rasters.bintime        = params.bindur;
end              
function [lcif_stim]            = subR_lcifstim_fittedGLM(pstar, GLMType,GLMPars,fitmovie,inputstats,glm_cellinfo)

if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;
clear bin_size basis_params
% PREPARE PARAMETERS
[paramind] =  prep_paramindGP(GLMType, GLMPars); 

% ORGANIZE STIMULUS COVARIATES
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
end
function [lcif, obj_val]        = subR_lcifdecomp_fittedGLM(pstar, GLMType,GLMPars,fitspikes,fitmovie,inputstats,glm_cellinfo)
if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
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
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins  = home_spbins(find(home_spbins < bins) );
if GLMType.PostSpikeFilter
    basis         = ps_basis';
    PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
end
if GLMType.CouplingFilters;
    basis = cp_basis';
    display('figure out coupling here!  CP_bin');
end
if GLMType.TonicDrive
    MU_bin = ones(1,bins);
end


% PREPARE PARAMETERS
[paramind] =  prep_paramindGP(GLMType, GLMPars); 



% ORGANIZE STIMULUS COVARIATES
center_coord       = glm_cellinfo.slave_centercoord;
WN_STA             = double(glm_cellinfo.WN_STA);
[X_frame,X_bin]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord, WN_STA);
clear WN_STA center_coord



%% Run through optimization .. get out pstart, fstar, eflag, output
% CONVEXT OPTIMIZATION
if GLMType.CONVEX
    glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
    % Maybe move this inside of the stimulus preparation % 
    bpf         = GLMPars.bins_per_frame;
    shifts      = 0:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    if isfield(GLMPars.stimfilter,'frames_negative')
        shifts = -(GLMPars.stimfilter.frames_negative)*bpf:bpf:(GLMPars.stimfilter.frames-1)*bpf;
    end
    X_bin_shift = prep_timeshift(X_bin,shifts);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isfield(paramind, 'MU')
        glm_covariate_vec( paramind.MU , : ) = MU_bin;
    end
    if isfield(paramind, 'X')
        glm_covariate_vec( paramind.X , : ) = X_bin_shift;
    end
    if isfield(paramind, 'PS')
        glm_covariate_vec( paramind.PS , : ) = PS_bin;
    end
    pstar = pstar';
    lcif.mu   = pstar(paramind.MU)*glm_covariate_vec(paramind.MU,:);
    lcif.stim = pstar(paramind.X) *glm_covariate_vec(paramind.X,:); 
    total_lcif = lcif.mu + lcif.stim;
    
    if GLMType.PostSpikeFilter
        lcif.ps   = pstar(paramind.PS)*glm_covariate_vec(paramind.PS,:);    
        lcif.ps_unoptimized.glm_covariate_vec = glm_covariate_vec(paramind.PS,:);
        lcif.ps_unoptimized.basis             = ps_basis;
        total_lcif = total_lcif + lcif.ps;
    end
    
    cif        = exp(total_lcif);
    obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));

end

end
function [f grad Hess log_cif]  = subr_optimizationfunction(linear_params,covariates,spikebins,bin_duration)
p = linear_params;
COV = covariates;
dt = bin_duration;
spt = spikebins;
% Find Conditional Intensity and its log
lcif = p' * COV;
cif  = exp(lcif);
% Evaluate the objective function (monotonic in log-likelihood)
f_eval = sum( lcif(spt) ) - dt * sum(cif);
% Evaluate the gradient
g_eval = sum(COV(:,spt),2)  - dt * ( COV * (cif') );
% Evaluate the hessian
hessbase = zeros(size(COV));
for i_vec = 1:size(COV,1)
    hessbase(i_vec,:) = sqrt(cif) .* COV(i_vec,:) ;
end
H_eval = -dt * (hessbase * hessbase');
% Switch signs because using a minimizer  fmin
f       = -f_eval;
grad    = -g_eval;
Hess    = -H_eval;
log_cif = lcif;
end
function spikesconcat           = subR_concat_fitspikes_fromorganizedspikes(blockedspikes, FitPars)
% AKHeitman 2014-04-14
% Concatenate Spikes from different blocks to a single spike train
% blocekdspikes: needs
%   .t_sp_withinblock
%
% FitPars needs
%   .fittest_skipseconds
%   .tstim
%   .fitframes
%   .FitBlocks


t_start   = FitPars.fittest_skipseconds;
tstim     = FitPars.computedtstim;
fitframes = FitPars.fitframes;
FitBlocks = FitPars.FitBlocks;


T_SP = []; blk_count = 0;
dur = tstim * length(fitframes);
for k = FitBlocks
	blk_count = blk_count + 1;
	t_sp_full = blockedspikes.t_sp_withinblock{k} ; % unit of time: sec, 0 for the onset of the block
	t_sp      = t_sp_full(find(t_sp_full >  t_start));
	t_sp = t_sp - t_start;
	t_spcontext = t_sp + ( blk_count -1 )*dur;
	T_SP = [T_SP ; t_spcontext];
end
spikesconcat = T_SP;
end
function raster_spiketimes      = subR_createraster(blockedspikes, TestPars)
% AKHeitman 2014-04-14
% Make a raster which takes into account GLM processing
% blocekdspikes: needs
%   .t_sp_withinblock
%
% TestPars needs
%   .fittest_skipseconds
%   .TestBlocks

rasterblocks = TestPars.TestBlocks;
t_start      = TestPars.fittest_skipseconds;

raster_spiketimes = cell(length(rasterblocks),1);

for i_blk = 1 : length(rasterblocks)
	blknum = rasterblocks(i_blk);
	sptimes = blockedspikes.t_sp_withinblock{blknum} - t_start;
	sptimes = sptimes(find(sptimes > 0 ) );
    % HACK NEEDED FOR 2013-10-10-0 and other long runs
    if isfield(TestPars, 'test_skipENDseconds')
        sptimes = sptimes(find(sptimes < (TestPars.test_skipENDseconds - TestPars.fittest_skipseconds - .1)));
    end
    
    raster_spiketimes{i_blk} = sptimes;
end 

end
function concat_fitmovie        = subR_concat_fitmovie_fromblockedcell(blockedmoviecell , FitPars)
% AKHeitman 2014-04-14
% Concatenate the fit movie (different blocks)
% FitPars needs
%   .width
%   .height
%   .FitBlocks
%   .novelblocks
%   .fitframes

height       = FitPars.height;
width        = FitPars.width;
fitblocks    = FitPars.FitBlocks;
fitframes    = FitPars.fitframes;
novelblocks  = FitPars.NovelBlocks;

fitframesperblock = length(fitframes) ;
totalframes       = length(fitblocks) * ( fitframesperblock) ;
concat_fullfitMovie = uint8(zeros(width, height, totalframes)) ;
for i_blk = fitblocks
        blkind = find(fitblocks == i_blk);
        framenums = ( (blkind -1)*fitframesperblock + 1 ) :  (blkind *fitframesperblock);  
        n_blkind = find(novelblocks == i_blk);
        concat_fullfitMovie(:,:,framenums) = blockedmoviecell{n_blkind}.matrix (:,:, fitframes);    
end

concat_fitmovie = concat_fullfitMovie;

end
function [center,sd]            = subR_visionSTA_to_xymviCoord(stafit_centercoord, stafit_sd, masterdim, slvdim)
% AKHeitman  2013-12-08
% Grab x, y coordinates of STA center of the master
% Convert to coordinates of the enslaved dataset 
x_coord   = round( stafit_centercoord(1)* (slvdim.width  /masterdim.width)  );
y_coord   = slvdim.height - round( stafit_centercoord(2)* (slvdim.height /masterdim.height) );

center.x_coord = x_coord;
center.y_coord = y_coord;

sd.xdir = round( stafit_sd(1)* (slvdim.width   / masterdim.width)  );
sd.ydir = round( stafit_sd(2)* (slvdim.height  / masterdim.height)  );

end
function obj_val = subR_rescale_stim(scalars, lcif_stim0, home_spbins,t_bin)

offset    = scalars(1);
lcif_stim = scalars(2) * lcif_stim0;

total_lcif = offset + lcif_stim;
cif        = exp(total_lcif);
obj_val    = -(sum( total_lcif(home_spbins) ) - t_bin * sum(cif));
end



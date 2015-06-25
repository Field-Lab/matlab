%{

%}
function read_NL_fits(exps,stimtypes,celltypes,cell_subset,baseGLM_settings,postfilterNL,runoptions)
% started: AKHeitman 2015-06-25


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
        baseGLM.Type.fitmoviefile  = origmatfile;
        if postfilterNL.debug
            display('shorten stimulus for post filter debugging mode')
            StimulusPars.slv.FitBlocks = StimulusPars.slv.FitBlocks(1:2);
        end
                
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



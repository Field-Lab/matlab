load('/Volumes/Analysis/nora/NSEM/GLM_Output/old_fits/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/NSEM_mapPRJ/2012-08-09-3/OFFPar_1471.mat');
fittedGLM.GLMType.Subunits=false;
% Load stim params
[StimulusPars, DirPars, ~, ~] = Directories_Params_v23(fittedGLM.cellinfo.exp_nm, fittedGLM.GLMType.fit_type, fittedGLM.GLMType.map_type);

% Load spikes
eval(sprintf('load %s/organizedspikes_%s.mat organizedspikes', '/Volumes/Analysis/nora/NSEM/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ', fittedGLM.cellinfo.cell_savename));

% Load the test movie
[testmovie] = loadmoviematfile(fittedGLM.cellinfo.exp_nm , fittedGLM.GLMType.fit_type, fittedGLM.GLMType.cone_model,'testmovie');

% Calculate
xvalperformance = eval_xvalperformance_PP(fittedGLM, StimulusPars.slv, organizedspikes,0,testmovie);

% Vis
time=fittedGLM.t_bin*(1:length(xvalperformance.rate));
plot(time,xvalperformance.rate)
xlabel('Time (seconds)')
ylabel('Poisson parameter')
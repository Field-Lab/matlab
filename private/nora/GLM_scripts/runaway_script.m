
%%
datapath='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/';
exp_names=['2012-08-09-3/CP_PCA/';'2013-08-19-6/CP_PCA/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

% Experiment 1
cells = 'ON';
Iexp = 1;

% Experiment 3
% Iexp = 2;
% cells = 'OFF';

% Just looking at NSEM
fittype = 2;

% Get file list
matfiles=dir([datapath fittypepath{fittype} exp_names(Iexp,:) cells '*.mat']);
n_cells=length(matfiles);

%Initialize matrices
runaway_trials = zeros(n_cells);

% Collect info from files
for file=1:n_cells
    disp(file)
    load([datapath fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
    runaway_trials(file) = sum(sum(fittedGLM.xvalperformance.rasters.glm_sim') > 4*sum(fittedGLM.xvalperformance.rasters.recorded'));
end
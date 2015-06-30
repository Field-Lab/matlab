
%%
datapath{1}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_CP_p8IDp8/standardparams/';
datapath{2}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/fixedSP_rk1_linear_MU_PS_noCP_p8IDp8/standardparams/';
datapath{3}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_CP_p8IDp8/standardparams/';
datapath{4}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk1_MU_PS_noCP_p8IDp8/standardparams/';

exp_names=['2012-08-09-3/CP_PCA/';'2013-08-19-6/CP_PCA/'];
fit_type{1}='WN';
fit_type{2}='NSEM';
fittypepath{2}=[fit_type{2} '_mapPRJ/'];
fittypepath{1}=[fit_type{1} '_mapPRJ/'];

%%
% Experiment 1
cells = 'ON';
Iexp = 1;

% Experiment 3
% Iexp = 2;
% cells = 'OFF';

% Just looking at NSEM
fittype = 2;

runaway_trials = zeros(40, 4);
BPS = zeros(40,4);

for fit = 1:4
% Get file list
matfiles=dir([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) cells '*.mat']);
n_cells=length(matfiles)

% Collect info from files
for file=1:n_cells
    load([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
    runaway_trials(file, fit) = sum(sum(fittedGLM.xvalperformance.rasters.glm_sim') > 4*sum(fittedGLM.xvalperformance.rasters.recorded'));
    BPS(file, fit) = fittedGLM.xvalperformance.glm_normedbits;
end

end

%%
hold on
[counts{1}, centers] = hist(runaway_trials(:,1));
counts{2} = hist(runaway_trials(:,2), centers);
counts{3} = hist(runaway_trials(:,3), centers);
counts{4} = hist(runaway_trials(:,4), centers);
hold off

hold on
for i = 1:4
    plot(centers, counts{i})
end

legend('fixedSP CP', 'fixedSP noCP', 'fitSP CP', 'fitSP no CP')
ylabel('Number of Cells')
xlabel('Number of Runaway Trials')

%% Find OFF runaway trials in rk1, CP
% Experiment 3
Iexp = 2;
cells = 'OFF';
% Just looking at NSEM
fittype = 2;
for fit = 3
% Get file list
matfiles=dir([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) cells '*.mat']);
n_cells=length(matfiles);

% Collect info from files
for file=1:n_cells
    load([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
    runaway_trials = sum(sum(fittedGLM.xvalperformance.rasters.glm_sim') > 4*sum(fittedGLM.xvalperformance.rasters.recorded'));
    if runaway_trials>0
        disp(matfiles(file).name)
    end
end

end

%%
datapath{5}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_CP_p8IDp8/standardparams/';
datapath{6}='/Volumes/Lab/Users/Nora/NSEM_Home/GLMOutput_Raw/rk2_MU_PS_noCP_p8IDp8/standardparams/';

% Experiment 3
Iexp = 2;
cells = 'OFF';

% Just looking at NSEM
fittype = 2;

runaway_trials = zeros(120, 6);

for fit = 5:6
% Get file list
matfiles=dir([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) cells '*.mat']);
n_cells=length(matfiles)

% Collect info from files
for file=1:n_cells
    load([datapath{fit} fittypepath{fittype} exp_names(Iexp,:) matfiles(file).name]);
    runaway_trials(file, fit) = sum(sum(fittedGLM.xvalperformance.rasters.glm_sim') > 4*sum(fittedGLM.xvalperformance.rasters.recorded'));
    cid(file) = fittedGLM.cellinfo.cid;
end

end

hold on
counts{1}, centers] = hist(runaway_trials(:,1));
counts{2} = hist(runaway_trials(:,2), centers);
hold off

hold on
for i = 1:4
    plot(centers, counts{i})
end

legend('fixedSP CP', 'fixedSP noCP', 'fitSP CP', 'fitSP no CP')
ylabel('Number of Cells')
xlabel('Number of Runaway Trials')
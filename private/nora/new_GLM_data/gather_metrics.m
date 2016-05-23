clear
runs = {...
'2012-08-09-3',...
'2015-05-27-11',...
'2016-02-17-1',...
'2014-09-10-1',...
'2015-10-06-0',...
'2015-05-27-3',...
'2012-09-27-3',...
'2013-08-19-6'};

Isolated = [0 0 0 0 0 1 1 1];
n_runs = length(Isolated);
Scores = cell(n_runs,1);

for i_run = 1:n_runs
    folder = runs{i_run}(runs{i_run}~='-');
    if Isolated(i_run)
        type = 'Isolated';
    else
        type = 'RPE';
    end
    
    %%{
    files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/' files(i_cell).name])
        try
            FR = sum(fittedGLM.xval.rasters.recorded(:))/numel(fittedGLM.xval.rasters.recorded);
            Scores{i_run} = [Scores{i_run} FR];
        catch
             FR = sum(fittedGLM.xvalperformance.rasters.recorded(:))/numel(fittedGLM.xvalperformance.rasters.recorded);
             Scores{i_run} = [Scores{i_run} FR];
        end
    end
    %}
%%{
        files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnP*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
        try
            FR = sum(fittedGLM.xval.rasters.recorded(:))/numel(fittedGLM.xval.rasters.recorded);
            Scores{i_run} = [Scores{i_run} FR];
        catch
            FR = sum(fittedGLM.xvalperformance.rasters.recorded(:))/numel(fittedGLM.xvalperformance.rasters.recorded);
             Scores{i_run} = [Scores{i_run} FR];
        end
    end
    %}
    %%{
        files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/ONP*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
        try
            FR = sum(fittedGLM.xval.rasters.recorded(:))/numel(fittedGLM.xval.rasters.recorded);
            Scores{i_run} = [Scores{i_run} FR];
        catch
            FR = sum(fittedGLM.xvalperformance.rasters.recorded(:))/numel(fittedGLM.xvalperformance.rasters.recorded);
             Scores{i_run} = [Scores{i_run} FR];
        end
    end
    %}
end
%%
for i = 1:n_runs
    run_scores = Scores{i};
    run_scores = run_scores(run_scores > 0);
    exp_std(i) = std(run_scores);
    exp_means(i) = mean(run_scores);
end
b = bar(1:2, exp_means(1:2), 'w', 'EdgeColor', 'r', 'LineWidth', 5);
hold on;
b = bar(3:5, exp_means(3:5), 'w');
b = bar(find(Isolated), exp_means(find(Isolated)), 'y');
errorbar(exp_means, exp_std, '.k', 'LineWidth', 5)
ylabel('Correlation Coefficient')
xlabel('Experiment')


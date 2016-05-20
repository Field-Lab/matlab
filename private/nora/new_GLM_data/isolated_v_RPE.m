clear
runs = {...
'2012-08-09-3',...
'2015-05-27-11',...
'2016-02-17-1',...
'2014-09-10-1',...
'2016-01-05-0',...
'2015-10-06-0',...
'2015-05-27-3',...
'2012-09-27-3',...
'2013-08-19-6'};

Isolated = [0 0 0 0 0 0 1 1 1];
n_runs = length(Isolated);
Scores = cell(n_runs,1);

for i_run = 1:n_runs
    folder = runs{i_run}(runs{i_run}~='-');
    if Isolated(i_run)
        type = 'Isolated';
    else
        type = 'RPE';
    end
    
    %{
    files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/' files(i_cell).name])
        try
            Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
        catch
             Scores{i_run} = [Scores{i_run} fittedGLM.xval.corr];
        end
    end
    %}
%%{
        files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OFFP*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
        try
            Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
        catch
             Scores{i_run} = [Scores{i_run} fittedGLM.xval.corr];
        end
    end
    %}
    %%{
        files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OffP*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
        try
            Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
        catch
             Scores{i_run} = [Scores{i_run} fittedGLM.xval.corr];
        end
    end
    %}
end
%%
for i = 1:n_runs; exp_std(i) = std(cell2mat(Scores(i)')); end
for i = 1:n_runs; exp_means(i) = mean(cell2mat(Scores(i)')); end
b = bar(1:2, exp_means(1:2), 'w', 'EdgeColor', 'r', 'LineWidth', 5);
hold on;
b = bar(3:5, exp_means(3:6), 'w');
b = bar(6:8, exp_means(Isolated), 'y');
errorbar(exp_means, exp_std, '.k', 'LineWidth', 5)
ylabel('Correlation Coefficient')
xlabel('Experiment')


runs = {...
'2012-08-09-3',...
'2015-10-06-0',...
'2016-02-17-1',...
'2015-05-27-11',...
'2015-05-27-3',...
'2012-09-27-3',...
'2013-08-19-6'};

Isolated = [0 0 0 0 1 1 1];
n_runs = length(Isolated);
Scores = cell(n_runs,1);

for i_run = 1:n_runs
    folder = runs{i_run}(runs{i_run}~='-');
    if Isolated(i_run)
        type = 'Isolated';
    else
        type = 'RPE';
    end
    files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/' files(i_cell).name])
        Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
    end
%%{
        files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/*.mat']);
    for i_cell = 1:length(files)
        load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
        Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
    end
    %}
end


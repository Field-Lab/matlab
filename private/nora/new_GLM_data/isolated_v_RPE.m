clear
runs = {...
    '2012-08-09-3',...
    '2015-05-27-11',...
    '2016-02-17-1',...
    '2014-09-10-1',...
    '2015-01-05-0',...
    '2015-10-06-0',...
    '2015-05-27-3',...
    '2012-09-27-3',...
    '2013-08-19-6'};

Isolated = [0 0 0 0 0 0 1 1 1];
n_runs = length(Isolated);
Scores = cell(n_runs,1);
Firing_Rates = cell(n_runs,1);
Variances = cell(n_runs,1);
FR_max = cell(n_runs, 1);
Sparsity = cell(n_runs, 1);
HS_exp = cell(n_runs, 1);
for i_run = 1:n_runs
    folder = runs{i_run}(runs{i_run}~='-');
    if Isolated(i_run)
        type = 'Isolated';
    else
        type = 'RPE';
    end
    
    %%{
    for i = 1:3
        switch i
            case 1
                files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnP*.mat']);
            case 2
                files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/ONP*.mat']);
            case 3
                files = dir(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/*.mat']);
        end
        for i_cell = 1:length(files)
            if i == 3
                load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/OnPar/' files(i_cell).name])
            else
                load(['/Volumes/Lab/Users/Nora/GLMFits/' type '/' folder '/WN/' files(i_cell).name])
            end
            try
                Scores{i_run} = [Scores{i_run} fittedGLM.xvalperformance.corr];
                FR = sum(fittedGLM.xvalperformance.rasters.recorded(:))/numel(fittedGLM.xvalperformance.rasters.recorded);
                FRmax = max(sum(fittedGLM.xvalperformance.rasters.recorded));%/numel(fittedGLM.xvalperformance.rasters.recorded);
                FR_max{i_run} = [FR_max{i_run}, FRmax];
                var = std(sum(fittedGLM.xvalperformance.rasters.recorded));%/numel(fittedGLM.xvalperformance.rasters.recorded);
                Firing_Rates{i_run} = [Firing_Rates{i_run} FR];
                Variances{i_run} = [Variances{i_run} var];
                S = sum((sum(fittedGLM.xvalperformance.rasters.recorded))==0);
                Sparsity{i_run} = [Sparsity{i_run} S];
                %HS = expfit(sum(fittedGLM.xvalperformance.rasters.recorded));
                HS = mean(abs(diff(sum(fittedGLM.xvalperformance.rasters.recorded))));
                HS_exp{i_run} = [HS_exp{i_run} HS];
            catch
                Scores{i_run} = [Scores{i_run} fittedGLM.xval.corr];
                FR = sum(fittedGLM.xval.rasters.recorded(:))/numel(fittedGLM.xval.rasters.recorded);
                var = std(sum(fittedGLM.xval.rasters.recorded));%/numel(fittedGLM.xval.rasters.recorded);
                Firing_Rates{i_run} = [Firing_Rates{i_run} FR];
                Variances{i_run} = [Variances{i_run} var];
                FRmax = max(sum(fittedGLM.xval.rasters.recorded));%/numel(fittedGLM.xvalperformance.rasters.recorded);
                FR_max{i_run} = [FR_max{i_run}, FRmax];
                S = sum((sum(fittedGLM.xval.rasters.recorded))==0);
                Sparsity{i_run} = [Sparsity{i_run} S];
                %HS = expfit(sum(fittedGLM.xval.rasters.recorded));
                HS = mean(abs(diff(sum(fittedGLM.xval.rasters.recorded))));
                HS_exp{i_run} = [HS_exp{i_run} HS];
            end
        end
    end
    %}
    
end
%%
for i = 1:n_runs
    run_scores = Scores{i};
    run_scores = run_scores(run_scores > 0.2);
    exp_std(i) = std(run_scores);
    exp_means(i) = mean(run_scores);
    % firing rate
    FR_std(i) = std(Firing_Rates{i}(run_scores>0.2));
    FR_means(i) = mean(Firing_Rates{i}(run_scores>0.2));
    % structure in PSTH
    var_std(i) = std(Variances{i}(run_scores>0.2));
    var_means(i) = mean(Variances{i}(run_scores>0.2));
    
    
    max_std(i) = std(FR_max{i}(run_scores>0.2));
    max_mean(i) = mean(FR_max{i}(run_scores>0.2));
    
    % sparsity
    S_std(i) = std(Sparsity{i}(run_scores>0.2));
    S_mean(i) = mean(Sparsity{i}(run_scores>0.2));
    
    
    HS_std(i) = std(HS_exp{i}(run_scores>0.2));
    HS_mean(i) = mean(HS_exp{i}(run_scores>0.2));
    
    % trial to trial differences
end

%%
figure;
b = bar(1:2, exp_means(1:2), 'w', 'EdgeColor', 'r', 'LineWidth', 5);
hold on;
b = bar(3:6, exp_means(3:6), 'w');
b = bar(find(Isolated), exp_means(find(Isolated)), 'y');
errorbar(exp_means, exp_std, '.k', 'LineWidth', 5)
ylabel('Correlation Coefficient')
xlabel('Experiment')

%%
figure;
plot(FR_means*10^3, exp_means, '.', 'MarkerSize', 10)
text(FR_means*10^3, exp_means, {'1', '2', '3', '4', '5', '6', '7', '8', '9'}, 'FontSize', 20)
ylim([0.3 1])
ylabel('Correlation Coefficient')
xlabel('Average Firing Rate')

%%
figure;
plot(HS_mean, exp_means, '.', 'MarkerSize', 10)
%errorbar(exp_means, S_mean, S_std, '.')
text(HS_mean, exp_means, {'1', '2', '3', '4', '5', '6', '7', '8', '9'}, 'FontSize', 20)
ylim([0.3 1])
ylabel('Mean Correlation Coefficient')
xlabel('Std of Firing Rate Across Trials')

clear 

params_201602171
n_reg = length(reg);
n_groups = length(cells);

Analysis_Path = ['/Volumes/Analysis/' piece '/map-from-' class '/'];
fig_save = ['/Users/Nora/Desktop/Fig_Output/' piece];
mkdir(fig_save);

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])

% classification run
class_datarun = load_data(['/Volumes/Analysis/' piece '/streamed/' class '/' class]);
class_datarun = load_params(class_datarun);
class_datarun = load_neurons(class_datarun);

%%
for i_group =1:n_groups
    
    n_subgroups = length(cell_idx{1});
    
    %%{
    % load reg repeats
    for i_reg = 1:n_reg
        datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
        datarun{i_reg} = load_neurons(datarun{i_reg});
        reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells{i_group}, 'datarun_class', class_datarun, 'visual_check', 0);
        clear datarun 
    end
    %}
    
    %{
    % plot reg repeats for stability check
    for i_cell = 1:size(reg_data{1}.testspikes,2)
        figure; hold on
        for i_reg = 1:n_reg
            IDP_plot_PSTH(reg_data{i_reg}, i_cell);
        end
        title(['cell ' num2str(cells{i_group}(i_cell))]);
    end
    %}
   
    datarun = load_data([ Analysis_Path masking{i_group} '/' masking{i_group}]);
    datarun = load_neurons(datarun);
    mask_data = interleaved_data_prep(datarun, 1100, n_masks{i_group}*2*30, 'cell_spec', cells{i_group}, 'datarun_class', class_datarun, 'visual_check', 1);
    pause()
    clear datarun  
%{
    % replace with cellfindered data if an option
datarun_cf = load_data([ Analysis_Path masking{i_group} '-cf/' masking{i_group} '-cf']);
datarun_cf = load_neurons(datarun_cf);
mask_data_cf = interleaved_data_prep(datarun_cf, 1100, 480, 'cell_spec', 'all', 'visual_check', 0);

for i_cell = 1:length(cells{i_group})
    try
        mask_data.testspikes(:,i_cell) = mask_data_cf.testspikes(:,cells{i_group}(i_cell) == datarun_cf.cell_ids);
    catch
        disp(['cell ' num2str(cells{i_group}(i_cell)) ' not cell findered'])
    end
end
clear datarun_cf mask_data_cf
%}
    
    % separate each condition
    idx = 1:30;
    for count = 1:(n_masks{i_group}*2)
        condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
    end
    
    % plot the PSTH for each condition for each cell
    for subgroup = 1:n_subgroups
        for i_cell = cell_idx{i_group}{subgroup}
            for i = 1:length(mask_conditions{i_group}{subgroup})
                color = [1 1 1]*(1-sigmas{i_group}(mask_conditions{i_group}{subgroup}(i))/12);
                figure(1); hold on
                IDP_plot_PSTH(condition{mask_conditions{i_group}{subgroup}(i)},i_cell, color);
                figure(2); hold on
                IDP_plot_PSTH(condition{comp_conditions{i_group}{subgroup}(i)},i_cell, color);
                %exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_cond' num2str(i) '_sigma' num2str(sigmas(i))], 'Bounds', 'loose', 'Color', 'rgb')
            end
            figure(1); IDP_plot_PSTH(reg_data{1}, i_cell, [0 0 0]); title('Mask')
            figure(2); title('Complement');
            pause()
            figure(1); clf; figure(2); clf;
        end
    end
    close all
end
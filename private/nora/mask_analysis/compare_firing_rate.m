clear
total_mask = [];
total_reg = [];
celltype_idx = [];

set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 8 4])

params_list

%sigma_idx = [0 1 0 2 3 4 0 5 0 6];

%%
for exp = exps
    disp(exp{1})
    eval(exp{1})
    
    n_reg = length(reg);
    mkdir(fig_save);

    n_subgroups = length(cell_idx);
    
    %%
    %%{
    % load reg repeats
    for i_reg = 1:n_reg
        datarun{i_reg} = load_data([ Analysis_Path reg{i_reg} '/' reg{i_reg}]);
        datarun{i_reg} = load_neurons(datarun{i_reg});
        reg_data{i_reg} = interleaved_data_prep(datarun{i_reg}, 1100, 30, 'cell_spec', cells_reg{i_reg}, 'visual_check', 0);
        clear datarun
    end
    %}
    
    %{
    % plot reg repeats for stability check
    for i_cell = 1:size(reg_data{1}.testspikes,2)
        figure; hold on
        for i_reg = 1:n_reg
            IDP_plot_PSTH(reg_data{i_reg}, i_cell, 'color', i_reg);
        end
        title(['cell ' num2str(i_cell)]);
    end
    %}
    % close all
    
    datarun = load_data([ Analysis_Path masking '/' masking]);
    datarun = load_neurons(datarun);
    mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
    clear datarun
    
    % separate each condition
    idx = 1:30;
    for count = 1:(n_masks*2)
        condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
    end
%
    % plot the PSTH for each condition for each cell
    for subgroup = 1:n_subgroups
        for i_cell = cell_idx{subgroup}
            disp(i_cell)
            reg_psth = IDP_plot_PSTH(reg_data{1},i_cell, 'color', 0);
            reg_spike_count(i_cell) = var(reg_psth)/length(reg_psth);
            for i = 1:length(mask_conditions{subgroup})
                s = sigmas(mask_conditions{subgroup}(i));
                mask_psth =IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0);
                mask_spike_count(i_cell, i) = var(mask_psth)/length(reg_psth); 
            end
        end
        if strcmp(cell_type{subgroup}, 'on')
             celltype_idx = [celltype_idx 1 1 1];
        else
            celltype_idx = [celltype_idx 0 0 0 ];
        end
    end
    
    total_mask = [total_mask; mask_spike_count];
    total_reg = [total_reg; reg_spike_count'];

%%
%{ 
set(0, 'defaultFigurePosition', [2 2 4 4])
figure(10);
hold on; plot([2 4 5 6], [mask_spike_count']./repmat(reg_spike_count, [4 1]))
%}
%{
hold on; plot(mask_spike_count, '.r');
hold on; plot(mask_spike_count(~logical(celltype_idx)), reg_spike_count(~logical(celltype_idx)), '.b');
hold on; plot([0 300], [0 300], 'k')
axis([0 30 0 30]);
axis square
%}

end

%%
total_mask = total_mask([1:2, 4:end], :);
total_reg = total_reg([1:2, 4:end]);
celltype_idx = celltype_idx([1:2, 4:end]);

plot([2 4 5 6], total_mask'./repmat(total_reg, [1 4])','Color', [1 1 1]*0.75)
hold on; plot([2 4 5 6], mean(total_mask./repmat(total_reg, [1 4])),'k', 'LineWidth', 2)

%a = fitlm(mask_spike_count,reg_spike_count);
%plot(a)
%exportfig(gcf, [fig_save '/' 'cell' num2str(i_cell) '_comp'], 'Bounds', 'loose', 'Color', 'rgb')

ylim([0.5 2.5])
plot([2 6], [1 1], '--')

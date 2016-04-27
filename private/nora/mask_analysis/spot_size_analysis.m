clear
params_list
set(0, 'defaultFigurePaperPositionMode', 'auto')
set(0, 'defaultFigurePaperOrientation', 'landscape')
set(0, 'defaultFigureUnits', 'inches')
set(0, 'defaultFigurePosition', [2 2 15 4])
mean_cell_corr = cell(10,1);
mean_cell_surr = cell(10,1);

for exp = exps
    disp(exp{1})
    eval(exp{1})
    n_reg = length(reg);
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
    
    %%
    datarun = load_data([ Analysis_Path masking '/' masking]);
    datarun = load_neurons(datarun);
    mask_data = interleaved_data_prep(datarun, 1100, n_masks*2*30, 'cell_spec', cells_masking, 'visual_check', 0);
    clear datarun
    
    % separate each condition
    idx = 1:30;
    for count = 1:(n_masks*2)
        condition{count}.testspikes = mask_data.testspikes(idx,:); idx = idx+30;
    end
    
    % plot the PSTH for each condition for each cell
    for subgroup = 1:n_subgroups
        for i_cell = cell_idx{subgroup}
            PSTH_full = IDP_plot_PSTH(reg_data{i_reg}, i_cell, 'color', 0);
            for i = 1:length(mask_conditions{subgroup})
                PSTH_center = IDP_plot_PSTH(condition{mask_conditions{subgroup}(i)},i_cell, 'color', 0);
                NSEM_Corr(i_cell,i) = err(PSTH_full, PSTH_center);
                PSTH_surround = IDP_plot_PSTH(condition{comp_conditions{subgroup}(i)},i_cell, 'color', 0);
                surround_struct(i_cell,i) = var(PSTH_surround);
            end
        end
    end
        
    %%
    sigma = sigmas(mask_conditions{subgroup});
    
    for i = 1:length(cells_masking)
        MSE =  NSEM_Corr(i,:);
        surr = surround_struct(i,:);
        if ~any(diff(MSE)>20) && ~any(diff(surr)>100)
            figure(1); hold on; plot(sigma, MSE, 'Color', [1 1 1]*0.75); hold on
            figure(2); hold on; plot(sigma, surr,  'Color', [1 1 1]*0.75); hold on
            if strcmp(piece, '2016-04-21-1') && i == 1
                example_Cell = {sigma, MSE, surr};
            end
            for s = 1:length(sigma)
                mean_cell_corr{sigma(s)} = [mean_cell_corr{sigma(s)} MSE(s)];
                mean_cell_surr{sigma(s)} = [mean_cell_surr{sigma(s)} surr(s)];
            end
        end
    end
end

%%
default_colors = get(gca,'ColorOrder');
all_sigmas = [2 4 5 6 8 10];
for s = 1:6
   mean_corr(s) = mean(mean_cell_corr{all_sigmas(s)}); 
   mean_surr(s) = mean(mean_cell_surr{all_sigmas(s)}); 
end
figure(1); hold on; plot(all_sigmas, mean_corr, 'k', 'LineWidth', 2)
%hold on; plot(example_Cell{1}, example_Cell{2}, 'Color', default_colors(1,:), 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])

figure(2); hold on; plot(all_sigmas, mean_surr, 'k', 'LineWidth', 2)
%hold on; plot(example_Cell{1}, example_Cell{3}, 'Color', default_colors(1,:), 'LineWidth', 2)
set(gcf, 'Position', [2 2 4 4])
axis square
xlim([2 10])


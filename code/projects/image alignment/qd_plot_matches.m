
% select cells
switch 1
    case 1 % 2007-09-18-4

        axons_to_plot = [9 11 14 15 39 41 44 ];
        % not: 17 65
end

% parameters
fig_num = 2;
num_cells = 2;

axons_to_plot = 15;

% for each, plot best two matches
for axon_id = axons_to_plot

    [all_corr,cell_ids] = find_matches(datarun,axons,axon_id);

    if axon_id == 15; cell_ids = cell_ids([1 4:end]); end
    
    figure(fig_num);clf
    for cc=1:num_cells
        cell_id = cell_ids(cc);
        subplot(num_cells,1,cc)
        qd_plot_aligned_rgc(datarun,axons,cell_id,axon_id,'foa',gca,'ei_scale',1,...
            'ei_cutoff',-1,'title',sprintf('corr = %0.3f, ',all_corr(cc)))
    end
    print(fig_num,sprintf('/snle/home/gauthier2/Desktop/axon%03d',axon_id),'-dpdf')
end


function  IDP_plot_raster(prepped_data, cell_to_plot, firing_rate)

% Visual check
figure;
set(gcf, 'Position', [1 1 6 2])
hold on
if nargin == 2
    for i = 1:size(prepped_data.testspikes, 1)
        plot(prepped_data.testspikes{i, cell_to_plot}, i*ones(length(prepped_data.testspikes{i, cell_to_plot})), 'k.')
    end
else
    
    for i = 1:size(prepped_data.testspikes, 1)
        plot(prepped_data.testspikes{i, cell_to_plot}, i*ones(length(prepped_data.testspikes{i, cell_to_plot})), 'k.')
    end

        dt = (1/119.5);
        time = (dt/10:dt/10:(dt*length(firing_rate)))+(dt*30);
        trials = size(prepped_data.testspikes, 1);
        logicalsim = Poisson_spiking(firing_rate, trials, 10, 119.5);
        start = trials;
        default_colors = get(gca,'ColorOrder');
        for i_trial = 1:trials
            sim1 = time(find(logicalsim(i_trial,:)));
            plot(sim1, i_trial+start, '.', 'Color', default_colors(1,:))
        end
end

    
    xlim([2 8])
    
end
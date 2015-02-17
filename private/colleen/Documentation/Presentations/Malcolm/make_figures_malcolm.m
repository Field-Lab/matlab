results = dataset([],[],[],[],[],[],[],'VarNames',...
    {'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'});
cell_types = {'On parasol','Off parasol','On midget','Off midget'};
% options to make various figures
load_results = true; % this also produces a plot on midget vs on parasol median absolute deviation (MAD)
save_figures = false; % saves figures to file
make_estimate_histogram = true;
make_nearest_neighbor_histogram = true;
make_firing_rate_plot = true; % need to have loaded data as in motion_script for this to work
make_electrode_map = true;
load_results_downsample = true; % also produces a plot of on midget vs on parasol MAD

% load results from 2007-03-27-1 dataset
if load_results
    cd('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/mthielk/data/2007-03-27-1')
    nbins = 20;
    for my_run = 14:19

        % set params
        if(my_run==15) 
            my_config = 3;
        elseif(my_run==17)
            my_config = 2;
        else
            my_config = 1;
        end

        if(my_run==14 || my_run==17)
            my_speed=1;
        elseif(my_run==15 || my_run==18)
            my_speed=2;
        else
            my_speed=3;
        end

        load(sprintf('On parasol_data_run_%d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 1, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-03-27-1, on parasol, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-03-27-1/hist_on_parasol_run%d_config%d',my_run,my_config),'png');
            end
        end

        load(sprintf('Off parasol_data_run_%d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 2, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];     %#ok<*AGROW>
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-03-27-1, off parasol, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-03-27-1/hist_off_parasol_run%d_config%d',my_run,my_config),'png');
            end
        end

        load(sprintf('On midget_data_run_%d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 3, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-03-27-1, on midget, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-03-27-1/hist_on_midget_run%d_config%d',my_run,my_config),'png');
            end
        end

        load(sprintf('Off midget_data_run_%d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 4, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-03-27-1, off midget, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures    
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-03-27-1/hist_off_midget_run%d_config%d',my_run,my_config),'png');
            end
        end

    end

    % load results from 2007-08-24-4 dataset
    cd('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/mthielk/data/2007-08-24-4')
    for my_run = 4:11

        % set params
        if(my_run==5 || my_run==10) 
            my_config = 2;
        elseif(my_run==6)
            my_config = 3;
        else
            my_config = 1;
        end

        if(my_run<6)
            my_speed=3;
        elseif(my_run<8)
            my_speed=2;
        elseif(my_run<10)
            my_speed=1;
        else
            my_speed=1.5;
        end

        load(sprintf('On parasol_data_run_%2.2d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-08-24-4'}, my_run, my_speed, my_config, 1, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-08-24-4, on parasol, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-08-24-4/hist_on_parasol_run%d_config%d',my_run,my_config),'png');
            end
        end

        load(sprintf('Off parasol_data_run_%2.2d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-08-24-4'}, my_run, my_speed, my_config, 2, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-08-24-4, off parasol, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures    
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-08-24-4/hist_off_parasol_run%d_config%d',my_run,my_config),'png');
            end
        end

        load(sprintf('On midget_data_run_%2.2d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-08-24-4'}, my_run, my_speed, my_config, 3, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        if make_estimate_histogram
            figure();
            histfit(estimates,nbins);
            title({sprintf('2007-08-24-4, on midget, run %d, config %d',my_run,my_config),sprintf('mean=%0.2f, mad=%0.2f',mean(estimates),std(estimates))});
            if save_figures
                saveas(gcf,sprintf('/Users/vision/Dropbox/Lab/Development/matlab-standard/private/malcolm/figures/2007-08-24-4/hist_on_midget_run%d_config%d',my_run,my_config),'png');
            end
        end
    end
    
    % make plot of on midget MAD vs on parasol MAD
    sd = double(results(:,7));
    speed = double(results(:,3));
    onp = double(results(:,5))==1;
    onm = double(results(:,5))==3;
    markersize=20;
    figure();
    loglog(sd(onp & speed==1),sd(onm & speed==1),'.m','markersize',markersize);
    xlim([0.03 100]);
    ylim([0.03 100]);
    hold on;
    loglog(sd(onp & speed==1.5),sd(onm & speed==1.5),'.b','markersize',markersize);
    loglog(sd(onp & speed==2),sd(onm & speed==2),'.g','markersize',markersize);
    loglog(sd(onp & speed==3),sd(onm & speed==3),'.r','markersize',markersize);
    l=legend({'1 pix/frame','4 pix/frame', '8 pix/frame', '16 pix/frame'},'location','southeast');
    a = get(l,'children');
    set(a(1:3:end),'markersize',markersize);
    h=refline(1,0);
    set(h,'Color','k');
    title('Comparison of variation in motion estimation');
    xlabel('On parasol MAD');
    ylabel('On midget MAD');
    
end

cd('/Users/vision/Desktop/GitHub code repository/private/colleen');
if make_nearest_neighbor_histogram
    load_data_malcolm(1);
    for i = 1:4
        figure();
        d=nearest_neighbor_distances(i);
        hist(d,20);
        title({'Nearest neighbor distance',sprintf('2007-03-27-1 %s',cell_types{i})});
        xlabel('Distance');
        [h,x]=hist(d,20);
        x = x(h==max(h));
        lim = ylim; lim = lim(2);
        line([mean(d), mean(d)], [0, lim], 'Color', 'red');
        line([x, x], [0, lim], 'Color', 'green');
    end
    load_data_malcolm(2);
    for i = 1:3
        figure();
        d=nearest_neighbor_distances(i);
        hist(d,20);
        title({'Nearest neighbor distance',sprintf('2007-08-24-4 %s',cell_types{i})});
        xlabel('Distance');
        [h,x]=hist(d,20);
        x = x(h==max(h));
        lim = ylim; lim = lim(2);
        line([mean(d), mean(d)], [0, lim], 'Color', 'red');
        line([x, x], [0, lim], 'Color', 'green');
    end
end

% firing rate plot during moving bar stimuli
% need to have loaded data as in motion_script for this to work
if make_firing_rate_plot
    cell_indices = get_cell_indices(datarun{2},{run_opt.cell_type});
    N = length(cell_indices);
    firing_rates = nan(N,1);
    for i = 1:N
        firing_rates(i) = length(datarun{2}.spikes{cell_indices(i)})/1800;
    end
    histfit(firing_rates)
    title(sprintf('%s Run %d %s', run_opt.data_set, run_opt.data_run, run_opt.cell_type));
    xlabel('Hz');
    if save_figures
        saveas(gcf,sprintf('figures/%s/firing_rate_%s_run%d_config%d',run_opt.data_set,run_opt.cell_type,run_opt.data_run,run_opt.config_num),'png');
    end
end
    
% Make electrode map
if make_electrode_map
    [xCoords, yCoords] = getElectrodeCoords512;
    figure();
    xlim([min(xCoords)-10 max(xCoords)+10]);
    ylim([min(yCoords)-10 max(yCoords)+10]);
    for i=1:512
        text(xCoords(i),yCoords(i),num2str(i));
    end
    title('512 electrode map');
end

if load_results_downsample
    results = dataset([],[],[],[],[],[],[],'VarNames',...
    {'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'});

    cd('/Users/vision/Desktop/GitHub code repository/private/colleen/results/2007-03-27-1')
    nbins = 20;
    for my_run = 14:19

        if(my_run==15) 
            my_config = 3;
        elseif(my_run==17)
            my_config = 2;
        else
            my_config = 1;
        end

        if(my_run==14 || my_run==17)
            my_speed=1;
        elseif(my_run==15 || my_run==18)
            my_speed=2;
        else
            my_speed=3;
        end

        if my_run<19
            load(sprintf('On parasol_data_run_%d_config_%d', my_run, my_config));
            results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 1, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];

            load(sprintf('downsample_On midget_data_run_%d_config_%d', my_run, my_config));
            results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 3, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        else
            load(sprintf('downsample_On parasol_data_run_%d_config_%d', my_run, my_config));
            results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 1, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];

            load(sprintf('On midget_data_run_%d_config_%d', my_run, my_config));
            results = [results; dataset({'2007-03-27-1'}, my_run, my_speed, my_config, 3, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];
        end

    end

    cd('/Users/vision/Desktop/GitHub code repository/private/colleen/results/2007-08-24-4')
    for my_run = 4:11

        if(my_run==5 || my_run==10) 
            my_config = 2;
        elseif(my_run==6)
            my_config = 3;
        else
            my_config = 1;
        end

        if(my_run<6)
            my_speed=3;
        elseif(my_run<8)
            my_speed=2;
        elseif(my_run<10)
            my_speed=1;
        else
            my_speed=1.5;
        end

        load(sprintf('downsample_On parasol_data_run_%2.2d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-08-24-4'}, my_run, my_speed, my_config, 1, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];

        load(sprintf('On midget_data_run_%2.2d_config_%d', my_run, my_config));
        results = [results; dataset({'2007-08-24-4'}, my_run, my_speed, my_config, 3, mean(estimates), mad(estimates), 'VarNames',{'Experiment','Run','Speed','Config','CellType','MeanEstimate','MADEstimate'})];

    end

    sd = double(results(:,7));
    speed = double(results(:,3));
    onp = mod(1:28,2)';
    onm = ~onp;
    onp(1) = 0;
    onm(2) = 0;
    markersize=20;
    figure();
    loglog(sd(onp & speed==1),sd(onm & speed==1),'.m','markersize',markersize);
    hold on;
    loglog(sd(onp & speed==1.5),sd(onm & speed==1.5),'.b','markersize',markersize);
    loglog(sd(onp & speed==2),sd(onm & speed==2),'.g','markersize',markersize);
    loglog(sd(onp & speed==3),sd(onm & speed==3),'.r','markersize',markersize);
    l=legend({'1 pix/frame','4 pix/frame', '8 pix/frame', '16 pix/frame'},'location','southeast');
    a = get(l,'children');
    set(a(1:3:end),'markersize',markersize);
    h=refline(1,0);
    set(h,'Color','k');
    title({'Comparison of variation in motion estimation','(downsampled spikes)'});
    xlabel('On parasol MAD');
    ylabel('On midget MAD');

end
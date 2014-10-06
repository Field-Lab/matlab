function final_magic_number = choose_magic_number(datarun,bcf,bcf_params)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING STEP 4 OF 5    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% choose_magic_number     identify the magic number to use for saving cones
%
%
% usage:  final_magic_number = choose_msf(datarun,bcf)
%
% arguments:     datarun - datarun struct
%                    bcf - struct of results from bayesian cone finding
%
% outputs:     final_magic_number - the magic number to use
%
%
%
% 2010-02  gauthier
%



 

% get passed parameters
fn = fieldnames(bcf_params);
for ff=1:length(fn)
    eval(sprintf('%s = bcf_params.%s;',fn{ff},fn{ff}))
end

fn = fieldnames(bcf);
for ff=1:length(fn)
    eval(sprintf('%s = bcf.%s;',fn{ff},fn{ff}))
end






% PLOT WITH USER-SCALABLE MAGIC SCALE FACTOR


switch 1
    case 1
        % set magic scale factor
        the_magic = 'msf';
        magic_values = [1 2 5 10 15 20 25 30 40 50 70 85 100 120 140 160 200 240 300 380 460 600 900 1300];
    case 2
        % set magic offset
        mo = 35; the_magic = 'mo';
        magic_values = [0 5 10 15 20 40 70 100 200 400];
end

% make figure
figure(cones_fig);clf

% make slider
start_index=4; index_min=1; index_max=length(magic_values);
slider = make_loop_slider_list(start_index,index_min,index_max);
first_time = 1;

while 1

    % identify slider value
    k = round(get(slider,'Value'));
    
    % select axis
    subplot('Position',[.05 .05 .6 .9])

    % note xlim, ylim
    if first_time%isempty(get(gcf,'Child'))
        % if new figure, set to be relevant region
        curr_lim = [relevant_region_x; relevant_region_y];
        first_time = 0;
    else
        curr_lim = [xlim; ylim];
    end

    % plot dll
    imagesc((norm_image(   dll   )-0.5).^0.4,...
        'xdata',[1 datarun.stimulus.field_width]-(1-kernel_spacing),'ydata',[1 datarun.stimulus.field_height]-(1-kernel_spacing))
    axis image;hold on

    % plot previous algorithm cones
    if ~isempty(datarun.cones.centers)
        plot_cone_mosaic(datarun,'fig_or_axes',gca,'clear',0,'bg_color',[],'drawnow',false)
    end

    % identify cones to plot
    switch the_magic
        case 'msf'
            % get value based on slider
            msf = magic_values(k);
            % find which cones to plot
            plot_indices = (all_added_cones(:,7) + (msf/magic_number)*all_added_cones(:,8)) > 0;
            % note in title
            plot_title = sprintf('magic scale factor, %0.1f',msf);
        case 'mo'
            % get value based on slider
            mo = magic_values(k);
            % find which cones to plot
            plot_indices = (all_added_cones(:,7) + all_added_cones(:,8)) > mo;
            % note in title
            plot_title = sprintf('magic offset, %0.1f',mo);
    end

    % plot them
    for cc=1:size(kernel_plot_colors,1)
        cns = (all_added_cones(:,3) == cc) & plot_indices;
        plot(all_added_cones(cns,4),all_added_cones(cns,5),'o','Color',kernel_plot_colors(cc,:),'MarkerSize',12)
    end

    % plot ROI boundaries
    if 1
        for rr = 1:size(rois_x,1)
            roi_x = rois_x(rr,:) + [padding_x-1 -padding_x];
            roi_y = rois_y(rr,:) + [padding_y-1 -padding_y];
            plot([roi_x roi_x(2) roi_x(1) roi_x(1)],[roi_y(1) roi_y(1) roi_y(2) roi_y(2) roi_y(1)],'w')
        end
    end

    % restore bounds
    xlim(curr_lim(1,:))
    ylim(curr_lim(2,:))

    % add title
    title([datarun.names.nickname ', ' plot_title])
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plot DRPs
    
    % set up axes
    drp_plot_axes = subplot_axes(cones_fig,[.7 .1 .28 .9],.25,.23,1,4);
    
    % put data into convenient variables
    centers = all_added_cones(plot_indices,[4 5]);
    Lc = centers(all_added_cones(plot_indices,3)==1,:);
    Mc = centers(all_added_cones(plot_indices,3)==2,:);
    Sc = centers(all_added_cones(plot_indices,3)==3,:);
    
    
    %%%%%%%%%%%%%%%%%%%%%  DRP  %%%%%%%%%%%%%%%%%%%%%

    % parameters
    num_bins = 15;
    bin_width = 0.5;

    [drp,bin_centers,extras] = density_recovery_profile(centers,num_bins,bin_width);
    axes(drp_plot_axes{1});cla
    bar(bin_centers,drp)
    xlabel('distance (pixels)');ylabel('density (cones/pixel^2)')

    % add mean density, and effective radius
    hold on
    xlim_ = get(gca,'XLim');
    plot(xlim_,[1 1]*extras.density,'r',[1 1]*extras.eff_rad,[0 extras.density],'r')
    text(xlim_(2)*0.9,extras.density*1.1,sprintf('mean density: %0.3f\nradius: %0.1f',extras.density,extras.eff_rad),...
        'HorizontalAlignment','Right','VerticalAlignment','Bottom')

    % set ylim
    best_ylim = [0 max(drp)*1.3];
    set(gca,'YLim',best_ylim)


    
    %%%%%%%%%%%%%%%%%%%%%  DRP, grouped by cone type  %%%%%%%%%%%%%%%%%%%%%

    % from L cones
    [drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width);
    drp_M = density_recovery_profile(Mc,num_bins,bin_width,'reference_centers',Lc);
    drp_S = density_recovery_profile(Sc,num_bins,bin_width,'reference_centers',Lc);

    data = [drp_S; drp_M; drp_L]';

    axes(drp_plot_axes{2});cla
    bar(bin_centers(1:end),data(1:end,:),'stacked')
    set(gca,'YLim',best_ylim)
    xlabel('distance from L cones (pixels)');ylabel('density (cones/pixel^2)')
    d1=data;


    % from M cones
    [drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width,'reference_centers',Mc);
    drp_M = density_recovery_profile(Mc,num_bins,bin_width);
    drp_S = density_recovery_profile(Sc,num_bins,bin_width,'reference_centers',Mc);

    data = [drp_S; drp_M; drp_L]';

    axes(drp_plot_axes{3});cla
    bar(bin_centers(1:end),data(1:end,:),'stacked')
    set(gca,'YLim',best_ylim)
    xlabel('distance from M cones (pixels)');ylabel('density (cones/pixel^2)')
    d2=data;



    % from S cones
    [drp_L,bin_centers] = density_recovery_profile(Lc,num_bins,bin_width,'reference_centers',Sc);
    drp_M = density_recovery_profile(Mc,num_bins,bin_width,'reference_centers',Sc);
    drp_S = density_recovery_profile(Sc,num_bins,bin_width);

    data = [drp_S; drp_M; drp_L]';

    axes(drp_plot_axes{4});cla
    bar(bin_centers(1:end),data(1:end,:),'stacked')
    set(gca,'YLim',best_ylim)
    xlabel('distance from S cones (pixels)');ylabel('density (cones/pixel^2)')
    d3=data;


    uiwait;
end


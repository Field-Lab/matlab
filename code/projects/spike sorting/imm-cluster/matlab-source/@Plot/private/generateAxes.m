function current_axes = generateAxes(fig,num_clusters)
% GENERATEAXES current_axes = generateAxes(fig,num_clusters)
% Create 3 plot axes and acf axes.
%
% tamachado@salk.edu 1/31/08


figure(fig)
clf

set(gcf,'MenuBar','none')
current_axes.projection_axes{1}=subplot('Position',[.05 .69 .4 .30 ]);
current_axes.projection_axes{2}=subplot('Position',[.05 .36 .4 .30 ]);
current_axes.projection_axes{3}=subplot('Position',[.05 .03 .4 .30 ]);


if num_clusters > 1

    % parameters for setting up acf axis
    n_c = num_clusters;
    clust_left = 0.5;
    clust_right = 0.99;
    clust_top = .99;
    clust_bot = .03;
    cluster_plot_height = (clust_top - clust_bot)/ n_c;
    clust_width = clust_right - clust_left;


    for cc = 1: n_c
        %set up acf axis
        y_offset = clust_bot + (clust_top - clust_bot) * (n_c - cc)/n_c;
        current_axes.acf_axis{cc} = subplot('Position',[clust_left    y_offset clust_width cluster_plot_height ]);
    end

end
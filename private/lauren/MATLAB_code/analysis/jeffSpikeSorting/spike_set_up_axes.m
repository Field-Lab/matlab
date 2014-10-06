function current_axes = spike_set_up_axes(fig,num_clusters)

% create 3 plot axes and some cell stats: sta, average spike, and acf

num_proj_axes = 3;

switch num_clusters
    case 1
        figure(fig)
        clf
        %set(gcf,'MenuBar','none')

        current_axes.projection_axes{1}=subplot('Position',[.05 .75+.03 .4 .22 ]);
        current_axes.projection_axes{2}=subplot('Position',[.05 .5+.03 .4 .22 ]);
        current_axes.projection_axes{3}=subplot('Position',[.05 .25+.03 .4 .22 ]);
        current_axes.avg_spike_axis{1}=subplot('Position',[.05 .05 .9 0.2]);
        current_axes.acf_axis{1}=subplot('Position',[0.5+.03 .75 .45 .23]);
        current_axes.sta_axis{1}=subplot('Position',[0.5+.03 .55 .45 .15]);
        
    otherwise
        figure(fig)
        clf
        %set(gcf,'MenuBar','none')
        
        proj_left = 0.03;
        proj_right = 0.27;
        proj_top = 0.99;
        proj_bottom = 0.03;
        proj_gap = 0.04;
       
        proj_width = proj_right - proj_left;
        proj_height = (proj_top - proj_bottom + proj_gap) / num_proj_axes - proj_gap;
        
        for pp = 1:num_proj_axes
            current_axes.projection_axes{num_proj_axes - pp + 1}=subplot('Position',...
                [proj_left  (proj_height+proj_gap)*(pp-1)+proj_bottom  proj_width  proj_height ]);
        end
        
        % if too few clusters are selected, space the plots as if there
        % were more
        if num_clusters < 4
            n_c = 4;
        else
            n_c = num_clusters;
        end
        
        clust_left = 0.3;
        clust_right = 0.99;
        clust_gap = 0.02;
        
        cluster_plot_height = 0.8 * 1/n_c;
        clust_width = (clust_right - clust_left + clust_gap) / 3 - clust_gap;
        clust_spacing = clust_width + clust_gap;
        
        for cc = 1:num_clusters
            y_offset = (n_c - cc)/n_c;
            current_axes.avg_spike_axis{cc} = subplot('Position',[clust_left                  y_offset+.03 clust_width cluster_plot_height ]);
            current_axes.acf_axis{cc} =       subplot('Position',[clust_left+clust_spacing    y_offset+.03 clust_width cluster_plot_height ]);
            current_axes.sta_axis{cc} =       subplot('Position',[clust_left+2*clust_spacing  y_offset+.03 clust_width cluster_plot_height ]);
        end
        
end

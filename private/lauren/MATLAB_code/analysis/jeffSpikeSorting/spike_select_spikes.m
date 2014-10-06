function [selected_indices] = spike_pca_select_spikes(proj_struct,num_clusters_to_select,cluster_colors)

% select spikes of a projection structure (proj_struct) AFTER they have been plotted
% these fields are used:
%
% proj_struct.plot_axes.projection_axes  (axes where the projections are plotted)
%            .spike_projections  (rows = points, columns = dimensions)
%            .dimensions_to_plot  (which dimensions to put in each plot, e.g. [1 2;1 3;2 3] )
%            .dimension_prefix  (for the figure label, 'PC' for pca, 'LD' for lda, etc )
%            .density_histogram_bins  (bins in each dimension of the density plot)
%            .histogram_centers  (locations of the bins in the density histogram plot)


% load parameters from the struct
spike_projections = proj_struct.spike_projections;
projection_axes = proj_struct.plot_axes.projection_axes;
bins = proj_struct.density_histogram_bins;
dimensions_to_plot = proj_struct.dimensions_to_plot;
histogram_centers = proj_struct.histogram_centers;


% SELECT SUBSET OF SPIKES

for ss = 1:num_clusters_to_select
    
    %ask user which plot
    sp = input('which plot? (must be 1, 2, or 3) ');  %"sp" = which subplot

    %recall which pca dimensions it shows
    projection_first = dimensions_to_plot(sp,1);
    projection_second = dimensions_to_plot(sp,2);

    %have user draw a polygon in the selected plot
    axes(projection_axes{sp});
    [bw,xi,yi] = roipoly;

    % set up color for polygon
    cluster_color = cluster_colors(mod(ss - 1,length(cluster_colors))+1);
    
    %draw user's polygon there
    axes(projection_axes{sp});
    hold on;plot(xi,yi,cluster_color); hold off

    %change coordinates of polygon vertices to PC axes
    c = histogram_centers{sp};
    xv = xi/bins*range(c{1})+min(c{1});
    yv = -1*(yi/bins*range(c{2})-max(c{2}));

    %figure out which spikes were inside the polygon
    selected_indices{ss} = find(inpolygon(spike_projections(:,projection_first),spike_projections(:,projection_second),xv,yv))  ;

    %plot spike points (as a sanity check that the selected spikes were actually used)
    % figure(200);plot(spike_projections(selected_indices,projection_first),spike_projections(selected_indices,projection_second),'.')


end


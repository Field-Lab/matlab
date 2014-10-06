function proj_struct = spike_density_gradient(proj_struct)

% compute projection of a projection structure (proj_struct)
% these fields are used:
%
% proj_struct.plot_axes.projection_axes
%            .spike_projections  (rows = points, columns = dimensions)
%            .dimensions_to_plot  (which dimensions to put in each plot, e.g. [1 2;1 3;2 3] )
%            .dimension_prefix  (for the figure label, 'PC' for pca, 'LD' for lda, etc )
%            .density_histogram_bins  (bins in each dimension of the density plot)



density_plot_axes = proj_struct.plot_axes.projection_axes;
spike_proj = proj_struct.spike_projections;
dim_to_plot = proj_struct.dimensions_to_plot;
bins = proj_struct.density_histogram_bins;
dimension_prefix = proj_struct.dimension_prefix;


for ff=1:length(dim_to_plot)
    % note which dimensions to plot (e.g. PC 1 and PC 2)
    first_dimension = dim_to_plot(ff,1);
    second_dimension = dim_to_plot(ff,2);
    
    %generate density gradient, and save histogram centers
    %make 2d histogram of spike_locations
    [temp,hist_centers{ff}]=hist3([spike_proj(:,first_dimension) spike_proj(:,second_dimension)],[1 1]*bins);

    % save histogram maximum value
    density_gradient_norm_factor{ff} = max(max(temp));
  
    %normalize the values
    temp=temp/density_gradient_norm_factor{ff};
    %density_gradient=-1*(erfc(-50*temp)+4*temp);
    %density_gradient=-1*(erfc(-2*temp)+0*temp);
    %density_gradient=-1*(log(temp+0.01)+2*temp);

    density_gradient{ff}=-1*((log(0.2*temp+0.01)+2*temp+6).*(temp>0));
    %density_gradient=-1*((5*temp+2).*(temp>0));

    density_gradient{ff}=rot90(density_gradient{ff});
end

% save values in the structure
proj_struct.density_gradient = density_gradient;
proj_struct.density_gradient_norm_factor = density_gradient_norm_factor;
proj_struct.histogram_centers = hist_centers;

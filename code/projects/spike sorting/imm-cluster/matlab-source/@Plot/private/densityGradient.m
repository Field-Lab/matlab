function image = densityGradient(proj_struct, projection_axes)
%DENSITYGRADIENT Density Gradient Plot function
%   image = densityGradient(proj_struct, projection_axes) generates
%   a "density gradient" that to plot PC projections.
%   These plots are generated within projection_axes that can be used
%   to also hold clusters, STAs, and ACFs.
%
%   tamachado@salk.edu 1/23/08
%   from spike_plot_density_gradient (jgauthier)

density_plot_axes = projection_axes;
spike_proj = proj_struct.projections';
dim_to_plot = [1 2;1 3;2 3];
bins = 230;
dimension_prefix = 'PC';

density_gradient = cell(1, length(dim_to_plot));
hist_centers = cell(1, length(dim_to_plot));
density_gradient_norm_factor = cell(1, length(dim_to_plot));

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
    density_gradient{ff}=-1*((log(0.2*temp+0.01)+2*temp+6).*(temp>0));
    density_gradient{ff}=rot90(density_gradient{ff});

    %plot the density gradient
    d = density_gradient{ff};
    axes(density_plot_axes{ff});
    imagesc(d);
    colormap('bone');
    
    %remove tick marks
    set(density_plot_axes{ff},'YTick',[],'XTick',[]);
    
    % label axes
    xlabel(sprintf('%s %d',dimension_prefix,first_dimension));
    ylabel(sprintf('%s %d',dimension_prefix,second_dimension));
end

% save values in the structure
image.density_gradient = density_gradient;
image.density_gradient_norm_factor = density_gradient_norm_factor;
image.histogram_centers = hist_centers;

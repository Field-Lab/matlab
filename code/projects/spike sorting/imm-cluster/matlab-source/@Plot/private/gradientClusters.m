function image = gradientClusters(image, clusters)
%GRADIENTCLUSTERS Add clusters and ACFs to cluster plot
%   image = gradientClusters(image, clusters) adds clusters to a
%   projections plot and displays them using a convex hull and contours.
%   This code is directly ported over from spike_show_same_spikes_fancy
%   which was written by Jeff Gauthier. Only minor modifications to make it
%   work in this new framework were made.
%   tamachado@salk.edu 1/29/08
%   from spike_show_same_spikes_fancy (jgauthier)

proj_struct.density_plot_axes = clusters.axes;
proj_struct.spike_projections = clusters.projections;
proj_struct.dimensions_to_plot = [1 2;1 3;2 3];
proj_struct.density_histogram_bins = 230;
proj_struct.dimension_prefix = 'PC';

proj_struct.density_gradient = image.density_gradient;
proj_struct.density_gradient_norm_factor = image.density_gradient_norm_factor ;
proj_struct.histogram_centers = image.histogram_centers;
proj_struct.selected_indices = clusters.selected_indices;
proj_struct.spike_times = clusters.spike_times;

proj_struct_orig = proj_struct;
proj_struct_new  = proj_struct;
clear proj_struct;

try
    % show spikes selected using one projection (new) when plotted in another
    % projection (original)

    proj_struct_orig.spike_projections = double(proj_struct_orig.spike_projections);
    proj_struct_new.spike_projections = double(proj_struct_new.spike_projections);

    which_selected_indices = 1:length(proj_struct_orig.selected_indices);
    
    
    % set up colors to use for various clusters
    cluster_colors = [
        1.0000         0         0;
        0    1.0000         0;
        0         0    1.0000;
        0    1.0000    1.0000;
        1.0000         0    1.0000;
        1.0000    1.0000         0;
        0         0         0;
        1.0000    0.8000    0.8000;
        0.8000    1.0000    0.8000;
        0.8000    0.8000    1.0000;
        0.6000    0.8000    1.0000;
        1.0000    0.8000    1.0000;
        1.0000    0.8000         0;
        0.5020    0.5020    0.5020;
        0.8471    0.1608         0;
        0    0.4980         0;
        0.0431    0.5176    0.7804;];


    % GET PARAMETERS OF NEW PROJECTION PLOT
    dim_to_plot = proj_struct_new.dimensions_to_plot;
    bins = 230;
    bin_centers = proj_struct_new.histogram_centers;
    dimension_prefix = proj_struct_new.dimension_prefix;
    acf_range = 30;
    cluster_display = 'convex hull'; % works for any number of clusters

    % SET UP AXES
    current_axes.acf_axis = clusters.acf_axis;
    current_axes.projection_axes = clusters.axes;
    proj_struct_new.plot_axes = clusters.axes;
    
    
    new_dg = cell(1, length(dim_to_plot));
    
    % START DENSITY PLOT
    for ff = 1:length(dim_to_plot)
        % use old gradient
        new_dg{ff} = zeros(bins,bins,3);
        for cc = 1:3
            new_dg{ff}(:,:,cc) = proj_struct_new.density_gradient{ff};
        end
    end

    % note which selected_indices don't exist in this projection
    dont_plot = zeros(1,length(which_selected_indices));  

    % GO THROUGH EACH CLUSTER (I.E. EACH SET OF SELECTED INDICES)
    for ss = 1:length(which_selected_indices)

        si = which_selected_indices(ss);

        % IDENTIFY SELECTED SPIKES IN THE ORIGINAL PROJECTION

        % note the spike times from the original projection
        orig_spike_times = double(proj_struct_orig.spike_times(proj_struct_orig.selected_indices{si}));

        % i_new is the list of indices in the new projection
        [foo i_orig i_new] = intersect(orig_spike_times, proj_struct_new.spike_times);
        
        % save them in the new projection
        proj_struct_new.selected_indices{ss} = i_new;

        % PLOT THE SPIKES IN THE NEW PROJECTION

        %make sure there is at least a miniumum number of spikes
        if (isempty(i_new) || length(i_new) < 100)
            dont_plot(ss) = 1;
        else

            % set the color for this cluster
            col = cluster_colors(mod(ss - 1,length(cluster_colors))+1,:);
            
            
            hull_points_x = cell(1, length(dim_to_plot));
            hull_points_y = cell(1, length(dim_to_plot));
            
            % density gradient plot
            for ff=1:length(dim_to_plot)

                % note which dimensions to plot (e.g. PC 1 and PC 2)
                first_dimension = dim_to_plot(ff,1);
                second_dimension = dim_to_plot(ff,2);

                %generate density gradient, and save histogram centers
                temp =hist3([proj_struct_new.spike_projections(i_new,first_dimension)...
                    proj_struct_new.spike_projections(i_new,second_dimension)],...
                    bin_centers{ff});
                
                temp = double(temp);

                % normalize the values
                temp=temp/proj_struct_new.density_gradient_norm_factor{ff};

                % replot on log scale
                density_gradient=-1*((log(0.2*temp+0.01)+2*temp+6).*(temp>0));

                % rotate by 90 degrees
                density_gradient=rot90(density_gradient);

                % compute the convex hull
                x = proj_struct_new.spike_projections(i_new,first_dimension);
                y = proj_struct_new.spike_projections(i_new,second_dimension);
                k = convhull(x,y);
                hull_points_x{ff}{ss}=x(k);
                hull_points_y{ff}{ss}=y(k);

                contours = .6;
                contour_line_widths = 4;

                % compute contour lines at this threshold
                cnt = contourc(double(bin_centers{ff}{1}),double(bin_centers{ff}{2}),double(-flipud(density_gradient)/max(max(-density_gradient))),double([contours contours]));

                % set up temporary variables to store contour line vertices
                contour_bands = find((cnt(1,:) == contours));
                contour_bands = [contour_bands length(cnt)+1];

                contour_union = [];

                % for each contour band at this threshold, save the vertices as a new item in contour_bands_observed
                for ii = 1:length(contour_bands)-1
                    % get the contours as observed
                    contour_bands_observed{si}{ff}{ii} = ...
                        [ cnt(1,contour_bands(ii)+1:contour_bands(ii+1)-1);...
                        cnt(2,contour_bands(ii)+1:contour_bands(ii+1)-1)];
                    % compute the union of points at this contour level
                    contour_union = [contour_union contour_bands_observed{si}{ff}{ii}];
                end

                % save the complex hull
                x = contour_union(1,:);
                y = contour_union(2,:);
                k = convhull(x,y);
                hull_contour_points_x{si}{ff}=x(k);
                hull_contour_points_y{si}{ff}=y(k);
            end

            % note the spike times and waveforms of this cluster
            selected_spike_times = orig_spike_times(i_orig);

            % plot acf and average spike
            computeACF(selected_spike_times,acf_range,current_axes.acf_axis{ss},col);
            
            %remove tick marks
            for ii = 1:length(current_axes.acf_axis)-1
                set(current_axes.acf_axis{ii},'YTick',[],'XTick',[]);
            end
            
            set(current_axes.acf_axis{length(current_axes.acf_axis)},'YTick',[]);

            %spike_plot_average_spike(selected_spike_waveforms,length(proj_struct_new.electrodes_for_projections),...
            %   proj_struct_new.window_length,current_axes.avg_spike_axis{ss},[0 0 0]);
        end

    end

    
    
    % go through each projection plot
    for ff = 1:length(dim_to_plot)

        % normalize the density gradient
        new_dg{ff} = new_dg{ff} - min(min(min(new_dg{ff})));
        new_dg{ff} = new_dg{ff} / max(max(max(new_dg{ff})));

        % plot it
        axes(current_axes.projection_axes{ff});
        ctrs=proj_struct_new.histogram_centers{ff};
        imagesc(new_dg{ff},'XData',[min(ctrs{1}) max(ctrs{1})],'YData',[max(ctrs{2}) min(ctrs{2})]);hold on;
        axis xy;

        % remove tick marks
        set(current_axes.projection_axes{ff},'YTick',[],'XTick',[])

        % label axes
        xlabel(sprintf('%s %d',dimension_prefix,dim_to_plot(ff,1)));
        ylabel(sprintf('%s %d',dimension_prefix,dim_to_plot(ff,2)))

        % add convex hull plot
        if strcmp(cluster_display,'convex hull')
            for ss = 1:length(which_selected_indices)
                if dont_plot(ss)

                else
                    % set the color for this cluster
                    col = cluster_colors(mod(ss - 1,length(cluster_colors))+1,:);

                    % plot convex hull
                    plot(current_axes.projection_axes{ff},hull_points_x{ff}{ss},hull_points_y{ff}{ss},'Color', col, 'LineWidth', 2.5);

                    % plot hull of contours
                    plot(current_axes.projection_axes{ff},hull_contour_points_x{ss}{ff},hull_contour_points_y{ss}{ff},...
                        'Color',col,'LineWidth',contour_line_widths);
                end
            end
        end
    end
catch
    disp(lasterror)
end
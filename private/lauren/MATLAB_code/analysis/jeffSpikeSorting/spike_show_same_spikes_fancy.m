function proj_struct_new = spike_show_same_spikes_fancy(proj_struct_orig,which_selected_indices,proj_struct_new,figure_to_use)
try
    % show spikes selected using one projection (new) when plotted in another
    % projection (original)

    proj_struct_orig.spike_projections = double(proj_struct_orig.spike_projections);
    proj_struct_new.spike_projections = double(proj_struct_new.spike_projections);


    % default to showing all selected indices
    if isempty(which_selected_indices)
        which_selected_indices = 1:length(proj_struct_orig.selected_indices);
    end

    % set up colors to use for various clusters
    %cluster_colors = 'rgbm'; cluster_colors = cluster_colors';
    cluster_colors = [0.2 0.2 0.2;0 1 0; 1 0.5 0; 0 0 1; 0 1 1;1 0 1; 1 .8 .8; 1 1 0];
    %cluster_colors = [1 0 1; 1 .8 .8; 1 1 0];
    acf_range = 30;

    %cluster_display = 'gradient'; % works only for three or fewer clusters
    cluster_display = 'convex hull'; % works for any number of clusters



    % GET PARAMETERS OF NEW PROJECTION PLOT
    dim_to_plot = proj_struct_new.dimensions_to_plot;
    bins = proj_struct_new.density_histogram_bins;
    bin_centers = proj_struct_new.histogram_centers;
    dimension_prefix = proj_struct_new.dimension_prefix;

    % SET UP AXES
    current_axes = spike_set_up_axes(figure_to_use,length(which_selected_indices));
    proj_struct_new.plot_axes = current_axes;

    % MAKE SCATTER PLOT
    % figure(figure_to_use + 1);clf
    % plot(proj_struct_new.spike_projections(:,1),proj_struct_new.spike_projections(:,2),'k.')

    % START DENSTIY PLOT
    for ff = 1:length(dim_to_plot)
        switch 1
            case 1
                % use old gradient
                new_dg{ff} = zeros(bins,bins,3);
                for cc = 1:3
                    new_dg{ff}(:,:,cc) = proj_struct_new.density_gradient{ff};
                end
            case 2
                % compute new density gradient
                [new_dg{ff},bin_centers{ff}]=hist3([proj_struct_new.spike_projections(:,dim_to_plot(ff,1))...
                    proj_struct_new.spike_projections(:,dim_to_plot(ff,1))],[bins bins]);
        end
    end

    % note which selected_indices don't exist in this projection
    dont_plot = zeros(1,length(which_selected_indices));

    % clear current selected indices from the new projection
    proj_struct_new.selected_indices = [];

    % GO THROUGH EACH CLUSTER (I.E. EACH SET OF SELECTED INDICES)
    for ss = 1:length(which_selected_indices)

        si = which_selected_indices(ss);

        % IDENTIFY SELECTED SPIKES IN THE ORIGINAL PROJECTION

        % note the spike times from the original projection
        orig_spike_times = double(proj_struct_orig.spike_times(proj_struct_orig.selected_indices{si}));
        [foo i_orig i_new] = intersect(orig_spike_times, proj_struct_new.spike_times);
        % i_new is the list of indices in the new projection

        % save them in the new projection
        proj_struct_new.selected_indices{ss} = i_new;

        % PLOT THE SPIKES IN THE NEW PROJECTION
        
        %make sure there is at least a miniumum number of spikes
        if (isempty(i_new) || length(i_new) < 100)
            dont_plot(ss) = 1;
        else

            % set the color for this cluster
            col = cluster_colors(mod(ss - 1,length(cluster_colors))+1,:);

            % 1. scatter plot
            %     figure(figure_to_use + 1);
            %     hold on;
            %     plot(proj_struct_new.spike_projections(i_new,1),proj_struct_new.spike_projections(i_new,2),'.','color',col);
            %     hold off;


            % 2. density gradient plot

            for ff=1:length(dim_to_plot)

                % note which dimensions to plot (e.g. PC 1 and PC 2)
                first_dimension = dim_to_plot(ff,1);
                second_dimension = dim_to_plot(ff,2);

                % generate the denstiy gradient for this cluster

                %generate density gradient, and save histogram centers
                %make 2d histogram of spike_locations
                [temp,nil]=hist3([proj_struct_new.spike_projections(i_new,first_dimension)...
                    proj_struct_new.spike_projections(i_new,second_dimension)],...
                    bin_centers{ff});

                temp = double(temp);

                % normalize the values
                temp=temp/proj_struct_new.density_gradient_norm_factor{ff};

                % replot on log scale
                density_gradient=-1*((log(0.2*temp+0.01)+2*temp+6).*(temp>0));

                % rotate by 90 degrees
                density_gradient=rot90(density_gradient);

                switch cluster_display
                    case 'gradient'

                        % add to the density gradient
                        for cc = 1:3
                            new_dg{ff}(:,:,cc) = new_dg{ff}(:,:,cc) + density_gradient * (- col(cc));
                        end

                    case 'convex hull'
                        % save density gradient
                        %new_dg{ff} = density_gradient;

                        % compute the convex hull
                        x = proj_struct_new.spike_projections(i_new,first_dimension);
                        y = proj_struct_new.spike_projections(i_new,second_dimension);
                        k = convhull(x,y);
                        hull_points_x{ff}{ss}=x(k);
                        hull_points_y{ff}{ss}=y(k);
                        [x(k) y(k)];

                        contours = [0.5 0.8];
                        contour_line_widths = [1 2];
                        num_contours = length(contours);

                        for tt = 1:num_contours

                            % compute contour lines at this threshold
                            cnt = contourc(double(bin_centers{ff}{1}),double(bin_centers{ff}{2}),double(-flipud(density_gradient)/max(max(-density_gradient))),double([contours(tt) contours(tt)]));

                            % set up temporary variables to store contour line vertices
                            contour_bands = find((cnt(1,:) == contours(tt)));
                            contour_bands = [contour_bands length(cnt)+1];

                            contour_union = [];

                            % for each contour band at this threshold, save the vertices as a new item in contour_bands_observed
                            for ii = 1:length(contour_bands)-1
                                % get the contours as observed
                                contour_bands_observed{si}{ff}{tt}{ii} = ...
                                    [ cnt(1,contour_bands(ii)+1:contour_bands(ii+1)-1);...
                                    cnt(2,contour_bands(ii)+1:contour_bands(ii+1)-1)];
                                % compute the union of points at this contour level
                                contour_union = [contour_union contour_bands_observed{si}{ff}{tt}{ii}];
                            end

                            % save the complex hull
                            x = contour_union(1,:);
                            y = contour_union(2,:);
                            k = convhull(x,y);
                            hull_contour_points_x{si}{ff}{tt}=x(k);
                            hull_contour_points_y{si}{ff}{tt}=y(k);

                        end

                        %disp('here')
                end
            end

            % note the spike times and waveforms of this cluster
            selected_spike_times = orig_spike_times(i_orig);
            selected_spike_waveforms = proj_struct_new.spike_waveforms(i_new,:);

            % plot acf and average spike
            spike_compute_and_plot_acf(selected_spike_times,acf_range,current_axes.acf_axis{ss},col);

            spike_plot_average_spike(selected_spike_waveforms,length(proj_struct_new.electrodes_for_projections),...
                proj_struct_new.window_length,current_axes.avg_spike_axis{ss},[0 0 0]);

            % plot STAs if they're calculated
            if isfield(proj_struct_new,'sta')
                if length(proj_struct_new.sta) >= si
                    if ~isempty(proj_struct_new.sta{si})
                        spike_plot_sta(proj_struct_new.sta{si},proj_struct_new.plot_axes.sta_axis{si});
                    end
                end
            end

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
                    plot(current_axes.projection_axes{ff},hull_points_x{ff}{ss},hull_points_y{ff}{ss},'Color',col);

                    % plot contours
                    %             for tt = 1:num_contours
                    %                 for ii = 1:length(contour_bands_observed{ss}{ff}{tt})
                    %                     plot(current_axes.projection_axes{ff},contour_bands_observed{ss}{ff}{tt}{ii}(1,:),contour_bands_observed{ss}{ff}{tt}{ii}(2,:),...
                    %                         'Color',col,'LineWidth',contour_line_widths(tt))
                    %                 end
                    %             end

                    % plot hull of contours
                    for tt = 1:num_contours
                        plot(current_axes.projection_axes{ff},hull_contour_points_x{ss}{ff}{tt},hull_contour_points_y{ss}{ff}{tt},...
                            'Color',col,'LineWidth',contour_line_widths(tt));
                    end
                end

            end
        end

    end
catch
    disp(lasterr)

end




% for ss = 1:length(which_selected_indices)
%     figure(502);hold on;
%     plot(hull_points_x{1}{ss},hull_points_y{1}{ss},'Color','b');
% end
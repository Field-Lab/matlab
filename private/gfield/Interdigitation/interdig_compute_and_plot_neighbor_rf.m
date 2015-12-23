function profiles = interdig_compute_and_plot_neighbor_rf(datarun,cell_types, varargin)
% plot average neighbor interaction


p = inputParser;
p.addParamValue('center_point_type', 'com') % com, obvius_fit, or vision_fit
p.addParamValue('profile_norm_type', 'norm') % 'norm' or 'no norm'
p.addParamValue('foa', 1)
p.addParamValue('clear', true)
p.addParamValue('colors', 'krbg')
p.addParamValue('profile_points', 100)
p.addParamValue('extension_factor', 2.0)
p.addParamValue('distance_cutoff_factor', 1.25)


p.parse(varargin{:})

params = p.Results;
colors = params.colors;
profile_points = params.profile_points;
extension_factor = params.extension_factor;
distance_cutoff_factor = params.distance_cutoff_factor;



% set up plot axes
plot_axes = set_up_fig_or_axes(params.foa,params.clear);
axes(plot_axes)
%set(1,'Color','w');



switch params.profile_norm_type
    case 'norm'
        normalize_average = 1;
    case 'no norm'
        normalize_average = 0;
end

plot_axes = set_up_fig_or_axes(params.foa,params.clear);


% plot for each cell type
for tt = 1:length(cell_types)

    % GATHER THE NEEDED DATA
    
    % get the cell type
    cell_type = cell_types(tt);
    % get the color
    color = colors(mod(tt - 1,length(colors))+1);

    % clear variables which store results
    clear avg_profile center_pts neighbor_distance

    % get center locations of all cells in this mosaic

    center_pts = rf_centers(datarun, cell_type, params.center_point_type);
    %for cc = 1:length(stas.cell_types.(cell_type))
        %center_pts(cc,1:2) = interdig_get_center_point(stas,center_point_type,stas.cell_types.(cell_type)(cc));
    %end
    
    % go through each cell and accumulate parameters in a local variable
    cell_indices = get_cell_indices(datarun, cell_type);
    num_rgcs = length(cell_indices);
    
    for cc = 1:num_rgcs

        % establish cell id number
        cell_id = datarun.cell_types{cell_type{:}}.cell_ids(cc);

        clear distances
        
        % compute distance to all neighbors
        distances = ipdm(center_pts, center_pts(cc,:));
        
        % set the self-distance to a large number
        distances(find(distances==0)) = max(distances);
        
        % find the cell number of the closest neighbor
        neighbor_distance(cc) = min(distances);
        neighbor_cell_id = datarun.cell_types{cell_type{:}}.cell_ids( find(distances==neighbor_distance(cc)));
        neighbor_cell_id = neighbor_cell_id(1);
        neighbor_index = get_cell_indices(datarun, neighbor_cell_id);
        
        % get center points for the cell and its nearest neighbor
        cell_ctr = rf_center(datarun, cell_id, 'com', 'sta');
        nbr_ctr = rf_center(datarun, neighbor_cell_id, 'com', 'sta');
        % make sure the center points are not empty, return error message if they are
        if isempty(cell_ctr) || isempty(nbr_ctr)
            if isempty(cell_ctr)
                temp_message = ['cell', num2str(cell_id), ' does not have a center point'];
            else
                temp_message = ['cell', num2str(neighbor_cell_ids), ' does not have a center point'];
            end
            error(temp_message)
        end
        
        % compute the matrix slices from each cell
        
        temp_rf = datarun.stas.rfs{cell_indices(cc)};
        temp_rf_nbr = datarun.stas.rfs{neighbor_index};
        if size(temp_rf,3) == 3
            temp_rf = temp_rf(:,:,2);
            temp_rf_nbr = temp_rf_nbr(:,:,2);
        end

        cell_profile{cc}{1} = matrix_slice(temp_rf, [cell_ctr; nbr_ctr], 'extension', extension_factor, 'samples', profile_points);
        cell_profile{cc}{2} = matrix_slice(temp_rf_nbr, [cell_ctr; nbr_ctr], 'extension', extension_factor, 'samples', profile_points);

        if 0 % show each cell
            figure(2)
            plot(1:profile_points,cell_profile{cc}{1},1:profile_points,cell_profile{cc}{2})
            pause
        end
        
    end
    
    maximum_neighbor_distance = distance_cutoff_factor*mean(neighbor_distance);
    
    % toss out cells which are beyond the maximum neighbor distance
    avg_profile{1} = zeros(profile_points,1);
    avg_profile{2} = zeros(profile_points,1);
    cell_counter = 0;
    for cc = 1:num_rgcs
        if neighbor_distance(cc) < maximum_neighbor_distance
            cell_counter = cell_counter +1;

            profile_a = cell_profile{cc}{1};
            nan_indices = find(isnan(profile_a));
            profile_a(nan_indices) = 0;

            profile_b = cell_profile{cc}{2};
            nan_indices = find(isnan(profile_b));
            profile_b(nan_indices) = 0;
 
            %size(profile_a)
            %size(avg_profile{1})
            avg_profile{1} = avg_profile{1} + profile_a;
            avg_profile{2} = avg_profile{2} + profile_b;
        end
    end
    
    % divide by number of cells
    avg_profile{1} = avg_profile{1} / cell_counter;
    avg_profile{2} = avg_profile{2} / cell_counter;
    
    % normalize variance
    if params.profile_norm_type
        avg_profile{1} = avg_profile{1} / std(avg_profile{1});
        avg_profile{2} = avg_profile{2} / std(avg_profile{2});
    end
    

    % PLOT FOR OVERLAY
    %figure(1); hold on
%plot([1:profile_points]/profile_points*(2*extension_factor+1),[avg_profile{2}],'Color',color);
    plot([1:profile_points]/profile_points*(2*extension_factor+1),[avg_profile{2} - mean(avg_profile{2}(1:10))],'Color',color);
    hold on
%plot([1:profile_points]/profile_points*(2*extension_factor+1),[avg_profile{1}],'Color',color);
    plot([1:profile_points]/profile_points*(2*extension_factor+1),[avg_profile{1} - mean(avg_profile{1}(1:10))],'Color',color);

    % change tick marks
    %set(gca,'TickDir','out');
    
    temp_name = datarun.cell_types{cell_type{:}}.name;
    temp_name = regexprep(temp_name, ' ', '_');  % replace a space with an underscore

    profiles.([temp_name, '_x']){1} = [1:profile_points]/profile_points*(2*extension_factor+1);
    profiles.([temp_name, '_y']){1} = avg_profile{1};
    profiles.([temp_name, '_x']){2} = [1:profile_points]/profile_points*(2*extension_factor+1);
    profiles.([temp_name, '_y']){2} = avg_profile{2};
    
end

% add black lines
%figure(1);
plot([0 2*extension_factor+1],[0 0],'k')
plot([extension_factor extension_factor],get(gca,'YLim'),'k')
plot([extension_factor+1 extension_factor+1],get(gca,'YLim'),'k')


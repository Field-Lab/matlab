% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = 'plantain';

% apricot
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2009-04-13-5/data005-s3600-s7200/data005/data005';
path_and_name{1,2} = 'apricot';

% peach
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-0/data001-0s-2400s/data001/data001';
path_and_name{1,2} = 'peach';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

%%

cell_type = 4;
cell_indices = get_cell_indices(datarun, {cell_type});
threshold = 0.1;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', threshold,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'U');   

connectivity = mosaic_weights .* selection; 
% the new datarun has S cone excluded from all fields of datarun.cone 
new_datarun = extras.new_datarun;
[num_cones, num_RGCs] = size(connectivity);


% cycle through cells and make convex hull for each
for rgc = 1:num_RGCs
    % get cone indices
    cone_indices = find(connectivity(:,rgc));
    % make sure there are enough cones to triangulate
    if length(cone_indices) > 2
        % get cone locations
        locations = new_datarun.cones.centers(cone_indices,:);
        % delauney triangulate the points
        temp_dt = DelaunayTri(locations(:,1), locations(:,2));
        % compute the convex hull
        temp_hull = convexHull(temp_dt);

        % store the hull in the cones field
        datarun.cones.convexhull{cell_indices(rgc)} = temp_hull;
        
%         figure
%         plot(locations(:,1), locations(:,2), 'r.');
%         hold on
%         plot(locations(temp_hull,1), locations(temp_hull,2), 'k-')
%         hold off
        
    end
end


%%

% compute purity indices
purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);

% get indices of cells with balanced L/M input
mean_purity = median(purity_indices);
relative_purities = purity_indices - mean_purity;
[s_rel_purities, sorted_pur_indices] = sort(abs(relative_purities), 'ascend'); 
balanced_indices = sorted_pur_indices(1:10);
%balanced_indices = find(purity_indices > -0.05 & purity_indices < 0.05);
[junk, ascend_list] = sort(purity_indices, 'ascend');
M_dom_indices = ascend_list(1:10);

[junk, descend_list] = sort(purity_indices, 'descend');
L_dom_indices = descend_list(1:10);

all_indices = [balanced_indices, M_dom_indices, L_dom_indices];

% get  indices of cells with mostly L input
%L_dom_indices = find(purity_indices > 0.7);
% get indicex of cells with mostly M input
%M_dom_indices = find(purity_indices < -0.5);



indices_of_interest = all_indices;

% do some plotting to visually compare these cells with their convex hulls

rgb = [0.5 0.5 0.5];
type_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0];
panel_num = ceil(sqrt(length(indices_of_interest)));

% generate background matrix
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));


for cll = 1:length(indices_of_interest)

    
    %subplot(panel_num, panel_num, cll)
    figure('Position', [1 1 1000 1000])
    image(plot_mat);axis image; hold on; axis off; axis square
    

    cell_id = datarun.cell_types{cell_type}.cell_ids(indices_of_interest(cll));
    cell_index = get_cell_indices(datarun,cell_id);


    % plot_center and surround spiders
    plot_cell_sampling(datarun,cell_id,'type','spider','fig_or_axes',-1,'clear',0,'plot_cones',0,...
       'line_width',[realmin 2.5],'cell_colors',[1 1 1],'plot_radius',[0 inf],'thresh',threshold,'polarity',1,'contiguity',true)
   
    % add white circles under connected cones
    % plot center cone halos
    cone_colors = [1 1 1]; 
    cone_colors = repmat(cone_colors,4,1);
    selection_params.radius = [0 inf];
    selection_params.contiguity = true;
    selection_params.thresh = threshold;
    selection_params.polarity = 1;
    [weights, center_selection, extras] = select_cone_weights(datarun, cell_id, selection_params);

    plot_cone_mosaic(datarun,'fig_or_axes',-1,'cone_size',15,'clear',0,'bg_color',[],...
                    'cone_roi', center_selection, 'roi_highlight', false, 'type_colors', cone_colors)

              
    % plot all cones
    plot_cone_mosaic(datarun,'fig_or_axes',-1,'cone_size',10,'clear',0,'bg_color',[], 'type_colors', type_colors)

    % plot convex hull of the cell
    temp_index = indices_of_interest(cll);
    temp_cone_indices = find(connectivity(:,temp_index));
    convergence = length(temp_cone_indices);
    temp_locations = new_datarun.cones.centers(temp_cone_indices,:);
    if convergence > 2
        temp_hull = datarun.cones.convexhull{cell_index};
        plot(temp_locations(temp_hull,1), temp_locations(temp_hull,2), 'k-');
    end
    
    
    % re-frame the image
    ws = 16;
    center_pt = datarun.cones.rf_fits{cell_index}.center;
    bg_frame_x = center_pt(1) - ws;
    end_frame_x = center_pt(1) + ws;
    bg_frame_y = center_pt(2) - ws;
    end_frame_y = center_pt(2) + ws;
    axis([bg_frame_x end_frame_x bg_frame_y end_frame_y])
    temp_fig_title = num2str(purity_indices(indices_of_interest(cll)));
    title(temp_fig_title)

end












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

%% Compute the purity indices from the data.
cell_type = 4;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0,...
                                            'remove_cones', 'S');   

connectivity = mosaic_weights .* selection; % keep weights continuous valued

% the new datarun has S cone excluded from all fields of datarun.cone 
new_datarun = extras.new_datarun;

[num_cones, num_RGCs] = size(connectivity);


% normalize cone weights
for RGC = 1:num_RGCs
    connectivity(:,RGC) = connectivity(:,RGC) ./ sum(connectivity(:,RGC));
end

% compute purity of cone mosaic
purity_indices = compute_opponency_index(connectivity, new_datarun.cones.types);
finited_indices = isfinite(purity_indices);
purity_sd = std(purity_indices(finited_indices));

[sorted_junk, invest_indices] = sort(abs(purity_indices), 'descend');
invest_indices = invest_indices(1:16);

%%

% plot it
figure;
% plot background
% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [0.5 0.5 0.5];
type_colors = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));


for cll = 1:16

    subplot(4,4,cll)
    image(plot_mat);axis image; hold on; axis off; axis square
    

    cell_id = datarun.cell_types{cell_type}.cell_ids(invest_indices(cll));
    cell_index = get_cell_indices(datarun,cell_id);


    % plot_center and surround spiders
    plot_cell_sampling(datarun,cell_id,'type','spider','fig_or_axes',-1,'clear',0,'plot_cones',0,...
       'line_width',[realmin 2.5],'cell_colors',[1 1 1],'plot_radius',[0 6],'thresh',0.1,'polarity',1,'contiguity',true)
   
    % plot all cones
    plot_cone_mosaic(datarun,'fig_or_axes',-1,'cone_size',10,'clear',0,'bg_color',[], 'type_colors', type_colors)

    % re-frame the image
    ws = 15;
    center_pt = datarun.cones.rf_fits{cell_index}.center;
    bg_frame_x = center_pt(1) - ws;
    end_frame_x = center_pt(1) + ws;
    bg_frame_y = center_pt(2) - ws;
    end_frame_y = center_pt(2) + ws;
    axis([bg_frame_x end_frame_x bg_frame_y end_frame_y])

end




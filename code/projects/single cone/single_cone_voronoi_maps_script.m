% specify the directory to write to
date_piece = '2010-10-18-1';
datarun = 'data003';
mosaic = 'on-midget';
full_map_folder = 'all-cones';
write_dir = '/snle/acquisition/maps/single-cone/';
cd(write_dir)
cone_folder = [date_piece,'-',datarun,'-',mosaic,];
mkdir([date_piece,datarun,'-',full_map_folder])
write_path = [write_dir, cone_folder,'/'];


% define data path
%datarun = load_data('/snle/acquisition/2010-09-24-1/data035/data035');
datarun = load_data('/snle/acquisition/2011-01-11-2/data005/data005');
path_and_name{1,2} = '_snle_acquisition_2011-01-11-2_data005_data005-bayes-msf_15.00-2011-01-11-2_data005_withparasols2';

% alternative data path
datarun = load_data('/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006');
path_and_name{1,2} = '2010-09-24-1_streamed_data006_data006-bayes-msf_15.00-15-foo';


% load data
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

width = datarun.stimulus.field_width;
height = datarun.stimulus.field_height;

%% make tessalation
centers = datarun.cones.centers;
[V,C] = voronoin(centers);

%% make cone map
verbose = false;
clear cone_map
cone_map = zeros(width, height);
figure(101); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    clear temp_mask
    if all(C{cn} ~=1)
        temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
        if length(find(temp_mask ==1)) <= 40
            if verbose
                patch(V(C{cn},1),V(C{cn},2), 'r')
                axis([1 320 1 320])
                axis square
                drawnow
                pause(0.01)
            end
            new_counter = new_counter +1;
            cone_map = cone_map + (new_counter * temp_mask);
        end
    end
end
% transpose cone map to get placement correct
write_cone_map = cone_map';
dlmwrite([write_dir, date_piece,'-', datarun,'-all-cones'], write_cone_map, 'delimiter', '\t', 'newline', 'pc')

%figure(2)
%imagesc(cone_map)
%print(101, '~/Desktop/cone_map_orange.pdf','-dpdf')


%% get the strongest cone feeding each cell of a given type

% define cell type of interest
cell_type = 4;


% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% get indices to strongest cones
[num_cones, num_rgcs] = size(connectivity);
peak_cone_indices = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    [max_val, temp_cn_index] = max(connectivity(:,cc));
    peak_cone_indices(cc) = temp_cn_index;
end

% make and plot the map of the strongest cones
peak_cone_map = zeros(width, height);
verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(peak_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                new_counter = new_counter +1;
                peak_cone_map = peak_cone_map + temp_mask;
                if new_counter == 2
                    continue
                end
                
            end
        end
    end
end
peak_cone_map = peak_cone_map';
dlmwrite([cone_folder,'/map-0000.txt'], peak_cone_map, 'delimiter', '\t', 'newline', 'pc')


%% get the 2nd strongest cone feeding each cell of a given type

cell_type = 4;


% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', false,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% get second strongest cone
[num_cones, num_rgcs] = size(connectivity);
off_peak_cone_indices = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    temp_weights = connectivity(:,cc);
    [sorted_weights, sorted_indices] = sort(temp_weights, 'descend');
    off_peak_cone_indices(cc) = sorted_indices(2);
end

off_peak_cone_map = zeros(width, height);
verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(off_peak_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                new_counter = new_counter +1;
                off_peak_cone_map = off_peak_cone_map + temp_mask;
                if new_counter == 2
                    continue
                end
            end
        end
    end
end
off_peak_cone_map = off_peak_cone_map';
dlmwrite([cone_folder,'/map-0001.txt'], off_peak_cone_map, 'delimiter', '\t', 'newline', 'pc')


%% get the 1st and 2nd strongest cone feeding each cell of a given type
cell_type = 4;


% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% get strongest and second strongest cone
[num_cones, num_rgcs] = size(connectivity);
both_peak_cone_indices = zeros((2*num_rgcs),1);
for cc = 1:num_rgcs
    temp_weights = connectivity(:,cc);
    [sorted_weights, sorted_indices] = sort(temp_weights, 'descend');
    bg_index = 2*cc -1;
    en_index = 2*cc;
    both_peak_cone_indices(bg_index:en_index) = sorted_indices(1:2);
end

both_peak_cone_map = zeros(width, height);
verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(both_peak_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                new_counter = new_counter +1;
                both_peak_cone_map = both_peak_cone_map + temp_mask;
            end
        end
    end
end
both_peak_cone_map = both_peak_cone_map';
dlmwrite([cone_folder,'/map-0002.txt'], both_peak_cone_map, 'delimiter', '\t', 'newline', 'pc')



%% get the strongest cone feeding each cell of a given type (serial)

cell_type = 4;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued


[num_cones, num_rgcs] = size(connectivity);
peak_cone_indices = zeros(num_rgcs,1);
for cc = 1:num_rgcs
    [max_val, temp_cn_index] = max(connectivity(:,cc));
    peak_cone_indices(cc) = temp_cn_index;
end

peak_cone_map_serial = zeros(width, height);
verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(peak_cone_indices == cn,1))
        clear temp_mask
        if all(C{cn} ~=1)
            temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
            if length(find(temp_mask ==1)) <= 40
                if verbose
                    patch(V(C{cn},1),V(C{cn},2), 'r')
                    axis([1 320 1 320])
                    axis square
                    drawnow
                    pause(0.01)
                end
                new_counter = new_counter +1;
                peak_cone_map_serial = peak_cone_map_serial + (new_counter*temp_mask);
                if new_counter == 2
                    continue
                end
            end
        end
    end
end
peak_cone_map_serial = peak_cone_map_serial';
dlmwrite([cone_folder,'/first-serial.txt'], peak_cone_map_serial, 'delimiter', '\t', 'newline', 'pc')

%% serial cone stimulation from parasol cells
cone_folder = '2010-10-18-1-data003-on-parasol';
mkdir(cone_folder);
cell_list = [75 140 207];
number_cones = 20;
base_weight = 0.25;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_list,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

% get lists of cones that are > base_weight of peak
chosen_indices = zeros(length(cell_list), number_cones);
for cll = 1:length(cell_list)
    temp_indices = find(connectivity(:,cll) > base_weight);
    %cone_keeper{cll} = temp_indices
    temp_weights = connectivity(temp_indices,cll);
    [sorted_weights, sorted_indices] = sort(temp_weights, 'descend');
    temp_cone_num = length(temp_indices);
    step_size = temp_cone_num ./ number_cones;
    new_indices = round([1:step_size:temp_cone_num]);
    chosen_indices(cll,:) = temp_indices(sorted_indices(new_indices));
end


for nm = 1:number_cones
    serial_cone_map = zeros(width, height);
    verbose = true;
    figure(nm); clf; axis ij
    new_counter = 0;
    for cn = 1:length(C)
        if ~isempty(find(chosen_indices(:,nm) == cn,1))
            clear temp_mask
            if all(C{cn} ~=1)
                temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
                if length(find(temp_mask ==1)) <= 40
                    if verbose
                        patch(V(C{cn},1),V(C{cn},2), 'r')
                        axis([1 320 1 320])
                        axis square
                        drawnow
                        pause(0.01)
                    end
                    new_counter = new_counter +1;
                    serial_cone_map = serial_cone_map + temp_mask;
                    if new_counter == 2
                        continue
                    end
                end
            end
        end
    end
    serial_cone_map = serial_cone_map';
    dlmwrite([cone_folder,'/serial-',num2str(nm),'.txt'], serial_cone_map, 'delimiter', '\t', 'newline', 'pc')
end

%% Get max cone and list of cones with weights closest to that specified.
cell_type = {4};
%cone_folder = '2011-01-11-2-data003-cone-contrasts';
%mkdir(cone_folder);
weight_list = [1 0.75 0.5 0.25];

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_type,...
                                            'thresh', 0.2,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

[num_cones, num_rgcs] = size(connectivity);
% get lists of cones that are > base_weight of peak

cone_lists = zeros(num_rgcs, length(weight_list));
for wt = 1:length(weight_list)
    for cc = 1:num_rgcs
        temp_weights = connectivity(:,cc);
        max_weight = max(temp_weights);
        temp_weights = abs(temp_weights - (max_weight*weight_list(wt)));
        [junk, temp_cone_index] = min(temp_weights);
        cone_lists(cc, wt) = temp_cone_index;
    end
end


    
for wt = 1:length(weight_list)
    serial_cone_map = zeros(width, height);
    verbose = true;
    figure(wt); clf; axis ij
    new_counter = 0;
    for cn = 1:length(C)
        if ~isempty(find(cone_lists(:,wt) == cn,1))
            clear temp_mask
            if all(C{cn} ~=1)
                temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
                if length(find(temp_mask ==1)) <= 40
                    if verbose
                        patch(V(C{cn},1),V(C{cn},2), 'r')
                        axis([1 320 1 320])
                        axis square
                        drawnow
                        pause(0.01)
                    end
                    new_counter = new_counter +1;
                    serial_cone_map = serial_cone_map + temp_mask;
                    if new_counter == 2
                        continue
                    end
                end
            end
        end
    end
    serial_cone_map = serial_cone_map;
    dlmwrite(['map-000',num2str(wt-1),'.txt'], serial_cone_map, 'delimiter', '\t', 'newline', 'pc')
end




















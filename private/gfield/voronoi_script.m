
% define data path
datarun = load_data('/snle/acquisition/2010-09-24-1/data035/data035');
path_and_name{1,2} = 'orange-35';

% alternative data path
%datarun = load_data('/snle/lab/Experiments/Array/Analysis/2010-09-24-1/streamed/data006/data006');
%path_and_name{1,2} = '2010-09-24-1_streamed_data006_data006-bayes-msf_15.00-15-foo';

% load data
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


%% make cone map
verbose = false;

centers = datarun.cones.centers;
% an alternative tesselation
clear cone_map
[V,C] = voronoin(centers);
width = datarun.stimulus.field_width;
height = datarun.stimulus.field_height;
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
dlmwrite('~/Desktop/full_cone_map.txt', write_cone_map, 'delimiter', '\t', 'newline', 'pc')

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
dlmwrite('~/Desktop/off-midget-primary-cone.txt', peak_cone_map, 'delimiter', '\t', 'newline', 'pc')


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
dlmwrite('~/Desktop/off-midget-secondary.txt', off_peak_cone_map, 'delimiter', '\t', 'newline', 'pc')


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
dlmwrite('~/Desktop/off-midget-primary-secondary.txt', both_peak_cone_map, 'delimiter', '\t', 'newline', 'pc')



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
dlmwrite('~/Desktop/off-midget-primary-sequential.txt', peak_cone_map_serial, 'delimiter', '\t', 'newline', 'pc')

%% shared cones

cell_type = 4;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, {cell_type},...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   
connectivity = mosaic_weights .* selection; % keep weights continuous valued

[num_cones, num_rgcs] = size(connectivity);
cn_counter = 0;
clear shared_cone_index
for cn = 1:num_cones
    temp_cone_proj = connectivity(cn,:);
    if length(find(temp_cone_proj)) > 1
        cn_counter = cn_counter + 1;
        shared_cone_index(cn_counter) = cn;
    end
end
    

shared_cone_map_serial = zeros(width, height);
verbose = true;
figure(102); clf; axis ij
new_counter = 0;
for cn = 1:length(C)
    if ~isempty(find(shared_cone_index == cn,1))
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
                shared_cone_map_serial = shared_cone_map_serial + (new_counter*temp_mask);
            end
        end
    end
end
shared_cone_map_serial = shared_cone_map_serial';
dlmwrite('~/Desktop/on-midget-primary.txt', peak_cone_map, 'delimiter', '\t', 'newline', 'pc')











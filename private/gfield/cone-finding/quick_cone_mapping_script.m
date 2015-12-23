%data_path = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
%datarun = load_data(data_path);
%datarun = load_sta(datarun,'load_sta',[]);
%new_datarun = datarun;

% data paths
data_path{1} = '/snle/acquisition/2010-03-05-2/data012-0/data012-0';
data_path{2} = '/snle/acquisition/2010-03-05-2/data012-1/data012-1';
data_path{3} = '/snle/acquisition/2010-03-05-2/data012-2/data012-2';
data_path{4} = '/snle/acquisition/2010-03-05-2/data012-3/data012-3';

data_path{1} = '/snle/lab/Experiments/Array/Analysis/2010-03-05-2/data012/data012';

% load information from a data set
clear datarun
for dset = 1:length(data_path)
    datarun{dset} = load_data(data_path{dset});
    datarun{dset} = load_sta(datarun{dset},'load_sta',[]);
end



% collapse (join) data sets
clear new_datarun
new_datarun.cell_ids = [];
total_rgcs = 0;
for dset = 1:length(data_path)
    temp_cell_ids = datarun{dset}.cell_ids;
    new_datarun.cell_ids = [new_datarun.cell_ids, temp_cell_ids];
    total_rgcs = length(temp_cell_ids) + total_rgcs;
end
new_datarun.stas.stas = cell(1,total_rgcs);
new_datarun.stas.rfs = cell(1,total_rgcs);
new_datarun.stas.marks = cell(1,total_rgcs);
new_datarun.stas.rf_coms = cell(1,total_rgcs);
new_datarun.stas.time_courses = cell(1,total_rgcs);
new_datarun.stas.java_sta = datarun{1}.stas.java_sta;
new_datarun.stimulus = datarun{1}.stimulus;



clear spat_sens_params

% how to combine RGB values of the RF
spat_sens_params.strength = 'vector length';
%spat_sens_params.strength = {'inner or',...
%   [0.4044 0.8854 0.2292;0.1483 0.9161 0.3726;0.0174 0.0792 0.9967]};

% how to filter the RF before looking for significant stixels
spat_sens_params.filter =  struct('type','given','filt',...
           make_gaussian('center_radius',0.6,'x_size',5,'y_size',5,'center',[3 3]));

% how to find significant stixels
spat_sens_params.selection_params = struct('type','thresh','thresh',5);
%params.selection_params = struct('type','max','thresh',5);

% how to combine stixels between cells
spat_sens_params.combine_stixels = 'sum';
%cone_loc_params.combine_stixels = 'max';

% online readout of what's going on
spat_sens_params.verbose = true;
spat_sens_params.fig_single_cell = [];
spat_sens_params.foa_spat_sens = 105;


% accumulate the spatial sensitivity across all STAs
tic
spatial_sensitivity = zeros(datarun{1}.stimulus.field_width, datarun{1}.stimulus.field_height);
for dset = 1:length(data_path)
    [temp_spatial_sensitivity,all_sig_stixels,spatial_cell_ids] = compute_spatial_sensitivity(datarun{dset}, 'all', spat_sens_params);
    spatial_sensitivity = spatial_sensitivity + temp_spatial_sensitivity;
end
toc

figure(15)
imagesc(spatial_sensitivity)
print(15,'~/Desktop/sensitivity2.pdf', '-dpdf')


cone_centers = find_local_maxima(spatial_sensitivity);

% find indices to get cone locations
num_cones = length(find(cone_centers == 1));
[x_indices, y_indices] = find(cone_centers == 1);
centers = zeros(num_cones,2);
centers = [x_indices, y_indices];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONE WAY TO TESSELATE, BUT THE NEXT ONE IS RECOMMENDED OVER THIS ONE.
% to a voronoi tesselation on the cone centers
figure(2)
[V,C] = voronoin(centers);
for i = 1:length(C)
    if all(C{i} ~=1)
        patch(V(C{i},1), V(C{i},2),i)   
    end
end
axis([0,320,0,320])

% assign pixels to cones (make a pixelated voronoi tesselation)

% make pixel grid
field_width = new_datarun.stimulus.field_width;
field_height = new_datarun.stimulus.field_height;
pixel_grid = zeros(field_width*field_height,2);
tmp_counter_width = zeros(field_width*field_height,1);
tmp_counter_height = zeros(field_width*field_height,1);
for cntr = 1:field_height
    begin_index_width = 1+((cntr-1)*field_width);
    end_index_width = cntr*field_width;
    tmp_counter_width(begin_index_width:end_index_width) = cntr;
end
for cntr = 1:field_width
    begin_index_height = 1+((cntr-1)*field_height);
    end_index_height = cntr*field_height;
    tmp_counter_height(begin_index_height:end_index_height) = [1:1:field_height];
end
pixel_grid = [tmp_counter_width, tmp_counter_height];

% assign index to each pixel corresponding to closest cone
pixel_to_cone_indices = dsearchn(centers,pixel_grid);

% check distribution of tesselated areas
pixel_areas = zeros(1,num_cones); 
for cn = 1:num_cones
    pixel_areas(cn) = length(find(pixel_to_cone_indices == cn));
end
figure(10)
hist(pixel_areas,[0:2:110])
axis([0 100 0 300])
area_threshold = 30;
big_area_indices = find(pixel_areas > area_threshold);

for cn = big_area_indices
    temp_indices = find(pixel_to_cone_indices == cn);
    pixel_to_cone_indices(temp_indices) = 0;
end


% reshape to make a matrix the same dims as the display with index entries
cone_map = reshape(pixel_to_cone_indices, field_width,field_height);
 
% view image of this indexed matrix
imagesc(cone_map)

% write matrix to a file.
dlmwrite('~/Desktop/cone_map.txt', cone_map, '\t')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% an alternative tesselation
clear cone_map
[V,C] = voronoin(centers);
width = new_datarun.stimulus.field_width;
height = new_datarun.stimulus.field_height;
cone_map = zeros(width, height);
verbose = true;
figure(101); clf;
new_counter = 0;
for cn = 1:length(C)
    clear temp_mask
    if all(C{cn} ~=1)
        temp_mask = poly2mask(V(C{cn},1),V(C{cn},2), height, width);
        if length(find(temp_mask ==1)) <= 40
            if verbose
                patch(V(C{cn},2),V(C{cn},1), 'r')
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
cone_map = cone_map';
figure(2)
imagesc(cone_map)
dlmwrite('~/Desktop/cone_map.txt', cone_map, 'delimiter', '\t', 'newline', 'pc')

print(101, '~/Desktop/cone_map3.pdf','-dpdf')




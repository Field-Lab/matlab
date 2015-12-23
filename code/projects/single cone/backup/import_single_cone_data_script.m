

% load files
file_prefix = '/snle/lab/Experiments/Array/Shared/one/kiwi/';
rgcs_file = dlmread([file_prefix 'rgcs.txt'],'\t');
cones_file = dlmread([file_prefix 'cones.txt'],'\t');
connections_file = dlmread([file_prefix 'connections.txt'],'\t');

fprintf('rgcs:        [%d x %d] matrix.\n',size(rgcs_file,1),size(rgcs_file,2))
fprintf('cones:       [%d x %d] matrix.\n',size(cones_file,1),size(cones_file,2))
fprintf('connections: [%d x %d] matrix.\n',size(connections_file,1),size(connections_file,2))


% put file data into variables


% rgc info

% id number for each RGC
cell_ids = rgcs_file(:,1);

% which cell type (0 = none, 1 = on-parasol, 2 = off-parasol, 3 = on-midget, 4 = off-midget, 5 = sbc, 6+ = other)
cell_types = rgcs_file(:,2);

% whether a cell is in the region of interest (1 = in ROI, 0 = not in ROI)
cells_in_roi = rgcs_file(:,3);

% RGC center, based on center of mass of the cones
% first column is x coordinate, second column is y coordinate
rgc_COMs = rgcs_file(:,[4 5]);

% RGC center, based on gaussian fit
rgc_gaussian_fit_centers = rgcs_file(:,[6 7]);

% RGC radiu, based on gaussian fit
rgc_gaussian_fit_radii = rgcs_file(:,8);



% cone info

% id number for each cone
cone_ids = cones_file(:,1);

% center location of each cone
cone_centers = cones_file(:,[2 3]);

% color of each cone (L,M,S, or U for unknown)
cone_colors_temp = cones_file(:,4);
if 1
    cone_colors = char;
    cone_colors(cone_colors_temp == 1,1) = 'L';
    cone_colors(cone_colors_temp == 2,1) = 'M';
    cone_colors(cone_colors_temp == 3,1) = 'S';
    cone_colors(cone_colors_temp == 4,1) = 'U';
    clear cone_colors_temp
end

% measured red-green-blue color sensitivity for each cone
cone_rgb_values = cones_file(:,[5 6 7]);


    
% weights info

% C by R matrix (C = # of cones, R = # of ganglion cells), giving cone weights for each RGC
cone_weights = connections_file;



% clean up workspace
clear rgcs_file cones_file connections_file

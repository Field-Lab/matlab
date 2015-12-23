path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data003/data003/data003';
path_and_name{1,2} = '2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard';


% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);


% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [.5 .5 .5];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(10);clf;image(plot_mat);axis image; hold on


plot_cone_mosaic(datarun,'fig_or_axes',10,'cone_size',8,'clear',0,'bg_color',[])


%% remove S cones from data files.

new_datarun = datarun;

S_cone_indices = find(datarun.cones.types == 'S');

num_cones = length(datarun.cones.types);

cone_flags(1:num_cones) = true;

cone_flags(S_cone_indices) = false;

new_datarun.cones.centers = new_datarun.cones.centers(cone_flags,:);
new_datarun.cones.types = new_datarun.cones.types(cone_flags);
new_datarun.cones.rgb = new_datarun.cones.rgb(cone_flags,:);
new_datarun.cones.weights = new_datarun.cones.weights(cone_flags,:);

%% Get weights for new cones

load manuel_cone_marks

num_new_cones = length(manuel_cone_marks(:,1));

% parameters for new cone matrix W
rf_size = [320 320];
for cn = 1:num_new_cones
    kernel_spec(cn).radius = 0.65;
    kernel_spec(cn).center = manuel_cone_marks(cn,:);
    kernel_spec(cn).type = 'S';
end
kernel_colors.L = [20.8575   66.9626   16.4906];
kernel_colors.M = [ 8.1011   71.6135   25.1743];
kernel_colors.S = [1.0147   2.2233   18.2843];

[W,kernel_norms] = make_cone_weights_matrix(rf_size,kernel_spec,kernel_colors);


% get weights for new cones
% get list of cell indices
cell_indices = get_cell_indices(datarun,'all');

% initialize
cone_weights = zeros(num_new_cones,length(cell_indices));

% go through list of cells
for cc = 1:length(cell_indices)
    
    fprintf('.')

    % get summary frame
    rf = get_rf(datarun,datarun.cell_ids(cell_indices(cc)));

    if isempty(rf)
        continue
    end

    % reshape for the regression
    rf = reshape(rf,[],1);

    % put in units of SNR
    rf = rf / robust_std(rf);

    % regress to get the cone weights
    cone_weights(:,cc) = W\rf;

end

% display how long it took
%fprintf('\n   done (%0.1f seconds)\n',etime(clock,start_time_regress));
    

% clean up 
clear rf cc cell_indices


%% Add clicked cones

% add center points of new cones
new_datarun.cones.centers = [new_datarun.cones.centers; manuel_cone_marks];

% add type info for new cones
cone_labels = repmat('S', num_new_cones,1);
new_datarun.cones.types = [new_datarun.cones.types; cone_labels];

% add weights for new cones
new_datarun.cones.weights = [new_datarun.cones.weights; cone_weights];


% get size and color
y = datarun.stimulus.field_height;
x = datarun.stimulus.field_width;
rgb = [.5 .5 .5];

% generate background matrix
plot_mat = cat(3,repmat(rgb(1),y,x),repmat(rgb(2),y,x),repmat(rgb(3),y,x));

% plot it
figure(11);clf;image(plot_mat);axis image; hold on
plot_cone_mosaic(new_datarun,'fig_or_axes',11,'cone_size',8,'clear',0,'bg_color',[])





















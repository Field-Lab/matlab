% plantain
path_and_name{1,1} = '/snle/lab/Experiments/Array/Analysis/2008-08-27-5/data001/data001';
path_and_name{1,2} = 'plantain-1';
%path_and_name{1,2} = 'erroneous_normalization/2008-08-27-5_data003_data003_data003-bayes-msf_70.00--standard-old';

% load data
datarun = load_data(path_and_name{1,1});
datarun = load_params(datarun,struct('verbose',1));  
datarun = load_sta(datarun,'load_sta',[]);
datarun = set_polarities(datarun);
datarun = import_single_cone_data(datarun, path_and_name{1,2});    
datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);

% get indices for cell type

cell_type = {4};

cell_indices = get_cell_indices(datarun, cell_type);

% get cone indices for different types
l_indices = find(datarun.cones.types == 'L');
m_indices = find(datarun.cones.types == 'M');
s_indices = find(datarun.cones.types == 'S');

num_l_cones = length(l_indices);

% get cone connectivity to OFF midget cells
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_type,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);   

connectivity = mosaic_weights .* selection;


num_rgcs = length(cell_indices);


cone_type = 'M';
cone_indices = find(datarun.cones.types == cone_type);

cone_clr_tc = zeros(num_frames,num_clrs);
for rgc = 1:num_rgcs

    rgc % counter

    % get sta for rgc of interest
    temp_sta = get_sta(datarun, datarun.cell_ids(cell_indices(rgc)));
    num_frames = size(temp_sta,4);
    num_clrs = size(temp_sta,3);

    % get indices to cone with sig weights
    temp_cone_indices = find(connectivity(:,rgc));
    % find those that are L cones
    temp_cone_index = intersect(cone_indices,temp_cone_indices);

    if ~isempty(temp_cone_index)

        % get rgb triplets for these l cones
        for cn = 1:length(temp_cone_index)

            %get cones weighting function
            temp_wc = Wc(:,temp_cone_index(cn));
            temp_wc = reshape(full(temp_wc), [320, 320, 3]);
            temp_max = max(max(temp_wc(:,:,2)));
            temp_wc = temp_wc(:,:,2) ./ temp_max;
            temp_wc = reshape(temp_wc,[],1);

            project_cone_sta = zeros(num_frames, num_clrs);
            for clr = 1:size(temp_sta,3);
                temp_clr = reshape(squeeze(temp_sta(:,:,clr,:)),[],num_frames);
                temp_clr = temp_wc' *temp_clr;
                project_cone_sta(:,clr) = temp_clr;
            end

            cone_clr_tc = cone_clr_tc + project_cone_sta;
        end
    end
end

 
base_step = 1/120;
step_end = base_step*(num_frames-1);
time_stamps = 0:base_step:step_end;
time_stamps = time_stamps * -1;
time_staps = sort(time_stamps, 'ascend');
figure
hold on
plot(time_staps, cone_clr_tc(:,1), 'r')
plot(time_staps, cone_clr_tc(:,2), 'g')
plot(time_staps, cone_clr_tc(:,3), 'b')
hold off
print(5, '~/Desktop/s_sbc_cone.pdf', '-dpdf')
   


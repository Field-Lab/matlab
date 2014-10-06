% LNP model of RGC

% 
% keep_datarun = datarun;
% 
% % load data structure
% cd /Analysis/gfield/
% load datarun
% 
% old_datarun = datarun;
% datarun = keep_datarun;
% clear keep_datarun
% 
% datarun.stas.snls{147} = old_datarun.stas.snls{147};
% 
% % load cone weights matrix
% load([single_cone_path datarun.names.nickname '/Wc.mat'])
% 
% 

cone_weight_params.thresh = 0.1;
cone_weight_params.radius = [0 inf];
cone_weight_params.polarity = 1;
cone_weight_params.contiguity = true;
cone_weight_params.scale = 3.0;

dataname = 'peach';

%% load data
[datarun_junk, cone_info] = load_data_and_cones(dataname, 'sta_summaries', false);

%% load datarun from server
cd([single_cone_path, 'saved/']);
load(dataname)


rgc_tcs = cone_info.rgc_tcs;
spike_times = cone_info.spike_times;
cone_inputs = cone_info.cone_inputs;
cone_ids = cone_info.cone_ids;
rgc_ids = cone_info.rgc_ids;
spike_rate = cone_info.spike_rate;
Wc = cone_info.Wc;
height = datarun.stimulus.field_height;
width = datarun.stimulus.field_width;
datarun.cones = datarun_junk.cones;
datarun = set_polarities(datarun);
clear cone_info datarun_junk



% get single cone data
%datarun = import_single_cone_data(datarun, '2008-08-27-0_data001-0s-2400s_data001_data001-bayes-msf_20.00--standard');
%datarun.cones.mosaic = make_mosaic_struct(datarun.cones.centers);
%datarun = load_sta(datarun,'load_sta',[]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% STIMULUS INFORMATION %%%
% refresh time
% refresh_time = datarun.stimulus.java_movie.getRefreshTime/1000;
% 
% % compute start and stop times in terms of stimulus number
% start_stim = floor(1+start_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
% end_stim = floor(1+end_time/(datarun.stimulus.java_movie.getRefreshTime/1000));
% stims = start_stim:end_stim;
% 
% % note stimulus size
% field_width = datarun.stimulus.java_movie.getWidth;
% field_height = datarun.stimulus.java_movie.getHeight;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define cell of interest
cell_id = 3334;
grid_points = [165 140; 163 150; 170 138; 171 144; 169 153; 178 135; 185 148; 195 140; 193 147;...
    180 142; 181 151; 175 156; 188 135; 190 142; 193 155; 187 160; 200 145; 201 150; 205 145; 175 149]; 
grid_points = [grid_points(:,2) grid_points(:,1)];
    
temp_cell_index = get_cell_indices(datarun, cell_id);
template_length = 6;

% extract connectivity
[mosaic_weights, selection, extras] = select_cone_weights(datarun, cell_id,...
                                            'thresh', 0.1,...
                                            'radius', [0 inf], 'polarity', 1,...
                                            'contiguity', true,'scale', 3.0);
                                        
% get image frame info
cell_com = datarun.cones.rf_fits{temp_cell_index}.center;
window = 35;
frame_begin_x = floor(cell_com(2) - window);
frame_end_x = floor(cell_com(2) + window);
frame_begin_y = floor(cell_com(1) - window);
frame_end_y = floor(cell_com(1) + window);
matrix_scale_factor = 10;


                                        
temp_cone_indices = find(selection);

% get time course from cell of interest
temp_tc = datarun.stas.time_courses{temp_cell_index};
num_frames = size(temp_tc,1);

if size(temp_tc,2) == 3
    temp_tc = temp_tc(num_frames - template_length+1:num_frames,:);
    temp_tc = sum(temp_tc, 2);
elseif size(temp_tc,2) ==1
    % cut out region of interest from TC
    temp_tc = temp_tc(num_frames - template_length+1:num_frames);
else
    error('number of color channels does not equal either 1 or 3')
end

temp_tc = temp_tc - mean(temp_tc);
temp_tc = temp_tc ./ norm(temp_tc);

% flip the impulse response to estimate the temporal filter
reverse_indices = template_length:-1:1;
impulse_filter = temp_tc(reverse_indices);
%impulse_filter = temp_tc;

cone_stim = filter(impulse_filter, 1, cone_inputs(:,temp_cone_indices));


% collect the RFs of the cones feeding the RGC of interest
cone_rfs = Wc(:, temp_cone_indices);
for cn = 1:size(cone_rfs,2);
    cone_rfs(:,cn) = cone_rfs(:,cn) ./ norm(cone_rfs(:,cn));
end

% LINEARLY SUM the CONE INPUTS;
% get cone weights to "gaussian" weight summation
% normalize so that weights sum to unity
weights = mosaic_weights(temp_cone_indices)' ./ sum(mosaic_weights(temp_cone_indices));

%fit_params = datarun.stas.snls{temp_cell_index}.fit_params;

nl_model = 'subunit';
overlapping = false;
subunit_NL_form = 'heaviside';
num_clusters = 5; % only used in 'subunit' model

switch nl_model
    case 'cone_nl'
        cone_generator = exp((fit_params.b) + (fit_params.a) * cone_stim);
        cone_generator = cone_generator .* repmat(weights,size(cone_stim,1),1);
    case 'end_nl'
        cone_generator = cone_stim .* repmat(weights,size(cone_stim,1),1);
    case 'cone_saturation'
        cone_generator = atan(cone_stim);
    case 'subunit'
        locations = datarun.cones.centers(temp_cone_indices,:);
        
        if overlapping
            grid_convergence = 12;  % number of cones converging to a grid point
            num_grid_points = size(grid_points, 1);
            cones_to_grid = zeros(length(locations), num_grid_points);
            for gp = 1:num_grid_points
                % calculate distance between grid point and cone locations
                grid_dists = ipdm(grid_points(gp,:), locations);
                % sort the distances, ascending
                [sorted_grid_dists, grid_dist_indices] = sort(grid_dists, 'ascend');
                % get X number of closest cones
                cones_to_grid(grid_dist_indices(1:grid_convergence),gp) = 1;
            end
            
            % make sure each cone connects to at least one grid point
            for cn = 1:size(locations,1)
                if isempty(find(cones_to_grid(cn,:) > 0, 1))
                    grid_dists = ipdm(locations(cn,:), grid_points);
                    [min_val, min_index] = min(grid_dists);
                    cones_to_grid(cn,min_index) = 1;
                    fprintf('cone %0.1f connected to grid point %0.1f \n',cn, min_index);
                end
            end
            num_clusters = num_grid_points;
            [U,S,V] = svd(cones_to_grid);
            figure
            plot(diag(S))
            
        else
            % if overlapping = false, make orthogonal, non-overlapping, groups of cones
            IDX = kmeans(locations, num_clusters);
            % plot the clusters as a check on clustering
            if 1
                figure(10); clf; 
                color_pal = {[1 1 1],[1 0 1], [1 1 0], [0 1 1], [1 0 0], [0 1 0], [0 0 1], [0.5 1 0], [1 0.5 0], [1 0.5 0.5], [0 1 0.5], [0 0.5 1]};
                rf_image = zeros((320*320*3),1);
                for clust = 1:num_clusters
                    clust_indices = find(IDX == clust);
                    sub_cone_rfs = cone_rfs(:,clust_indices);
                    for rf = 1:size(sub_cone_rfs,2)
                        clr = color_pal{mod(clust,length(color_pal))+1};
                        tmp_rf = reshape(full(sub_cone_rfs(:,rf)),[320,320,3]);
                        [tmp_rw, tmp_cl] = find(squeeze(tmp_rf(:,:,2)));
                        px_weights = tmp_rf(tmp_rw, tmp_cl, 2);
                        tmp_rf(tmp_rw,tmp_cl,1) = clr(1)*px_weights;
                        tmp_rf(tmp_rw,tmp_cl,2) = clr(2)*px_weights;
                        tmp_rf(tmp_rw,tmp_cl,3) = clr(3)*px_weights;
                        sub_cone_rfs(:,rf) = reshape(tmp_rf, [],1);
                    end
                    rf_image = rf_image + sum(sub_cone_rfs,2);
                end
                nw_rf = reshape(rf_image, [320 320 3]);
                imagesc(matrix_scaled_up(nw_rf,matrix_scale_factor))
                axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
                axis off; axis square
                
            end
        end
            
        % apply NL to output of subunits
        sub_gen = zeros(size(cone_stim,1),num_clusters);
        
        for clust = 1:num_clusters
            
            if overlapping
                clust_indices = find(cones_to_grid(:,clust) == 1);
            else
                clust_indices = find(IDX == clust);
            end
            
            sub_gen_ = sum(cone_stim(:,clust_indices), 2);
            sub_gen_ = sub_gen_ ./ (abs(ext(sub_gen_)));
            %figure(20)
            %hist(sub_gen_,100)    

            switch subunit_NL_form
                
                case 'exp'
                    % exponential NL
                    %sub_gen(:,clust) = exp(-fit_params.b + fit_params.a * sub_gen_)./sqrt(num_clusters);
                    sub_gen(:,clust) = exp(-2 + 4 * sub_gen_);

                case 'heaviside'
                    % heaviside NL
                    squashed_indices = find(sub_gen_ < 0.4);
                    amped_indices = find(sub_gen_ >= 0.4);
                    sub_gen(squashed_indices,clust) = 0;
                    sub_gen(amped_indices, clust) = 5;

                case 'threshold'
                    % thresholing NL
                    low_indices = find(sub_gen_ < 0.5);
                    sub_gen_(low_indices) = 0;
                    sub_gen(:,clust) = sub_gen_*2;
                
                case 'rectified_linear'
                    thresh = 0;
                    y_int = 0;
                    slp = 5;
                    sub_gen_ = (slp*sub_gen_) + y_int;
                    sub_gen_(sub_gen_ < 0) = 0;
                    sub_gen(:,clust) = sub_gen_;
                    
                    hist(sub_gen_)
                    pause
                    
                case 'squaring'
                    sub_gen(:,clust) = sub_gen_.^2;
                    
            end
            
           %figure(21)
           %hist(sub_gen(:,clust),100)
            
        end
        cone_generator = sub_gen;

    otherwise
        error('start time and end time not set')
end



if strcmp(nl_model,'end_nl')
    net_generator_sig = sum(cone_generator, 2);
    % rescale so that the range of generator signals spans -0.5 to 0.5 
    desired_range = 1.0;
    gen_range = range(net_generator_sig);
    net_generator_sig = net_generator_sig .* (desired_range/gen_range);
    NL_generator_sig = exp(fit_params.b + fit_params.a * net_generator_sig);
    poiss_train = poissrnd(NL_generator_sig);
else 
    net_cone_input = sum(cone_generator, 2);
    % rescale the net cone input so that the generator signal produces some
    % multiple of the firing rate of the neuron
    expected_spike_number = mean(net_cone_input) * length(net_cone_input);
    desired_spike_number_factor = 5; % relative to the firing rate of the neuron used for the model;
    data_spike_number = length(datarun.spikes{temp_cell_index});
    input_rescale_factor = data_spike_number * desired_spike_number_factor ./expected_spike_number ;
    generator_sig = net_cone_input * input_rescale_factor;
    % produce spike train from linear->nonlinear generator signal
    poiss_train = poissrnd(generator_sig);
end

sum(poiss_train) % report the number of spikes


% IDENTIFY SPIKES IN STIMULUS BINS
% translate to spike times (with duplicates for multiple spikes per time bin)
spike_times_ = [];
for nn = 1:max(poiss_train)
    spike_times_ = [spike_times_; find( poiss_train > (nn-1) )];
end
model_cell_spike_times = sort(spike_times_);


%%
% CONSTRUCT THE SPIKE-TRIGGERED ENSEMBLE
STE = cone_stim(model_cell_spike_times,:);

% get image frame info
cell_com = datarun.cones.rf_fits{temp_cell_index}.center;
window = 35;
frame_begin_x = cell_com(2) - window;
frame_end_x = cell_com(2) + window;
frame_begin_y = cell_com(1) - window;
frame_end_y = cell_com(1) + window;

% --- COMPUTE STA ---
temp_polarity = datarun.stas.polarities{temp_cell_index};
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;

% --- PLOT STA ---
figure(102)
image_STA = Wc(:,temp_cone_indices) * STA';
imagesc(matrix_scaled_up(norm_image(reshape(image_STA, [320 320 3])), matrix_scale_factor))
axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
axis off; axis square
%print(102, '~/Desktop/sta.pdf', '-dpdf')

% --- PLOT STA with subunits ---
if overlapping
    figure(101); clf;
    image_STA = Wc(:,temp_cone_indices) * STA';
    imagesc(matrix_scaled_up(norm_image(reshape(image_STA, [320 320 3])), matrix_scale_factor))
    axis(matrix_scale_factor * [frame_begin_y frame_end_y frame_begin_x frame_end_x])
    hold on
    axis off; axis square
    % plot grid points
    plot(grid_points(:,1), grid_points(:,2), 'k.');
    % plot the lines between grid points and cones
    for gd = 1:num_grid_points
        temp_indices = find(cones_to_grid(:,gd) == 1);
        for cn = 1:length(temp_indices)
            plot(matrix_scale_factor*[locations(temp_indices(cn),1),grid_points(gd,1)],matrix_scale_factor*[locations(temp_indices(cn),2),grid_points(gd,2)], 'k')
        end
    end
end

% --- COMPUTE STC ---
STA_ = mean(STE)';
z_STE = STE';
z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
%STE = STE - (STE*STA_') * STA_ ./ norm(STA_)^2;  % project out the mean
z_filtered_cone_inputs = cone_stim';
z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean

STC = cov(z_STE') - cov(z_filtered_cone_inputs');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% get indices to meaningful dimensions
keep_indices = find(abs(eig_vals) > 1e-10); % note this number might need to be changed depending on data
% store eigvals;
eig_vals = eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
PCs = PCs(:,keep_indices);

% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = num_clusters-1;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,dm);
    figure(dm)
   image(norm_image(reshape(image_PC, [320,320,3])))
    axis([frame_begin_y frame_end_y frame_begin_x frame_end_x])

%   image(matrix_scaled_up(norm_image(reshape(image_PC, [320,320,3])), matrix_scale_factor))
%   axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off; axis square
    %save_path = ['~/Desktop/PC',num2str(dm),'.pdf'];
    %print(dm, save_path, '-dpdf', '-painters')
end
% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(num_dims + 1); clf;
plot(eig_vals, 'ko')

%%
num_sig_dims = 4;

num_iters = 100;
time_offsets = ceil(10000*rand(num_iters,1));
shifted_eig_vals = zeros(num_iters, length(temp_cone_indices) - num_sig_dims - 1);
for iter = 1:num_iters
    
    % bin up spikes for entire duration
    s_spike_times = 1+mod(model_cell_spike_times+time_offsets(iter), size(cone_stim,1));
    s_STE = cone_stim(s_spike_times,:); 
    
    
    s_STA_ = mean(s_STE)';
    z_STE = s_STE';
    z_STE = z_STE - STA_*(STA_'*z_STE)./ norm(STA_)^2;  % project out the mean
    z_filtered_cone_inputs = cone_stim';
    z_filtered_cone_inputs = z_filtered_cone_inputs - STA_*(STA_'*z_filtered_cone_inputs)./ norm(STA_)^2;  % project out the mean
    for dm = 1:num_sig_dims
        z_STE = z_STE - PCs(:,dm)*(PCs(:,dm)'*z_STE)./ norm(PCs(:,dm))^2;
        z_filtered_cone_inputs = z_filtered_cone_inputs - PCs(:,dm)*(PCs(:,dm)'*z_filtered_cone_inputs)./ norm(PCs(:,dm))^2;
    end
    s_STC = cov(z_STE') - cov(z_filtered_cone_inputs');
    
    % --- FACTORIZE ---
    [s_PCs, s_eig_val_matrix] = eig(s_STC);
    s_eig_vals = diag(s_eig_val_matrix);
    [s_eig_vals, s_sorted_eig_indices] = sort(s_eig_vals, 'descend');
    s_PCs = s_PCs(:, s_sorted_eig_indices);  

    % get indices to meaningful dimensions
    keep_indices = find(abs(s_eig_vals) > 1e-10); % note this number might need to be changed depending on data
    % store eigvals;
    shifted_eig_vals(iter,:) = s_eig_vals(keep_indices); % ./ abs(sum(s_eig_vals));
end

num_sds = 3;
mean_shifted_egvals = mean(shifted_eig_vals);
sd_shifted_egvals = std(shifted_eig_vals);
high_shifted_egvals = mean_shifted_egvals + num_sds*sd_shifted_egvals;
low_shifted_egvals = mean_shifted_egvals - num_sds*sd_shifted_egvals;

figure(103); clf;
max_eig_val = max(eig_vals);
plot(eig_vals./max_eig_val, 'k.', 'MarkerSize', 20)
hold on
temp_counter = (num_sig_dims+1):1:length(eig_vals);
plot(temp_counter,high_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
plot(temp_counter,low_shifted_egvals./max_eig_val, 'r-', 'LineWidth', 2)
axis([0, size(PCs,1)+1, -0.5 1.1])
title(['PC spectrum ', num2str(num_sig_dims)])
hold off
%print(103, '~/Desktop/spectrum.pdf', '-dpdf')

% --- PLOT SHIFTED EIGENVECTORS OF THE STC ---
for dm = 1:2
    image_shifted_PCs = Wc(:,temp_cone_indices) * s_PCs(:,dm);
    figure(dm+50)
    image(matrix_scaled_up(norm_image(reshape(image_shifted_PCs, [320,320,3])),matrix_scale_factor))
    axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off
    axis square
    save_path = ['~/Desktop/PC_s',num2str(dm),'.pdf'];
    %print(dm+50, save_path, '-dpdf')
end

%%



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% try nnmf
num_factors = 10;
[W, H] = nnmf(spike_frames'+1.75,num_factors);

%W_orth = orth(W);

for comp = 1:num_factors;
    figure(comp)
    cone_egvec = cone_rfs * W(:,comp);
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    title_text = ['NNF ', num2str(comp)];
title(title_text)
end

% try projecting the spike frames into the low-D space determined by PCA + STA and compute the subunits
low_d_spike_frames = spike_frames * PCs_plus(:,1:num_clusters);
[W, H] = nnmf(low_d_spike_frames'+3,3);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try k-means clustering to look for vectors that aptly cluster the data
num_VQs = 8;

% make response subspace from sig PCs and STA
response_subspace = [PCs(:,1:num_VQs-1), (STA'./norm(STA))];
% project stimulus into this subspace
stim_weights = STE * response_subspace;
% cluster the stimulus in this subspace

% white data along each dim
for dm = 1:size(stim_weights,2)
    stim_weights(:,dm) = stim_weights(:,dm) ./ std(stim_weights(:,dm));
end

[IDX, VQs] = kmeans(stim_weights, num_VQs, 'start', 'uniform', 'MaxIter', 200, );


% plot(stim_weights(:,1), stim_weights(:,2), 'k.')
% hold on
% plot(VQs(:,1), VQs(:,2), 'r.')
% hold off
% 
% scatter3(stim_weights(:,1),stim_weights(:,2), stim_weights(:,3), 'k.')
% hold on
% scatter3(VQs(:,1),VQs(:,2),VQs(:,3), 'r.')
% hold off
% 
% 
% scatter3(stim_weights((IDX == 1),1),stim_weights((IDX == 1),2), stim_weights((IDX == 1),3), 'b.')
% hold on
% scatter3(stim_weights((IDX == 2),1),stim_weights((IDX == 2),2), stim_weights((IDX == 2),3), 'g.')
% scatter3(stim_weights((IDX == 3),1),stim_weights((IDX == 3),2), stim_weights((IDX == 3),3), 'r.')
% scatter3(VQs(:,1),VQs(:,2),VQs(:,3), 'k.')
% hold off

cone_space_kernels = response_subspace * VQs';

matrix_scale_factor = 10;

for comp = 1:num_VQs;
    figure(comp)
    cone_egvec = cone_rfs * cone_space_kernels(:,comp);
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(matrix_scaled_up(norm_image(cone_egvec), matrix_scale_factor))
    axis off; axis image
    axis([frame_begin_y frame_end_y frame_begin_x frame_end_x] * matrix_scale_factor)
    title_text = ['k ', num2str(comp)];
    title(title_text)
    save_path = ['~/Desktop/cluster',num2str(comp),'.pdf'];
    %print(comp, save_path, '-dpdf', '-painters', '-r600')
    
end




%%
% set starting weights
start_weights = diag(ones(5,1));


% try to minimize
[temp, fval] = fminsearch(@(init_weights) sum_norms(response_subspace, init_weights), start_weights);

fval

new_space = response_subspace * temp;

for comp = 1:num_VQs;
    figure(comp)
    cone_egvec = cone_rfs * new_space(:,comp);
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    title_text = ['k ', num2str(comp)];
title(title_text)
end

%% try st-ICA

num_sig_dims = 5;
whiteMat = (diag(eig_vals(1:num_sig_dims)));
new_STE = (diag(eig_vals(1:num_sig_dims))).^0.5 * PCs(:,1:num_sig_dims)' * z_STE;

    
ICs = fastica(z_STE, 'g', 'gauss', 'whiteSig', new_STE, 'whiteMat', whiteMat, 'dewhiteMat', inv(whiteMat));


ICs = fastica(z_STE, 'g', 'gauss', 'lastEig', 5);


[ICs, A, W] = fastica(new_STE, 'g', 'gauss', 'pcaE', PCs(:,1:num_sig_dims)', 'pcaD', whiteMat);



[A, W] = fastica(new_STE, 'g', 'gauss');

temp_sub_space = [STA_, PCs(:,1:num_sig_dims)];

[ICs, A, W] = fastica(STE, 'g', 'gauss', 'pcaE', temp_sub_space);

num_dims = 5;
for dm = 1:num_dims
    image_PC = Wc(:,temp_cone_indices) * PCs(:,1:num_sig_dims) * W(dm,:)';
    image_PC = reshape(image_PC, [height, width, 3]);
    image_PC = image_PC(frame_begin_x:frame_end_x,frame_begin_y:frame_end_y,:);
    figure(dm+10)    
    image(norm_image(image_PC))
    axis off; axis square
    title(['D',num2str(dm)])
    drawnow
    norm(sum(norm_image(image_PC),3),1)
end


%%
%%%  Histogram the weights
% 
% for dm = 1:4
%     figure(10+dm)
%     hist(PCs(:,dm), 60)
% end

IDX = kmeans(PCs(:,1:(num_clusters-1)),num_clusters);


% --- PLOT EIGENVECTORS OF THE STC with increased variance---
num_dims = num_clusters;
for dm = 1:num_dims

    temp_IDs = find(IDX == dm);
    ID_indices = zeros(1,134);
    ID_indices(temp_IDs) = 1;

    
    subunit_image = Wc(:,temp_cone_indices) * ID_indices';
    figure(dm+20)
   %image(norm_image(reshape(subunit_image, [320,320,3])))
    %axis([frame_begin_y frame_end_y frame_begin_x frame_end_x])

   image(matrix_scaled_up(norm_image(-1*reshape(subunit_image, [320,320,3])), matrix_scale_factor))
   axis(matrix_scale_factor*[frame_begin_y frame_end_y frame_begin_x frame_end_x])
    axis off; axis square
    save_path = ['~/Desktop/PC',num2str(dm),'.pdf']
    print(dm+20, save_path, '-dpdf', '-painters')
end







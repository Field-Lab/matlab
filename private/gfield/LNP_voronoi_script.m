% LNP model of RGC

weights = STA;

cone_stim = filtered_patches;


nl_model = 'end_nl';
overlapping = false;
subunit_NL_form = 'threshold';
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
            
        else
            % if overlapping = false, make orthogonal, non-overlapping, groups of cones
            IDX = kmeans(locations, num_clusters);
            % plot the clusters as a check on clustering
            if 1
                figure(10); clf; 
                color_pal = {[1 0 1], [1 1 0], [0 1 1], [1 1 1], [1 0 0], [0 1 0], [0 0 1], [0.5 1 0], [1 0.5 0], [1 0.5 0.5], [0 1 0.5], [0 0.5 1]};
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
                imagesc(nw_rf)
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
            sub_gen_ = sub_gen_ ./ (abs(ext(sub_gen_)) *2);
            %hist(sub_gen_)    

            switch subunit_NL_form
                
                case 'exp'
                    % exponential NL
                    sub_gen(:,clust) = exp(-fit_params.b + fit_params.a * sub_gen_)./sqrt(num_clusters);
                    sub_gen(:,clust) = exp(-1 + 3 * sub_gen_);
            
                case 'heaviside'
                    % heaviside NL
                    squashed_indices = find(sub_gen_ < 0.1);
                    amped_indices = find(sub_gen_ >= 0.1);
                    sub_gen(squashed_indices,clust) = 0;
                    sub_gen(amped_indices, clust) = 1;

                case 'threshold'
                    % thresholing NL
                    low_indices = find(sub_gen_ < 0.1);
                    sub_gen_(low_indices) = 0;
                    sub_gen(:,clust) = sub_gen_*1;
            end

           % hist(sub_gen(:,clust),100)
            
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
    NL_generator_sig = exp(-1 + 3 * net_generator_sig);
    poiss_train = poissrnd(NL_generator_sig);
else 
    net_cone_input = sum(cone_generator, 2);
    generator_sig = net_cone_input;
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

% get the stimulus frames that produced spikes
spike_frames = cone_stim(model_cell_spike_times,:);

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
temp_polarity = -1;
STA = mean(STE);
STA = STA ./ max(STA) ./2;
STA = STA * temp_polarity;


% --- PLOT STA ---
figure(101)
for mark = 1:num_marks
    map_indices = find(cone_map == mark_indices(mark));
    temp_cone_map(map_indices) = STA(mark);
end
image(norm_image(repmat(temp_cone_map, [1 1 3])))



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


% --- PLOT EIGENVECTORS OF THE STC ---
num_dims = 4;
for dm = 1:num_dims
    temp_map = zeros(size(cone_map));
    temp_pc = PCs(:,dm); % ./ (2 * abs(ext(PCs(:,dm)))) + 0.5; 
    for mark = 1:num_marks
        map_indices = find(cone_map == mark_indices(mark));
        temp_map(map_indices) = temp_pc(mark);
    end
    figure(dm)
    image(norm_image(repmat(temp_map, [1,1,3])))
    axis off
end

% --- PLOT THE SORTED EIGENVALUE SPECTRUM ---
figure(num_dims + 1); clf;
plot(eig_vals, 'ko')

%%
num_sig_dims = 0;

num_iters = 100;
time_offsets = ceil(10000*rand(num_iters,1));
shifted_eig_vals = zeros(num_iters, length(mark_indices) - num_sig_dims - 1);
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
temp_counter = (num_sig_dims+2):1:length(mark_indices);
plot(temp_counter,high_shifted_egvals, 'r-')
hold on
plot(temp_counter,low_shifted_egvals, 'r-')
plot(eig_vals, 'k.')
%plot(eig_vals./abs(sum(eig_vals)), 'ko')
%axis([0 100 -0.003 0.01])
hold off

% --- PLOT SHIFTED EIGENVECTORS OF THE STC ---
for dm = 1:num_dims
    image_shifted_PCs = Wc(:,temp_cone_indices) * s_PCs(:,dm);
    figure(dm+50)
    image(norm_image(reshape(image_shifted_PCs, [320,320,3])))
    axis([frame_begin_y frame_end_y frame_begin_x frame_end_x])
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% try nnmf
num_factors = 5;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try k-means clustering to look for vectors that aptly cluster the data
num_VQs = 10;
[IDX, VQs] = kmeans(spike_frames, num_VQs);

%[IDX, VQs] = kmeans(shifted_spike_frames, num_VQs);

for comp = 1:num_VQs;
    figure(comp)
    cone_egvec = cone_rfs * VQs(comp,:)';
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    title_text = ['k ', num2str(comp)];
title(title_text)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate noise eigenvalues by creating shifting spike times relative to stimulus
num_iters = 50;
time_offsets = -6-round(10000*rand(num_iters,1));
shifted_egvals = zeros(num_iters, length(cone_ids));
for iter = 1:num_iters

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    shifted_spike_rate_ = circshift(poiss_train(start_stim:end_stim),time_offsets(iter));

    % translate to spike times (with duplicates for multiple spikes per time bin)
    shifted_spike_times = [];
    for nn = 1:max(shifted_spike_rate_)
        shifted_spike_times = [shifted_spike_times; find( shifted_spike_rate_ > (nn-1) )];
    end
    shifted_spike_times = sort(shifted_spike_times);

    shifted_spike_frames = cone_stim(shifted_spike_times,:);

    [shifted_PCs, shifted_weights, shifted_egvals_] = princomp(shifted_spike_frames);

    shifted_egvals(iter,:) = shifted_egvals_'./sum(shifted_egvals_);

end

ave_egval = mean(shifted_egvals);
std_egval = std(shifted_egvals);

figure(20)
clf; hold on; 
plot(eigvals./sum(eigvals), 'ro')
plot(ave_egval+(3*std_egval), 'k')
plot(ave_egval-(3*std_egval), 'k')
xlabel('PCs')
ylabel('relative variance')
title('red = data, black = +/- 3-sigma noise')



% get clumping estimates on the eigenvectors
[shifted_PC_indices, temp_CIs] = compute_clumping_on_spatial_PC(datarun, cone_ids, shifted_PCs);

figure(4)
cone_egvec = cone_rfs * shifted_PCs(:,shifted_PC_indices(1));
cone_egvec = reshape(cone_egvec, [320, 320, 3]);
imagesc(norm_image(cone_egvec))
title('PI w/ high clumping')

num_PCs = 8;
for comp = 1:num_PCs;
    figure(comp)
    cone_egvec = cone_rfs * shifted_PCs(:,comp);
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    title_text = ['PC ', num2str(comp)];
title(title_text)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%


ICAs = fastica(spike_frames', 'numOfIC', 10);



















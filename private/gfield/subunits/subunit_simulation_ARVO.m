% THIS SCRIPT SETS UP A SIMULATION BASED ON A REAL OFF PARASOL GANGLION CELL
% FROM THE PRIMATE RETINA.  

cd /Analysis/gfield/for-pillow/

% ---------------------------------------------------------------
% --- LOAD VARIOUS INFORMATION FROM DISK ---

% LOAD THE CONE WEIGHTS -- MANY ELEMENTS ARE ZERO.  THESE CORRESPOND TO CONES 
% THAT DID NOT CONTRIBUTE TO THE RF OF THIS GAGNLION CELL
load cone_weights
cone_indices = find(cone_weights);
% normaize the such that they sum to 1
cone_weights = cone_weights(cone_indices)' ./ sum(cone_weights(cone_indices));


% LOAD THE CENTER LOCATIONS OF ALL THE CONES
load cone_centers


% LOAD THE RECEPTIVE FIELDS OF ALL THE CONES -- INCLUDING THOSE WHICH HAVE
% A ZERO WEIGHT ON THE GANGLION CELL RF
load all_cone_rfs
% collect the RFs of the cones feeding the RGC of interest
cone_rfs = all_cone_rfs(:, cone_indices);
% normalize each cone RF
for cn = 1:size(cone_rfs,2);
    cone_rfs(:,cn) = cone_rfs(:,cn) ./ norm(cone_rfs(:,cn));
end


% LOAD THE STIMULUS: THE ORITINAL BINARY RGB WHITE NOISE STIMULUS HAS BEEN 
% PROJECTED THROUGH THE CONES.  NOW THE STIMULUS CAN BE THOUGHT OF AS A
% GENERATOR SIGNAL FROM EACH CONE.
% store stimulus as "cone_inputs"
% COMMENTED OUT TO MAKE LOADED INFO SMALLER
cone_fid = fopen('cone_input.bin','r', 'b');
cone_inputs = fread(cone_fid, [48015, 2100], 'float32');
fclose(cone_fid);

% get cone stim for the RF of the cell of interest
% COMMENTED OUT TO MAKE LOADED INFO SMALLER
%rf_cone_inputs = cone_inputs(:, cone_indices);

% skip loading the stimuli for all cones, and only load for the cones
% with non-zero weights to the cell of interest
cone_fid = fopen('rf_cone_input.bin','r', 'b');
rf_cone_inputs = fread(cone_fid, [48015, 134], 'float32');
fclose(cone_fid);
%------------------------------------------------------


% --- SIMULATE NEURON AND GENERATE SPIKE TRAIN ---
% simulation parameters
num_subunits = 5; % only used in 'subunit' model
desired_spike_count = 1e5;
verbose = true;

% image parameters
stimulus_width = 320;
stimulus_height = 320;
image_colors = 3;

% nonlinear model parameters
nl_model = 'exponential';  % can be "exponential" or "heaviside" nonlinearity
exp_params.a = 3.8394;
exp_params.b = -0.0277;
hvsd_params.thresh = 0.3;
hvsd_params.gain = 1;

% USE K-MEANS CLUSTERING ON CONE LOCATIONS TO SPECIFY SUBUNITS
% NOTE: THIS MAKES SUBUNITS NON-OVERLAPPING, WHICH MAY BE INCORRECT
locations = cone_centers(cone_indices,:);
[ids, centroids] = kmeans(locations, num_subunits);


% plot the clusters as a check on clustering and subunit structure -- clumsy code
if verbose
    figure(11); clf; 
    color_pal = {[1 0 1], [1 1 0], [0 1 1], [0 1 0], [1 0 0], [1 1 1], [0 0 1], [0.5 1 0], [1 0.5 0], [1 0.5 0.5], [0 1 0.5], [0 0.5 1]};
    rf_image = zeros((stimulus_width*stimulus_height*image_colors),1);
    for clust = 1:num_subunits
        clust_indices = find(ids == clust);
        sub_cone_rfs = cone_rfs(:,clust_indices);
        for rf = 1:size(sub_cone_rfs,2)
            clr = color_pal{mod(clust,length(color_pal))};
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
    nw_rf = reshape(rf_image, [stimulus_width stimulus_height image_colors]);
    imagesc(nw_rf)
    axis([ 110 180 150 220])
    axis square
    axis off

end


% initialize the matrix to store the generator signals for each subunit
sub_gen = zeros(size(rf_cone_inputs,1), num_subunits);

switch nl_model

    case 'exponential'
        
        % apply NL to output of subunits

        for clust = 1:num_subunits
            clust_indices = find(ids == clust);
            sub_gen_ = sum(rf_cone_inputs(:,clust_indices), 2);
            sub_gen_ = sub_gen_ ./ abs(ext(sub_gen_));
            sub_gen(:,clust) = exp(-exp_params.b + exp_params.a * sub_gen_);
        end

    case 'heaviside'

        for clust = 1:num_subunits
            % heaviside function
            clust_indices = find(ids == clust);
            sub_gen_ = sum(rf_cone_inputs(:,clust_indices), 2);
            sub_gen_ = sub_gen_ ./ abs(ext(sub_gen_));
            amped_indices = find(sub_gen_ >= hvsd_params.thresh);
            sub_gen(amped_indices, clust) = hvsd_params.gain * sub_gen_(amped_indices);

        end
end

% sum over subunit to generate an nut generator function
net_generator = sum(sub_gen, 2);

% Scale the generator signal such that approximately the disired number 
% of spikes are generated by the poisson process
mean_generator = mean(net_generator);
scale_factor = mean_generator ./ (desired_spike_count ./ length(net_generator))';
scaled_generator = net_generator ./ scale_factor;


% generate spike train from linear->nonlinear generator signal
poiss_train = poissrnd(scaled_generator);

% calculate and report the total number of spikes generated
total_spike_count = sum(poiss_train);
sprintf('total spike count is %d \n', total_spike_count)



% IDENTIFY SPIKES IN STIMULUS BINS
% translate to spike times (with duplicates for multiple spikes per time bin)
spike_times_ = [];
for nn = 1:max(poiss_train)
    spike_times_ = [spike_times_; find( poiss_train > (nn-1) )];
end
model_cell_spike_times = sort(spike_times_);


% CHECK THE STA
% compute the average set of frames that coincided with a spike
if verbose
    full_spike_frames = rf_cone_inputs(model_cell_spike_times,:);
    ave_spike_frame = mean(full_spike_frames, 1);
    % project spatial weights into cone RF space
    STA = all_cone_rfs(:,cone_indices) * ave_spike_frame';
    STA = STA - min(STA);
    STA = STA ./ max(STA) ./2;
    baseline = 0.5 * ones((stimulus_width * stimulus_height * image_colors),1);
    STA = STA + baseline;
    figure(1); clf;
    imagesc(reshape(STA, [stimulus_width, stimulus_height, image_colors]))
    axis([ 110 180 150 220])
    axis square
    axis off
end

%print(1, '~/Desktop/basic_sta.pdf', '-dpdf')

% get the stimulus frames that produced spikes across the cones that comprise
% the RGC receptive field
STE = rf_cone_inputs(model_cell_spike_times,:);

stim_cov = cov(rf_cone_inputs);

STA_ = mean(STE);
STC = STE - (STE*STA_') * STA_ ./ norm(STA_)^2;  % project out the mean
STC = cov(STC) - stim_cov;


% --- FACTORIZE THE STC ---
%[PCs, eig_vals_matrix] = eig(STC);
%eig_vals = diag(eig_vals_matrix);
%[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
%PCs = PCs(:, sorted_eig_indices);

%--- Princomp STE ---
[PCs, weights, eig_vals] = princomp(STE);


% --- PLOT EIGENVECTORS OF THE STC ---
num_dims = 6;
for dm = 1:num_dims
    figure(dm); clf;
    cone_egvec = all_cone_rfs(:,cone_indices) * PCs(:,dm); 
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    axis([ 110 180 150 220])
    axis square
    axis off

    %title_text = ['PC ', num2str(comp)];
    %title(title_text)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate noise eigenvalues by creating shifting spike times relative to stimulus
num_iters = 50;
time_offsets = -6-round(10000*rand(num_iters,1));
shifted_egvals = zeros(num_iters, length(cone_indices));
for iter = 1:num_iters

    % take spikes just from the relevant subset and time-shift to align peak frame with stimulus
    shifted_spike_rate_ = circshift(poiss_train,time_offsets(iter));

    % translate to spike times (with duplicates for multiple spikes per time bin)
    shifted_spike_times = [];
    for nn = 1:max(shifted_spike_rate_)
        shifted_spike_times = [shifted_spike_times; find( shifted_spike_rate_ > (nn-1) )];
    end
    shifted_spike_times = sort(shifted_spike_times);

    shifted_STE = rf_cone_inputs(shifted_spike_times,:);
    
    % princomp the STE
    [shifted_PCs, shifted_weights, shifted_egvals_] = princomp(shifted_STE);    

    % more careful massaging of the STE
%     shifted_STA_ = mean(STE);
%     shifted_STC = shifted_STE - (shifted_STE*shifted_STA_') * shifted_STA_ ./ norm(shifted_STA_)^2;  % project out the mean
%     shifted_STC = cov(shifted_STC) - stim_cov;
%     [shifted_PCs, shifted_eig_vals_matrix] = eig(shifted_STC);
%     shifted_egvals_ = diag(shifted_eig_vals_matrix);
%     [shifted_egvals_, shifted_sorted_eig_indices] = sort(shifted_egvals_, 'descend');
%     shifted_PCs = shifted_PCs(:, shifted_sorted_eig_indices);

    shifted_egvals(iter,:) = shifted_egvals_'./abs(sum(shifted_egvals_));

end

ave_egval = mean(shifted_egvals);
std_egval = std(shifted_egvals);

figure(20)
clf; 
plot(eig_vals./abs(sum(eig_vals)), 'r.','MarkerSize', 40)
hold on; 
plot(ave_egval+(3*std_egval), 'k')
plot(ave_egval-(3*std_egval), 'k')
xlabel('PCs')
ylabel('relative variance')
%title('red = data, black = +/- 3-sigma noise')
axis([0 20 0.008 0.01])

print(20,'~/Desktop/eig-val-spec.pdf', '-dpdf')

% --- PLOT EIGENVECTORS OF THE STC ---
num_dims = 6;
for dm = 1:num_dims
    figure(dm); clf;
    cone_egvec = all_cone_rfs(:,cone_indices) * PCs(:,dm); 
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    axis([ 110 180 150 220])
    axis square
    axis off

    %title_text = ['PC ', num2str(comp)];
    %title(title_text)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_pts = zeros(size(STE,2), num_subunits);
for i = 1:num_subunits
    temp_su_indices = find(ids == i);
    start_pts(temp_su_indices, i) = 0.05; 
end
[IDX, VQs] = kmeans(STE, 5, 'Start', start_pts');


[IDX, VQs] = kmeans(STE, 5);
kernels = VQs';

[IDX, VQs] = kmeans(weights(:,1:5), 5);
kernels = PCs(:,1:5) * VQs';


% --- PLOT CLUSTER CENTROIDS ---
num_dims = 5;
for dm = 1:num_dims
    figure(dm); clf;
    cone_egvec = all_cone_rfs(:,cone_indices) * kernels(:,dm); 
    cone_egvec = reshape(cone_egvec, [320, 320, 3]);
    imagesc(norm_image(cone_egvec))
    axis([ 110 180 150 220])
    axis square
    axis off

    %title_text = ['PC ', num2str(comp)];
    %title(title_text)
end





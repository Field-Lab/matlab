function [mean_spectrum, sd_spectrum] =  test_stc_significance(spike_times, cone_inputs, sig_dimensions, varargin)

% Inputs:
% spike_times, a vector of times in which a particular cell fired spikes
% cone_inputs, a matrix that is num cones by frames of cone generator
%               signals

p = inputParser;

p.addRequired('spike_times', @isnumeric);
p.addRequired('cone_inputs', @isnumeric);
p.addRequired('sig_dimensions', @isnumeric);

p.addParamValue('iternation_number', 10, @isnumeric);
p.addParamValue('Wc', [], @isnumeric);

p.parse(spike_times, sig_dimensions, cone_inputs, varargin{:})
Wc = p.Results.Wc;


num_iters = p.Results.iternation_number;


% get the random number of temporal offsets (shuffles) between spike train
% and stimulus
time_offsets = ceil(10000*rand(num_iters,1));

% Initialize large variable
shifted_eig_vals = zeros(num_iters, size(cone_inputs,2) - sig_dimensions - 1);

% compute the STA 
STE = cone_inputs(spike_times,:);

STA = mean(STE)';

% compute the STCs 
z_STE = STE' - STA*(STA'*STE')./ norm(STA)^2;  % project out the STA
z_cone_inputs = cone_inputs' - STA*(STA'*cone_inputs')./ norm(STA)^2;  % project out the STA

STC = cov(z_STE') - cov(z_cone_inputs');

% --- FACTORIZE THE STC ---
[PCs, eig_vals_matrix] = eig(STC);
eig_vals = diag(eig_vals_matrix);
[eig_vals, sorted_eig_indices] = sort(eig_vals, 'descend');
PCs = PCs(:, sorted_eig_indices);

% Shuffle the spike times and compute the eigenvalue spectum
for iter = 1:num_iters
    
    % circularly shift the spike times
    s_spike_times = 1 + mod(spike_times + time_offsets(iter), size(cone_inputs,1));
    
    % get the spike triggered ensemble (shifted)
    s_STE = cone_inputs(s_spike_times,:); 
    
    % project the STA out of the STE and the raw stimulus ensemble
    z_STE = s_STE';
    z_STE = z_STE - STA*(STA'*z_STE)./ norm(STA)^2;  % project out the mean
    z_cone_inputs = cone_inputs';
    z_cone_inputs = z_cone_inputs - STA*(STA'*z_cone_inputs)./ norm(STA)^2;  % project out the mean

    % project out the putatively signficant covariance dimensions
    for dm = 1:sig_dimensions
        z_STE = z_STE - PCs(:,dm)*(PCs(:,dm)'*z_STE)./ norm(PCs(:,dm))^2;
        z_cone_inputs = z_cone_inputs - PCs(:,dm)*(PCs(:,dm)'*z_cone_inputs)./ norm(PCs(:,dm))^2;
    end
    s_STC = cov(z_STE') - cov(z_cone_inputs');
    
    % --- FACTORIZE ---
    [s_PCs, s_eig_val_matrix] = eig(s_STC);
    s_eig_vals = diag(s_eig_val_matrix);
    [s_eig_vals, s_sorted_eig_indices] = sort(s_eig_vals, 'descend');
    s_PCs = s_PCs(:, s_sorted_eig_indices);  

    % get indices to meaningful dimensions
    keep_indices = find(abs(s_eig_vals) > 1e-10); % note this number might need to be changed depending on data
    % store eigvals;
    shifted_eig_vals(iter,:) = s_eig_vals(keep_indices); 
end

mean_spectrum = mean(shifted_eig_vals);
sd_spectrum = std(shifted_eig_vals);



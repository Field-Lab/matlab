function sig_str = pop_motion_signal_no_integral(velocity, spikes, indices1, indices2, x_pos, trigger, trial_length, tau, tol)
%POP_MOTION_SIGNAL calculates the cumulative motion signal for a given velocity between all pairs in population of neurons.
%   sig_str = POP_MOTION_SIGNAL(velocity, spikes, indices1, indices2, x_pos, trigger, stop, tau)
%   caluclates the signal strength for the given velocity where:
%
%       spikes -> is a [num_neurons x 1] cell containing [num_spikes x 1] 
%                 vectors corresponding to all the spike times of that
%                 specific cell.  Contains all the neurons picked up by the
%                 array.
%       indices1 -> contains the indices (for x_pos) for the neurons
%                   identified using white noise to check receptive fields
%       indices2 -> contains the indices (for spikes) to identify the
%                   neurons during each trial run.  indices1(i) and
%                   indices2(i) correspond to the same neuron.
%       x_pos -> is a [num_neurons2 x 1] vector containing the x position
%                of all the neurons identified in the white noise stimuli.
%                i.e. the x coordinate of the gaussian fit of the receptive
%                field.
%       trigger ->the start time for the trial you are using
%       trial_length -> length of the trial
%       tau -> the tuning parameter to use.  The width of the gaussians in
%              the filtered response algorithm
%
%   sig_str = POP_MOTION_SIGNAL(..., tol)
%   same as above with:
%       tol -> tolerance to use for the integration algorithm in matlab.
%              Defaults to 1e-3
%
%   For examples on usage see motion_script.m
%
% Marvin Thielk 2013
% mthielk@salk.edu
if nargin < 9
    tol = 1e-3;
end

sig_str = 0;

% create a list of all possible pairs between the neurons
pairs = zeros(2, length(indices2) * (length(indices2) - 1) / 2);
counter = 1;
for i = 2:length(indices2)
    for j = 1:i-1
        pairs(:,counter) = [i; j];
        counter = counter + 1;
    end
end

for j = 1:length(pairs)
    spks_1 = spikes{indices2(pairs(1,j))};
    spks_2 = spikes{indices2(pairs(2,j))};
    dx = x_pos(indices1(pairs(2,j))) - x_pos(indices1(pairs(1,j)));
    sig_str = sig_str + motion_signal_no_integral(velocity, spks_1, spks_2, dx, trigger, trial_length, tau, tol);
end

end
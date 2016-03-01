function [sig_str, spks_out] = pop_motion_signal_getSpikes(velocity, spikes, indices1, indices2, x_pos, trigger, start, trial_length, tau, tol, datarun, direction, stx)
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
%
% Modified for speed by Malcolm Campbell 2014
% malcolmc@stanford.edu
%
% Modified for faster motion signal computation by Colleen Rhoades 2014
% rhoades@stanford.edu


if nargin < 9
    tol = 1e-3;
end

sig_str = 0;

% create a list of all possible pairs between the neurons
RelTol = tol;
AbsTol = tol * .1;
trial_length = trial_length - start;
numberOfT = trial_length*138; % determined to be good enough
t = linspace(start, trial_length, numberOfT);
valueAtEachT = zeros(length(t),2);


    % reference is set to 0 (still works if the 0 is changed to another
    % number)
    dx = x_pos(indices1)*stx/10;




spks_2 = spikes(indices2); % Cell array of spike trains of all types considering

% trial start @ t=0

spks_2=cellfun(@(x) x-trigger, spks_2,'UniformOutput',false);

% only consider spikes that occured in the trial
spks_2 = cellfun(@(y) y(y >= start & y <= start+trial_length), spks_2, 'UniformOutput', false);

% circularly shift spikes by dt

dx_cell = num2cell(dx);

if strcmp(direction, 'right') % bar moving right
spks_2_shiftedRight = cellfun(@(z, c) z - c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
spks_out = spks_2_shiftedRight;
spks_2_shiftedRight = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedRight, 'UniformOutput', false);
spks_2_shiftedLeft = cellfun(@(z, c) z + c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
spks_2_shiftedLeft = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedLeft, 'UniformOutput', false);

else % bar moving left
spks_2_shiftedRight = cellfun(@(z, c) z + c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
spks_out = spks_2_shiftedRight;
spks_2_shiftedRight = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedRight, 'UniformOutput', false);
spks_2_shiftedLeft = cellfun(@(z, c) z - c / velocity, spks_2,dx_cell, 'UniformOutput', false); %right shift
spks_2_shiftedLeft = cellfun(@(z) sort(mod(z, trial_length)), spks_2_shiftedLeft, 'UniformOutput', false);
end

% replicate spikes before and after trial to minimize artifacts of spikes
% shifting circularly across the border
idx=cellfun('isempty',spks_2_shiftedLeft);
spks_2_shiftedLeft(idx)={0}; %It replaces all empty cells with number 0

idx=cellfun('isempty',spks_2_shiftedRight);
spks_2_shiftedRight(idx)={0}; %It replaces all empty cells with number 0



% changed both sides to use ceil
spks_2_shiftedRight = cellfun(@(z) [z(ceil(end/2):end) - trial_length; z; z(1:ceil(end/2)) + trial_length], spks_2_shiftedRight, 'UniformOutput', false);
spks_2_shiftedLeft = cellfun(@(z) [z(ceil(end/2):end) - trial_length; z; z(1:ceil(end/2)) + trial_length], spks_2_shiftedLeft, 'UniformOutput', false);

% filter responses
flt_rsp2 = cellfun(@(x) filtered_response(x, tau), spks_2, 'UniformOutput', false);
flt_rsp2_shiftedRight = cellfun(@(x) filtered_response(x, tau), spks_2_shiftedRight, 'UniformOutput', false);
flt_rsp2_shiftedLeft = cellfun(@(x) filtered_response(x, tau), spks_2_shiftedLeft, 'UniformOutput', false);

for c= 1:length(t)
    time = t(c);
    % First column is every cell aligned to cell 1
    % First row is value of first pairing at t = 0
    [right, left] = summedCellsAtT(time, flt_rsp2_shiftedRight,  flt_rsp2_shiftedLeft, indices1,spks_2);
    valueAtEachT(c,1) = right;
    valueAtEachT(c,2) = left;
end


% sum from 2:end because first and last value are identical!
sig_str = sum(valueAtEachT(2:end,1).^2)-sum(valueAtEachT(2:end,2).^2);


end
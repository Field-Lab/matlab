function [ccf, time] = compute_ccf(s1, s2, options)
% COMPUTE_CCF     Calculate cross-correlation between two spike trains.
%
% usage:   [ccf, time] = compute_ccf(s1, s2);
%          [ccf, time] = compute_ccf(s1, s2, options)
%
% arguments:    s1, s2 - row vector of spike times from 2 neurons
%              options -  structure of options
%                         [ dt       -  bin size (s) [offset/128]
%                         [ offset   -  maximum offset for CCF (s) [100e-3]
%                         [ shuffle  -  string to determine type of shuffle
%                                        'none'  - standard CCF   [default]
%                                        'stim' - signal correlations
%                                        'noise' - noise correlations
%                                           Must include 'trial' if latter 
%                                           two options are used.
%                         [ trial    -  trial duration of shuffle
%
% outputs:    ccf - 1xN cross correlation between two neurons
%            time - 1xN time offsets for cross correlations
%
% example:
%  options = struct('offset',100e-3,'scale','ms','shuffle','none');
%  [ccf, time] = compute_ccf(spikes{1},spikes{2},options);
%
% See also: PLOT_CCF, COMPUTE_ACF, PLOT_ACF
%
% jds 2006-07-17
% jds 2005-09-13
%

% default number of samples in CCF
N = 63;
% sampling frequency of spike trains
sampling_frequency = 20000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% default arguments
default_options = struct('offset', 100e-3, 'shuffle', 'none', 'trial', 1);
if (nargin == 2)
  options = default_options;
end

% check if individual items are not specified
if ~isfield(options,'offset')
  options.offset = default_options.offset;
end
if ~isfield(options,'dt')
  options.dt = options.offset / N;
end
if ~isfield(options,'shuffle')
  options.shuffle = default_options.shuffle;
end
if ~isfield(options,'trial')
  options.trial = default_options.trial;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. TIME CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert to integer sampling
s1 = int32(s1 * sampling_frequency);
s2 = int32(s2 * sampling_frequency);

% convert time to samples
offset = int32(options.offset * sampling_frequency);
dt = int32(options.dt * sampling_frequency);
trial = int32(options.trial * sampling_frequency);

% determine time
time = [-offset:dt:offset]';
time = double(time) / sampling_frequency;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. SELECT CALCULATION TYPE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for speed, do each calculation individually

s1 = reshape(s1, length(s1),1);
s2 = reshape(s2, length(s2),1);

switch options.shuffle
 case 'none'
  ccf = xcorr_sparse(s1, s2, -offset, offset, dt);
 case 'stim'
  ccf = ccf_shuffle(s1, s2, offset, dt, trial);
 case 'poisson'
  ccf = ccf_poisson(s1, s2, offset, dt, trial);
 case 'noise'
  ccf_standard = xcorr_sparse(s1, s2, -offset, offset, dt);
  ccf_stim = ccf_shuffle(s1, s2, offset, dt, trial);
  ccf_poisson = ccf_poisson(s1, s2, offset, dt, trial);

  ccf = ccf_standard - ccf_stim + ccf_poisson;

 otherwise
  error('compute_ccf: unrecognized shuffle type');
end

% convert to double
ccf = double(ccf);

% convert to firing rate
ccf = ccf / (length(s1) * options.dt);

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccf = ccf_poisson(s1, s2, offset, dt, trial_duration);

duration = int32(max(s1));

% circularly shift spike train (to destroy all correlation)
s1 = sort(rem(s1 + int32(trial_duration/2), duration));

% compute ccf
ccf = xcorr_sparse(s1, s2, -offset, offset, dt);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ccf = ccf_shuffle(s1, s2, offset, dt, trial_duration);


%% 1. SET UP SHUFFLING
%%%%%%%%%%%%%%%%%%%%%%

% shuffling
if isempty(trial_duration)
  error('compute_ccf: need trial duration for shuffle');
end

% approximate duration of spike trains 
duration = max([max(s1) max(s2)]);
duration = ceil(duration / dt) * dt;

% determine number of trials
trials = floor(duration / trial_duration);

% just need a different index for shuffling
shuffle_trials = circshift([1:trials]',1);


% determine trials for each set of spike times
spikes1 = mod(s1, trial_duration);
trials1 = int32(floor(double(s1) / double(trial_duration)) + 1);
spikes2 = mod(s2, trial_duration);
trials2 = int32(floor(double(s2) / double(trial_duration)) + 1);
clear s1 s2

%% 2. CALCULATE CCF
%%%%%%%%%%%%%%%%%%%


% determine time
time = [-offset:dt:offset]';

% calculate ccf for each trial
ccf = int32(zeros(size(time)));

% calculate CCF in chunks (or trials)
for i=1:trials
  
  % find spikes within trial
  index1 = find(trials1 == i);
  index2 = find(trials2 == shuffle_trials(i));

  % error check
  if any(diff(spikes1(index1)) < 0) |  any(diff(spikes2(index2)) < 0)
    error('compute_ccf: round off error');
  end
  
  % compute cross-correlation
  ccf = ccf + xcorr_sparse(spikes1(index1), spikes2(index2), ...
			   -offset, offset, dt);  
end

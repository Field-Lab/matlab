function tai = compute_ta_ind(st, ter, tparams)
%COMPUTE_TA_IND Computes the indices of the triggered average of some 
%quantity.
% 
%  TAI = COMPUTE_TA_IND(ST, TIME_EVENTS_REFRESH) assumes that a series of
%  events happened such that the k-th event took place between the
%  TIME_EVENTS_REFRESH(k) and TIME_EVENTS_REFRESH(k+1)-1 samples of the
%  recording. Typically, these events would be the display of an image on a
%  monitor, for example. 
%  Then, for a spike time series ST, this function will return an n x k
%  matrix TAI such that the spike-triggered quantity associated with the 
%  events can be calculated by averaging each of the events with indices
%  TAI(:, kk) with 1 <= kk <= k.
%
%  TAI = COMPUTE_TA_IND(..., TPARAMS) lets you further speficy how far into
%  the future/past and with what time granularity the triggered average
%  should be computed. TSTEPS should be a 3x1 vector such that 
%  TPARAMS = [T_PRE, T_POST, STEP_SZ] with T_PRE the number of samples
%  preceeding each spike used for the computation (defaults to 1500), 
%  T_POST the number of samples following each spike (defaults to 0) and
%  STEP_SZ the STA step (defaults to 20 samples).

if nargin == 2
    t_pre = 1500;
    t_post = 0;
    step_sz = 20;
else
    assert(numel(tparams) == 3);
    t_pre = tparams(1);
    t_post = tparams(2);
    step_sz = tparams(3);
end
t_steps = (-t_pre):step_sz:t_post;

% Trim the spike train: remove spikes for which calculating in the TA
% would make us look at parts of the recording that took place before
% it started
st = st(st > t_pre);

% Reshape st and t_steps matrices
if size(st, 1) < size(st, 2)
    st = st.';
end
if size(t_steps, 1) > size(t_steps, 2)
    t_steps = t_steps.';
end

% Now we can fill in the triggered average indices matrix
% Start by calculating the time in samples of each of the events of
% interest
tai = repmat(st, 1, length(t_steps)) + repmat(t_steps, length(st), 1);

% Then, convert those times in samples to event indices
tai = arrayfun(@(x)(sum(ter<x)), tai);

end % compute_ta_indices
function cone_inputs = calc_cone_inputs(Wc, stimindices, stimobject, get_movie_frame)
% CALC_CONE_INPUTS      Reexpress movie in terms of cone signals
% usage: cone_inputs = calc_cone_inputs(Wc, stimindices, stimobject, get_movie_frame)
%
% WC is as saved out of cone finding.  STIMINDICES indices which frames of
% the stimulus to calculate on.  STIMOBJECT and GET_MOVIE_FRAME are a
% complementary pair: STIMOBJECT holds the stimulus frames, GET_MOVIE_FRAME
% is a function handle that gets the desired frame from STIMOBJECT.
%
% It may be beneficial to preprocess WC to match the native frame format of
% STIMOBJECT, rather than build a lot of processing into GET_MOVIE_FRAME
% that must be run on every movie frame.
%
% 2012-09, phli, abstracted out of cone_preprocess
%

% cycle through each stimulus
fprintf('Computing cone input in frames %d to %d...', stimindices(1), stimindices(end))
T = text_waitbar();
cone_inputs = zeros(length(stimindices), size(Wc,2), 'single');
for ss = 1:length(stimindices)
    T = text_waitbar(T, ss/length(stimindices) - 0.01);
    
    % note which stim
    this_stim = stimindices(ss);
    
    % get new frame
    new_frame = get_movie_frame(stimobject, this_stim);
    
    % convert to cone space
    cone_inputs(ss,:) = single((Wc' * new_frame) - 0.5);
end
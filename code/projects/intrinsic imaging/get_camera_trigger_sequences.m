function seqs = get_camera_trigger_sequences(triggers, varargin)
% GET_CAMERA_TRIGGER_SEQUENCES    Tries to parse triggers from the camera setup into useful groups
% usage: seqs = get_camera_trigger_sequences(triggers, [opts])
%
% When we record triggers from the camera, we typically end up with an
% experiment with sequences of triggers.  For example, the triggers might
% have 1000 fast, evenly spaced triggers, followed by a few scattered
% artifact trigger events, followed by another 1000 fast, evenly spaced
% triggers.  The goal of this function is to pull out the two useful trigger
% seqences.
%
% Note that the triggers from our SPOT cameras are not perfectly accurate;
% there is often several milliseconds of jitter.
%
% inputs: triggers
% options: accept_factor  []    How much longer (mutiplicative) an interval can
%                               be compared to the shortest interval before we
%                               detect it as a break in the sequence.  If
%                               blank, uses the accept_offset instead.
%
%          accept_offset  0.1   How much longer (additive) an interval can
%                               be before it is detected as a break in the
%                               sequence.  Overridden by accept_factor.
%
%          min_sequence   2     The minimum sequence length.  Sequences
%                               shorter than this are assumed to be
%                               artifactual.
%
% 2010-06 phli
%

opts = inputParser;
opts.addParamValue('accept_factor', []);
opts.addParamValue('accept_offset', 0.1);
opts.addParamValue('min_sequence', 5);
opts.parse(varargin{:});
opts = opts.Results;


intervals = diff(triggers);
min_interval = min(intervals);

if ~isempty(opts.accept_factor)
    accept_interval = min_interval * opts.accept_factor;
else
    accept_interval = min_interval + opts.accept_offset;
end
long_intervals = find(intervals > accept_interval);


sequence_endpoints = [1; (long_intervals(:) + 1); (length(triggers) + 1)];
seqs = {};
for i = 1:(length(sequence_endpoints) - 1)
    seq_start = sequence_endpoints(i);
    seq_end   = sequence_endpoints(i+1) - 1;
    
    % Useful for detecting two long intervals next to each other
    if seq_end - seq_start < opts.min_sequence
        continue;
    end
    
    seqs{end+1} = triggers(seq_start:seq_end); %#ok<AGROW>
end
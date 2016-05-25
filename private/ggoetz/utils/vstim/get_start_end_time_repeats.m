function [ time_repeats ] = get_start_end_time_repeats(fo, tf)
%GET_START_END_TIME_REPEATS calculates the start and end time of repeats in
%a movie.
%
%  [TIME_REPEATS] = GET_START_END_TIME_REPEATS(FO, TF) returns a nRepeats 
%   x 2 matrix whose first column contains the start times, and second
%   column the end times of the repeated frames in the frame sequence FO,
%   whose frames were displayed at times TF.

% Find which frames were repeats
fo_sorted = sort(fo);
% Frame 0 is special, remove it
fo_sorted(fo_sorted == 0) = [];
repeat_frames = zeros(length(fo_sorted)-1, 1);
for kk=1:(length(fo_sorted)-1)
    if (fo_sorted(kk+1) == fo_sorted(kk))
        repeat_frames(kk) = fo_sorted(kk);
    end
end
repeat_frames(repeat_frames == 0) = [];
repeat_frames = unique(repeat_frames);
is_repeat = double(ismember(fo, repeat_frames));

% Build the repeats matrix. Should be vectorized, someday...
time_repeats = [];
fi = 1;
while fi < length(tf)
    % find start of the repeat
    if (is_repeat(fi+1) - is_repeat(fi)) > 0
        start_repeat_time = tf(fi+1);
        
        % find end of the repeat
        while ~(is_repeat(fi+1) - is_repeat(fi) < 0) && (fi < length(tf))
            fi = fi +1;
        end
        end_repeat_time = tf(fi);
        
        % store repeat
        time_repeats = [time_repeats; start_repeat_time end_repeat_time];
    end
    fi = fi + 1;
end

end % get_start_end_time_repeats
function [ st ] = remove_repeats_from_st(tr, st)
%REMOVE_REPEATS_FROM_ST Removes the spikes that occured during a repeat
%from a spike train.
%
%  [ST_NEW] = REMOVE_REPEATS_FROM_ST(TIME_REPEATS, ST) trims the spike 
%  train ST according to a nx2 matrix TIME_REPEATS, whose first column
%  contains the time at which repeats start and second column the time at
%  which they end.

for kk=1:size(tr, 1)
    start_repeat_time = tr(kk, 1);
    end_repeat_time = tr(kk, 2);
    st(logical((st > start_repeat_time).*(st < end_repeat_time))) = [];
end

end % remove_repeats_from_st


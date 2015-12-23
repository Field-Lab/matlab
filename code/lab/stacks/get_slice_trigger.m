function trig = get_slice_trigger(stack, islice)
% GET_SLICE_TRIGGER
% usage: trig = get_slice_trigger(stack, islice)
%
% Get time of image trigger, usually relative to beginning of datarun.
% Note that the triggers from our SPOT cameras are not perfectly accurate;
% there is often several milliseconds of jitter.
%
% 2010-08 phli
%

if isfield(stack, 'triggers') && length(stack.triggers) >= islice
    trig = stack.triggers(islice);
else
    trig = [];
end
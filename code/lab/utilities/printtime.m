function p = printtime(rawsecs)

secs = mod(rawsecs,60);
rawmins = floor(rawsecs/60);

mins = mod(rawmins,60);
rawhours = floor(rawmins/60);

hours = mod(rawhours,24);
rawdays = floor(rawhours/60);

days = mod(rawdays,7);
rawweeks = floor(rawdays/7);

weeks = rawweeks;


if rawweeks > 0
    p = sprintf('%dW %dD %dH', weeks, days, hours);
elseif rawdays > 0
    p = sprintf('%dD %d:%d', days, hours, mins);
elseif rawmins > 0
    p = sprintf('%d:%02d:%02d', hours, mins, secs);
else
    p = sprintf('%d s', secs);
end
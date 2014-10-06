function max_frame = get_ei_max_frame(ei, start_frame, end_frame, negative)

if nargin < 5
    logscale=0;
end

if nargin < 4
    negative = true;
end

if nargin < 3
    end_frame = size(ei,2);
end

if nargin < 2
    start_frame = 1;
end


max_frame = zeros(size(ei,1),1);
for e = 1:size(ei,1)
    if ~negative
        [t,t] = max(abs(ei(e, start_frame:end_frame)));
        max_frame(e) = ei(e,t);
    else
        max_frame(e) = min(ei(e, start_frame:end_frame));
    end
    
end

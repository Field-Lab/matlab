function [firstposes, secondposes, firstnegs, negs, elec_colors] = flipflop_ei_cat(ei, varargin)

opts = inputParser;
opts.addParamValue('neg_thresh', -10);
opts.addParamValue('neg_width', 4);
opts.addParamValue('pos_width', 1);
opts.addParamValue('pos_width2', 2);
opts.parse(varargin{:});
opts = opts.Results;


baselines = mean(ei(:,1:10)');
tracelen = size(ei, 2);
ei_base = ei - repmat(baselines(:), [1 tracelen]);

[mins, minlocs] = min(ei_base');

negs = mins < opts.neg_thresh;
firstneg = min(minlocs(negs));
firstnegs = find(negs & (minlocs < firstneg+opts.neg_width));

poslocs = [];
for i = firstnegs
    trace    = ei_base(i,:);
    traceloc = 1:length(trace);
    
    postpos = find(traceloc > minlocs(i) & trace > 0);
    poslocs(end+1) = postpos(1);
end
firstpos = min(poslocs);
firstposes = firstnegs(poslocs < firstpos+opts.pos_width);
secondposes = firstnegs(poslocs < firstpos+opts.pos_width2);


% ToDo-PHL-2010-05 Move this to separate function...
numelec = size(ei, 1);
elec_colors = repmat([0 0 0], [numelec, 1]);
elec_colors(negs, 2) = 1;
elec_colors(firstnegs,:)  = repmat([0 0 1], [length(firstnegs),  1]);
elec_colors(secondposes,:) = repmat([1 0 1], [length(secondposes), 1]);
elec_colors(firstposes,:) = repmat([1 0 0], [length(firstposes), 1]);
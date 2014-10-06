function [poslocs, pos_tracex, pos_tracey] = flipflop_ei_cat2(ei, positions, varargin)

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

negs = find(mins < opts.neg_thresh);

poslocs = cell([1 tracelen]);
for i = negs
    trace    = ei_base(i,:);
    traceloc = 1:length(trace);
    
    postpos = find(traceloc > minlocs(i) & trace > 0);
    if ~isempty(postpos)
        poslocs{postpos(1)}(end+1) = i;
    end
end

pos_tracex = [];
pos_tracey = [];
% for i = 1:length(poslocs)
%     if ~isempty(poslocs{i})
%         for j = 1:length(poslocs{i})
%             electrode = poslocs{i}(j);
%             xposes(j) = positions(electrode, 1);
%             yposes(j) = positions(electrode, 2);
%         end
%         pos_tracex(end+1) = mean(xposes);
%         pos_tracey(end+1) = mean(yposes);
%     end
% end
for i = 1:length(poslocs)
    if ~isempty(poslocs{i})
        pos_tracex(:,i) = NaN;
        pos_tracey(:,i) = NaN;
        for j = 1:length(poslocs{i})
            electrode = poslocs{i}(j);
            pos_tracex(j,i) = positions(electrode, 1);
            pos_tracey(j,i) = positions(electrode, 2);
        end
    end
end

for i = 1:numel(pos_tracex)
    if pos_tracex(i) == 0 && pos_tracey(i) == 0
        pos_tracex(i) = NaN;
        pos_tracey(i) = NaN;
    end
end
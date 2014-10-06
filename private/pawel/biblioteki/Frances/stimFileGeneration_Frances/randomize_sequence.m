function newseq = randomize_sequence(oldseq)
% function newseq = randomize_sequence(oldseq)
% This function takes an array oldseq and randomize its order into newseq.
%
%   Frances Yeh 2012/11/29

newseq = zeros(size(oldseq));

for i = 1:length(oldseq)
    idx = ceil(rand*length(oldseq));
    newseq(find(newseq == 0, 1)) = oldseq(idx);
    oldseq(idx) = [];
end

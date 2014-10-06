function st = refracst(st, refrac, jitter)
% Cut spikes in st that are less than refrac after the preceding spike.
%
% Assumes ST is sorted.  Does only one pass, so may cut some spikes later
% in a train that would not be cut if the 2nd spike were taken out in a
% first pass...
%

dst = diff(st(:));

if nargin > 2 && jitter > 0
    dst = dst + jitter*randn(size(dst));
end

keep = [true; dst > refrac];
st = st(keep);
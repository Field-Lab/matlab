function [st comb] = refracst_ro(ro, st, refrac, jitter)
% Cut spikes in ST that are less than REFRAC after the preceding spike.
% First combines spikes in ST with spikes in RO, but RO is read-only so
% those spikes cannot be removed.
%
% Assumes ST is sorted.
%


nro = length(ro);
nst = length(st);
comb  = [ro(:);        st(:)];
combi = [zeros(1,nro), 1:nst]';
[comb sorti] = sort(comb);
combi = combi(sorti);

% First cut all generated spikes that come after another spike by too
% little.  Don't cut true spikes though!  Keep track of which generated
% (versus true spikes) were cut using combinedindices
dst = diff(comb);
if nargin > 3 && jitter > 0
    dst = dst + jitter*randn(size(dst));
end
keep = [true; dst > refrac];
keep(combi == 0) = true; % True spikes are read only!
comb  = comb(keep);
combi = combi(keep);

% Since we didn't cut true spikes, there will still be refractory
% violations.  Now cut any generated spikes that precede true spikes by
% too little; this should now eliminate all violations.
dst = diff(comb);
if nargin > 3 && jitter > 0
    dst = dst + jitter*randn(size(dst));
end
keep = [dst > refrac; true];
keep(combi == 0) = true;
comb  = comb(keep);
combi = combi(keep);

combi = setdiff(unique(combi), 0);
st = st(combi);

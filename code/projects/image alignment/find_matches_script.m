
% get length of each axon
axon_lengths = [];
for aa=1:length(axons);
    axon_lengths(aa)=sum(sqrt(sum(diff(axons{aa}).^2,2)),1);
end

% only look for matches for sufficiently long axons
long_axons = find(axon_lengths>120);

for aa=1:length(long_axons)
    find_matches(datarun,axons,long_axons(aa));
    pause
end

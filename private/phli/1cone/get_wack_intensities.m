function stimstruct = get_wack_intensities(stimstruct)
% GET_WACK_INTENSITIES      Loads intensities based on info in STIMSTRUCT
%
% For the old style wack, the one that was really wack because we used
% white noise.  This was only run a few times...
%
% STIMSTRUCT should have MAPIM field with a single map entry; this is used 
% to determine the number of regions to draw for.  STRIMSTRUCT must also 
% have a SEEDS field listing the seed for each trigger.
%
% 2011-09 phli
%

regions = setdiff(unique(stimstruct.mapims{1}), 0);
seeds = stimstruct.seeds;

stimstruct.intensities = zeros(length(regions), length(seeds));
for s = 1:length(seeds)
    seed = seeds(s);
    rng = edu.ucsc.neurobiology.vision.math.RandomJavaV2(seed);
    for r = 1:length(regions)
        stimstruct.intensities(r,s) = rng.nextBits(1);
    end
end
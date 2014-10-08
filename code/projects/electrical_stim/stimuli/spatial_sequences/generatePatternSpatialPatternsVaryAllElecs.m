function array = generatePatternSpatialPatternsVaryAllElecs(electrodes, fullAmps)
%
% pattern 1: all electrodes stim'd stimultaneously
% patterns 2:nElec+1: each electrode stim'd individually
%
% generates array of values representing relative amplitudes of pulses on each electrode,
% with first dimension corresponding to electrodes (first electrode = primary) and second dimension
% corresponding to different pattern numbers
%
% electrodes: vector of electrode numbers

nElec = length(electrodes);
array = zeros(nElec, nElec+1);

for j = 1:nElec
    array(j, 1) = fullAmps(j); % simultaneous stimulation
    array(j, j+1) = fullAmps(j); % individual stimulation
end
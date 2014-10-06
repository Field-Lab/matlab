function [array variedElec] = generatePatternSpatialPatternsVary1Elec(electrodes, fullAmps, threshAmps)
%
% pattern 1: all electrodes stim'd stimultaneously
% patterns 2-nElec+1: each electrode stim'd individually
%
% generates array of values representing relative amplitudes of pulses on each electrode,
% with first dimension corresponding to electrodes (first electrode = primary) and second dimension
% corresponding to different pattern numbers
%
% fullAmps: amplitudes that result in 100% response of corresponding cell
% threshAmps: amplitudes that result in 50% response of corresponding cell

nElec = length(electrodes);

if size(fullAmps, 1) == 1
    fullAmps = fullAmps';
end

%array = zeros(nElec, nElec+9); %four extra simultaneous, four extra single elec

varyAmps = cell(1,nElec);
for ii = 1:nElec
    varyAmps{ii} = [0.81*threshAmps(ii) 0.9*threshAmps(ii) threshAmps(ii) 1.1*threshAmps(ii) 1.21*threshAmps(ii)];
end

%individual electrodes at full amplitudes
array = diag(fullAmps);
variedElec = zeros(1,nElec);

for ii = 1:nElec %loops through which electrode is being varied
    %simultaneous stimulation
    for kk = 1:5 %each near-thresh amp
        array = [array fullAmps];
        array(ii, end) = varyAmps{ii}(kk);
        variedElec = [variedElec ii];
    end
    
    %varied electrode alone
    for kk = 1:5
        array = [array zeros(nElec,1)];
        array(ii, end) = varyAmps{ii}(kk);
        variedElec = [variedElec ii];
    end
end

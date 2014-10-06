function [electrodes array details] = generatePatternSpatiotempProbeNew(centerElecs, neighborRadius, amps, varargin)
%
% generates array of values representing relative amplitudes of pulses on each electrode,
% with first dimension corresponding to electrodes (first electrode = primary) and second dimension
% corresponding to different pattern numbers
%
% electrodes: vector of electrode numbers
% neighborRadius: maximum distance of prepulse elec from center electrode,
% in units of nearest-neighbor distance
% amps    .neighborAmps{clusterIndex}       set of amplitudes of pulses applied to all neighboring electrodes
%                                           not specified in 'specificSecondaryAmps'
%         .centerAmpPrepulses{clusterIndex} set of amplitudes of pulses applied to center electrode
%                                           when being used as a prepulse
%         .centerAmpsLims{clusterIndex}     [lowerAmpLimit upperAmpLimit] defines range of amplitudes
%                                           to be used on stimulating (center) electrode
%         .centerIncrement(clusterIndex)    defines amplitude step between subsequent amplitudes to
%                                           be used on stimulating electrodes
%
%
%


p = inputParser;

p.addRequired('centerElec', @isnumeric)
p.addRequired('neighborRadius', @isnumeric)
p.addRequired('amps', @isstruct)

p.addParamValue('specificSecondaryAmps', [], @iscell) %overrides secondary amplitudes specified in 'amps' on specified electrodes
% specify secondary amplitudes in this format: specificSecondaryAmps{clusterIndex} = [electrode1 electrode2 electrode3;
%                                                                                     amplitude  amplitude  amplitude;
%                                                                                     amplitude  amplitude  amplitude]


p.parse(centerElecs, neighborRadius, amps, varargin{:})

specifySecondaryAmps = p.Results.specificSecondaryAmps;




%%

neighborRadius = neighborRadius*2; %because in getElectrodeCoords61, the nearest-neighbor distance is 2


% calculates center amplitudes
centerAmps = cell(1, length(centerElecs));
for i = 1:length(centerElecs)
    centerAmps{i} = amps.centerAmpLims{i}(1);
    while centerAmps{i}(end) < amps.centerAmpLims{i}(2)
        centerAmps{i} = [centerAmps{i} centerAmps{i}(end)*(1 + amps.centerIncrement(i))];
    end
end


% finds clusters of center + surrounding electrodes

[xCoords yCoords] = getElectrodeCoords61();
clusters = cell(length(centerElecs), 1);
electrodes = [];
for i = 1:length(centerElecs)
    distances = zeros(64, 1);
    for j = 1:64
        distances(j) = norm([xCoords(centerElecs(i)) - xCoords(j), yCoords(centerElecs(i)) - yCoords(j)]);
    end
    clusters{i} = find(squeeze(distances)<= neighborRadius*1.01); %includes center electrode and neighbors within neighborRadius
    clusters{i}(clusters{i}==centerElecs(i)) = [];
    electrodes = [electrodes; centerElecs(i); clusters{i}]; %#ok<AGROW>
end
nElec = length(electrodes);

%checks for overlap between clusters
for i = length(electrodes):-1:1
    if sum(electrodes == electrodes(i)) ~= 1
        warndlg('there is some cluster overlap -- check to make sure centers aren''t too close to eachother')
        electrodes(i) = [];
       
        %repInd = find(electrodes == electrodes(i));
        %electrodes(repInd(2:end)) = [];
    end
end



array = [];
details.prePulse = [];
details.stimPulseAmps = [];
details.clusterID = [];
details.pairsCenterAmps = [];
for i = 1:length(centerElecs) %should change to using concatenatePatterns
    pInd = find(electrodes==centerElecs(i));
    sInds = [];
    for j = 1:length(clusters{i})
        sInds = [sInds find(electrodes==clusters{i}(j))]; %#ok<AGROW>
    end
    
    %center alone, all amplitudes
    for j = 1:length(centerAmps{i})
        array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        array(pInd, end) = centerAmps{i}(j);
        details.prePulse = [details.prePulse 0];
        details.stimPulseAmps = [details.stimPulseAmps j];
        details.pairsCenterAmps = [details.pairsCenterAmps 0];
        details.clusterID = [details.clusterID i];
    end
    
    %center alone, prepulse amplitudes
    for j = 1:length(amps.centerAmpPrepulses{i})
        array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
        array(pInd, end) = amps.centerAmpPrepulses{i}(j);
        details.prePulse = [details.prePulse 1];
        details.stimPulseAmps = [details.stimPulseAmps 0];
        details.pairsCenterAmps = [details.pairsCenterAmps 0];
        details.clusterID = [details.clusterID i];
    end
    
    %neighbor alone, each amplitude
    for j = 1:length(sInds)
        if ~isempty(specifySecondaryAmps{i})
            specSecElec = find(specifySecondaryAmps{i}(1,:) == electrodes(sInds(j)));
        else
            specSecElec = [];
        end
        if ~isempty(specSecElec)
            ampsToUse = specifySecondaryAmps{i}(2:end, specSecElec);
            ampsToUse = ampsToUse(ampsToUse ~= 0);
            for k = 1:length(ampsToUse)
                array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                array(sInds(j), end) = ampsToUse(k);
                details.prePulse = [details.prePulse 1];
                details.stimPulseAmps = [details.stimPulseAmps 0];
                details.pairsCenterAmps = [details.pairsCenterAmps 0];
                details.clusterID = [details.clusterID i];
            end
        else %use amplitudes specified in amps.neighborAmps{i}
            for k = 1:length(amps.neighborAmps{i})
                array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                array(sInds(j), end) = amps.neighborAmps{i}(k);
                details.prePulse = [details.prePulse 1];
                details.stimPulseAmps = [details.stimPulseAmps 0];
                details.pairsCenterAmps = [details.pairsCenterAmps 0];
                details.clusterID = [details.clusterID i];
            end
        end
    end
    
    %center + neighbor
    for j = 1:length(centerAmps{i})
        %jth center amplitude, all neighbors and all neighbor amplitudes
        for k = 1:length(sInds)
            if ~isempty(specifySecondaryAmps{i})
                specSecElec = find(specifySecondaryAmps{i}(1,:) == electrodes(sInds(k)));
            else
                specSecElec = [];
            end
            if ~isempty(specSecElec)
                ampsToUse = specifySecondaryAmps{i}(2:end, specSecElec);
                ampsToUse = ampsToUse(ampsToUse ~= 0);
                for m = 1:length(ampsToUse)
                    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                    array(pInd, end) = centerAmps{i}(j);
                    array(sInds(k), end) = ampsToUse(m);
                    details.prePulse = [details.prePulse 0];
                    details.stimPulseAmps = [details.stimPulseAmps 0];
                    details.pairsCenterAmps = [details.pairsCenterAmps j];
                    details.clusterID = [details.clusterID i];
                end
            else
                for m = 1:length(amps.neighborAmps{i})
                    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                    array(pInd, end) = centerAmps{i}(j);
                    array(sInds(k), end) = amps.neighborAmps{i}(m);
                    details.prePulse = [details.prePulse 0];
                    details.stimPulseAmps = [details.stimPulseAmps 0];
                    details.pairsCenterAmps = [details.pairsCenterAmps j];
                    details.clusterID = [details.clusterID i];
                end
            end
        end
    end
end


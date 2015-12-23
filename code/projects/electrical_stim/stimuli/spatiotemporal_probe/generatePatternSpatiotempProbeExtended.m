function [electrodes array details] = generatePatternSpatiotempProbeExtended(centerElec, amps, varargin)
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
p.addRequired('amps', @isstruct)

p.addParamValue('specificSecondaryAmps', [], @isnumeric) %overrides secondary amplitudes specified in 'amps' on specified electrodes
% specify secondary amplitudes in this format: specificSecondaryAmps{clusterIndex} = [electrode1 electrode2 electrode3;
%                                                                                     amplitude  amplitude  amplitude;
%                                                                                     amplitude  amplitude  amplitude]


p.parse(centerElec, amps, varargin{:})

specifySecondaryAmps = p.Results.specificSecondaryAmps;



%%



% calculates center amplitudes
centerAmps = amps.centerAmpLims(1);
while centerAmps(end) < amps.centerAmpLims(2)
    centerAmps = [centerAmps centerAmps(end)*(1 + amps.centerIncrement)]; %#ok<AGROW>
end



% finds clusters of center + surrounding electrodes

[xCoords yCoords] = getElectrodeCoords61();


distances = zeros(64, 1);
angles = zeros(64, 1);
for j = 1:64
    delx = abs(xCoords(centerElec) - xCoords(j));
    dely = abs(yCoords(centerElec) - yCoords(j));
    
    distances(j) = norm([delx, dely]);
    angles(j) = (180/pi)*atan(dely/delx);
end

valAngles = abs(angles - 30) < 0.01 | abs(angles - 90) < 0.01;
valDists = abs(distances - 2) < 0.01 | abs(distances - 4) < 0.01 | abs(distances - 8) < 0.01;
cluster = find(valAngles & valDists);
electrodes = [centerElec; cluster]; %#ok<AGROW>

nElec = length(electrodes);


array = [];
details.prePulse = [];
details.stimPulseAmps = [];
details.pairsCenterAmps = [];


pInd = find(electrodes==centerElec);
sInds = find(electrodes~=centerElec);

%center alone, all amplitudes
for j = 1:length(centerAmps)
    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    array(pInd, end) = centerAmps(j);
    details.prePulse = [details.prePulse 0];
    details.stimPulseAmps = [details.stimPulseAmps j];
    details.pairsCenterAmps = [details.pairsCenterAmps 0];
end

%center alone, prepulse amplitudes
for j = 1:length(amps.centerAmpPrepulses)
    array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
    array(pInd, end) = amps.centerAmpPrepulses(j);
    details.prePulse = [details.prePulse 1];
    details.stimPulseAmps = [details.stimPulseAmps 0];
    details.pairsCenterAmps = [details.pairsCenterAmps 0];
end

%neighbor alone, each amplitude
for j = 1:length(sInds)
    if ~isempty(specifySecondaryAmps)
        specSecElec = find(specifySecondaryAmps(1,:) == electrodes(sInds(j)));
    else
        specSecElec = [];
    end
    if ~isempty(specSecElec)
        ampsToUse = specifySecondaryAmps(2:end, specSecElec);
        ampsToUse = ampsToUse(ampsToUse ~= 0);
        for k = 1:length(ampsToUse)
            array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
            array(sInds(j), end) = ampsToUse(k);
            details.prePulse = [details.prePulse 1];
            details.stimPulseAmps = [details.stimPulseAmps 0];
            details.pairsCenterAmps = [details.pairsCenterAmps 0];
        end
    else %use amplitudes specified in amps.neighborAmps
        for k = 1:length(amps.neighborAmps)
            array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
            array(sInds(j), end) = amps.neighborAmps(k);
            details.prePulse = [details.prePulse 1];
            details.stimPulseAmps = [details.stimPulseAmps 0];
            details.pairsCenterAmps = [details.pairsCenterAmps 0];
        end
    end
end

%center + neighbor
for j = 1:length(centerAmps)
    %jth center amplitude, all neighbors and all neighbor amplitudes
    for k = 1:length(sInds)
        if ~isempty(specifySecondaryAmps)
            specSecElec = find(specifySecondaryAmps(1,:) == electrodes(sInds(k)));
        else
            specSecElec = [];
        end
        if ~isempty(specSecElec)
            ampsToUse = specifySecondaryAmps(2:end, specSecElec);
            ampsToUse = ampsToUse(ampsToUse ~= 0);
            for m = 1:length(ampsToUse)
                array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                array(pInd, end) = centerAmps(j);
                array(sInds(k), end) = ampsToUse(m);
                details.prePulse = [details.prePulse 0];
                details.stimPulseAmps = [details.stimPulseAmps 0];
                details.pairsCenterAmps = [details.pairsCenterAmps j];
            end
        else
            for m = 1:length(amps.neighborAmps)
                array = [array zeros(nElec,1)];%#ok<AGROW> %add column of zeros to array
                array(pInd, end) = centerAmps(j);
                array(sInds(k), end) = amps.neighborAmps(m);
                details.prePulse = [details.prePulse 0];
                details.stimPulseAmps = [details.stimPulseAmps 0];
                details.pairsCenterAmps = [details.pairsCenterAmps j];
            end
        end
    end
end



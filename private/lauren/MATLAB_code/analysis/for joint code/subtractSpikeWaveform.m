function subtractedTrace = subtractSpikeWaveform(trace, template, offset)

%trace = matrix (electrodes x samples)
%template = matrix (electrodes x samples)
%trace and template must have the same size in the first dimension (same number of electrodes)

templateCell{1} = template;

if ~mod(offset, 1) %offset (latency) is a whole number
    subtractedTrace = subtractWithShifts(trace, templateCell, offset, offset, 1);
else
    fraction = mod(offset, 1);
    stepSize = min(fraction, 1-fraction);
    i = 2;
    while mod(stepSize*i, 1) > 0.0001 && i < 100
        i = i+1;
    end
    stepSize = 1/i;   
    offsetWindow = [floor(offset) ceil(offset)];
    offsetInd = find(offsetWindow(1):stepSize:offsetWindow(2) == offset) + 1;
    subVectors = subtractWithShifts(trace, templateCell, offsetWindow(1), offsetWindow(2), stepSize);
    subtractedTrace = subVectors{offsetInd};
end

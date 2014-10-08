function errors = calcErrors(traces, modelType, startFit, endFit, estArtifact, flags, prevOffsetInd)

% Calculates the 2-norm of residuals of traces (raw data after subtraction of spike templates
%    at various latencies) and model artifact.
% If model artifact is exponential or sum of exponentials, determines best-fit exponential and uses
%    it as the model artifact.
%
%
% arguments
%   traces:  a nTemplates-dimensional cell array, in which each cell is nElectrodes x nSamples
%   modelType: a string indicating which artifact model type was used (determines whether the norm
%      is calculated without first mean-zeroing, and whether an exponential model is fit to the
%      traces)
%   startFit, endFit: specifies range (in samples) over which norm will be calculated
%   estArtifact (for 'linkage', 'ttx', 'prevArtifact' and 'hypArtifact' model types only): the
%   estimated artifact (nElec x nSamples)
%   flags: a vector (length = nTemplates) indicating whether to force certain templates to have a spike
%     2 means force that particular template to have a spike (set the errors corresponding to that
%     neuron not spiking to infinity)
%     3 means force one of the templates to have a spike (set error corresponding to no templates
%     spiking to infinity)
%
% returns
%   errors: an array with the same dimensions as traces, with each value corresponding to the 2-norm
%   of the residual of each corresponding trace in traces and the model artifact
%
% modified 2009-06-06 to handle using multiple electrodes
% modified 2009-06-20 to handle hand-picked artifact model and to not mean-zero residuals for
% prevArtifact model
% modified 2009-06-23 to handle flags argument

% = 0 means calculated errors normally
% = 1 means force trace to be as defined in prevLatencies (if this occurs, prevLatencies must be
% supplied)
% = 2 means force trace to not be failure for specific template
% = 3 means force trace to not be a failure for one of the possible templates (should = 4 for all
% neurons)
% = 4 means force trace to be a failure for specific template

nShifts = size(traces, 1) - 1;
if any(size(traces) == 1)
    nTemplates = 1;
else
    nTemplates = max(size(size(traces)));
end

if ~exist('flags', 'var')
    flags = zeros(nTemplates, 1);
end

tempArray = traces{ones(nTemplates)}; %kluge to get "size" to work
nElecs = size(tempArray, 1);

errors = zeros(size(traces));

for j = 1:nElecs
    for i1 = 1:nShifts + 1
        if nTemplates < 2
            if (flags(1)==2 && i1==1) || (flags(1)==4 && i1~=1) || (flags(1)==1 && i1~=prevOffsetInd(1) + 1)
                errors(i1) = inf;
            elseif any(strcmpi(modelType, {'ttx', 'linkage'}))
                temp = traces{i1}(j, startFit:endFit) - estArtifact(j, startFit:endFit);
                errors(i1) = errors(i1) + norm(temp - mean(temp));
            elseif strcmpi(modelType, 'singExp')
                errors(i1) = errors(i1) + expFitter(traces{i1}(j, startFit:endFit));
            elseif strcmpi(modelType, 'sumExp')
                errors(i1) = errors(i1) + expFitter2(traces{i1}(j, startFit:endFit));
            elseif any(strcmpi(modelType, {'hypArtifact', 'prevArtifact', 'nextArtifact', 'handPicked', 'currentArtifact', 'otherOccurrences'}))
                errors(i1) = errors(i1) + norm(traces{i1}(j, startFit:endFit) - estArtifact(j, startFit:endFit));
            end
        else
            for i2 = 1:nShifts + 1
                if nTemplates < 3
                    if      (flags(1)==2 && i1==1) || (flags(1)==4 && i1~=1) || (flags(1)==1 && i1~=prevOffsetInd(1) + 1) ||...
                            (flags(2)==2 && i2==1) || (flags(2)==4 && i2~=1) || (flags(2)==1 && i2~=prevOffsetInd(2) + 1)
                        errors(i1,i2) = inf;
                    elseif any(strcmpi(modelType, {'ttx', 'linkage'}))
                        temp = traces{i1,i2}(j, startFit:endFit) - estArtifact(j, startFit:endFit);
                        errors(i1,i2) = errors(i1,i2) + norm(temp - mean(temp));
                    elseif strcmpi(modelType, 'singExp')
                        errors(i1,i2) = errors(i1,i2) + expFitter(traces{i1,i2}(j, startFit:endFit));
                    elseif strcmpi(modelType, 'sumExp')
                        errors(i1,i2) = errors(i1,i2) + expFitter2(traces{i1,i2}(j, startFit:endFit));
                    elseif any(strcmpi(modelType, {'hypArtifact', 'prevArtifact', 'nextArtifact', 'handPicked', 'currentArtifact', 'otherOccurrences'}))
                        errors(i1,i2) = errors(i1,i2) + norm(traces{i1,i2}(j, startFit:endFit) - estArtifact(j, startFit:endFit));
                    end
                else
                    for i3 = 1:nShifts + 1
                        if nTemplates < 4
                            if      (flags(1)==2 && i1==1) || (flags(1)==4 && i1~=1) || (flags(1)==1 && i1~=prevOffsetInd(1) + 1) ||...
                                    (flags(2)==2 && i2==1) || (flags(2)==4 && i2~=1) || (flags(2)==1 && i2~=prevOffsetInd(2) + 1) ||...
                                    (flags(3)==2 && i3==1) || (flags(3)==4 && i3~=1) || (flags(3)==1 && i3~=prevOffsetInd(3) + 1)
                                errors(i1,i2,i3) = inf;
                            elseif any(strcmpi(modelType, {'ttx', 'linkage'}))
                                temp = traces{i1,i2,i3}(j, startFit:endFit) - estArtifact(j, startFit:endFit);
                                errors(i1,i2,i3) = errors(i1,i2,i3) + norm(temp - mean(temp));
                            elseif strcmpi(modelType, 'singExp')
                                errors(i1,i2,i3) = errors(i1,i2,i3) + expFitter(traces{i1,i2,i3}(j, startFit:endFit));
                            elseif strcmpi(modelType, 'sumExp')
                                errors(i1,i2,i3) = errors(i1,i2,i3) + expFitter2(traces{i1,i2,i3}(j, startFit:endFit));
                            elseif any(strcmpi(modelType, {'hypArtifact', 'prevArtifact', 'nextArtifact', 'handPicked', 'currentArtifact', 'otherOccurrences'}))
                                errors(i1,i2,i3) = errors(i1,i2,i3) + norm(traces{i1,i2,i3}(j, startFit:endFit) - estArtifact(j, startFit:endFit));
                            end
                        else
                            for i4 = 1:nShifts + 1
                                if nTemplates < 5
                                    if      (flags(1) == 2 && i1 == 1) || (flags(1) == 4 && i1 ~= 1) || (flags(1)==1 && i1~=prevOffsetInd(1) + 1) ||...
                                            (flags(2) == 2 && i2 == 1) || (flags(2) == 4 && i2 ~= 1) || (flags(2)==1 && i2~=prevOffsetInd(2) + 1) ||...
                                            (flags(3) == 2 && i3 == 1) || (flags(3) == 4 && i3 ~= 1) || (flags(3)==1 && i3~=prevOffsetInd(3) + 1) ||...
                                            (flags(4) == 2 && i4 == 1) || (flags(4) == 4 && i4 ~= 1) || (flags(4)==1 && i4~=prevOffsetInd(4) + 1)
                                        errors(i1,i2,i3,i4) = inf;
                                    elseif any(strcmpi(modelType, {'ttx', 'linkage'}))
                                        temp = traces{i1,i2,i3,i4}(j, startFit:endFit) - estArtifact(j, startFit:endFit);
                                        errors(i1,i2,i3,i4) = errors(i1,i2,i3,i4) + norm(temp - mean(temp));
                                    elseif strcmpi(modelType, 'singExp')
                                        errors(i1,i2,i3,i4) = errors(i1,i2,i3,i4) + expFitter(traces{i1,i2,i3,i4}(j, startFit:endFit));
                                    elseif strcmpi(modelType, 'sumExp')
                                        errors(i1,i2,i3,i4) = errors(i1,i2,i3,i4) + expFitter2(traces{i1,i2,i3,i4}(j, startFit:endFit));
                                    elseif any(strcmpi(modelType, {'hypArtifact', 'prevArtifact', 'nextArtifact', 'handPicked', 'currentArtifact', 'otherOccurrences'}))
                                        errors(i1,i2,i3,i4) = errors(i1,i2,i3,i4) + norm(traces{i1,i2,i3,i4}(j, startFit:endFit) - estArtifact(j, startFit:endFit));
                                    end
                                else
                                    error('Current state of algorithm can''t handle more than 4 templates')
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if flags(1) == 3
    switch nTemplates
        case 1
            errors(1) = inf;
        case 2
            errors(1,1) = inf;
        case 3
            errors(1,1,1) = inf;
        case 4
            errors(1,1,1,1) = inf;
    end
end

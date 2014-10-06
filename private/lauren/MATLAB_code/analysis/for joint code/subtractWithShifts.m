function [subtractedTraces offsets] = subtractWithShifts(trace, templates, shiftStart, shiftEnd, shiftStep)

% subtracts a template waveform from a data trace at various offsets
%
% arguments: 
%   trace (vector/2-dim. array) electrodes x samples
%   templates (cell array of vectors/2-dim. arrays: electrodes x samples) - the waveform to be 
%   subtracted at different offsets from 0
%   shiftStart (integar) - position of the first offset, where 0 indicates alignment of first sample of trace
%   and template, negative indicates that the template is shifted to before the trace (one or more
%   of the initial samples of the template aren't used) and positive indicates that the template is
%   shifted to after the start of the trace
%   shiftEnd (integar) - analagous to shiftStart, must be >= to shiftStart
%   shiftStep (integar or fraction of one) - step size between subsequent shifts (shiftEnd - 
%   shiftStart must be divisible by shiftStep)
%
%   ***to find subtracted trace for single offset, use only 3 arguments: (trace, templates,
%   offset) (only works if size(templates) == (1,1)
%
%
% returns:
%   subtractedTraces - a cell array of vectors representing the trace after subtraction of the
%   template at the various offsets, with dimensions (template 1 subtraction (different offsets) x template 2
%   subtraction x ... x template n subtraction)
%   offsets - array of offset values corresponding to the vectors in subtractedTraces
%   
%   ***when using single offset, subtractedTraces is an array (not a cell array)
%   
% author: Lauren Hruby, Jan 2009
% 
% updates: allowed for multiple template subtraction (up to 4)
% allowed for multiple electrodes (multiple traces per template and trace)
%     * number of electrodes must be less than the number of samples in template or trace
%
% July 2009: made option to find subtracted trace for single offset (only works for single template)
%




% makes sure first dimension of traces and templates corresponds to electrode
if size(trace,1) > size(trace,2)
    trace = trace';
end

if size(templates{1},1) > size(templates{1},2)
    for i = 1:length(templates)
        templates{i} = templates{i}';
    end
end


% generates shiftStart, shiftEnd and shiftStep values for single offset case
if nargin < 4 %calculate for single offset value
    if max(size(templates)) > 1
        warndlg('can''t use single offset option with multiple templates')
        uiwait(gcf)
        return
    end
    offset = shiftStart;
    %templateCell{1} = templates;
    
    if ~mod(offset, 1) %offset (latency) is a whole number
        shiftStart = offset;
        shiftEnd = offset;
        shiftStep = 1;
        offsetInd = 2;
    else
        fraction = mod(offset, 1);
        stepSize = min(fraction, 1-fraction);
        i = 2;
        while mod(stepSize*i, 1) > 0.0001 && i < 100
            i = i+1;
        end
        if i == 100
            warnh = warndlg('suitable step size could not be found');
            uiwait(warnh)
            return
        end
            
        shiftStep = 1/i;
        shiftStart = floor(offset);
        shiftEnd = ceil(offset);
        offsetInd = find(shiftStart:shiftStep:shiftEnd == offset) + 1;
    end
end


lTrace = length(trace);
lTemplate = length(templates{1});
nShifts = (shiftEnd - shiftStart)/shiftStep + 1;
nTemplate = length(templates);
nElecs = size(trace, 1);

if size(templates{1},1) ~= nElecs
    error('number of electrodes in template and trace are not consistent')
end

if nShifts ~= round(nShifts)
    error('shiftEnd - shiftStart must be divisible by shiftStep')
end

if shiftStep < 1 % requires interpolation of templates
    shiftDenom = 1/shiftStep;
    if shiftDenom - round(shiftDenom) < 0.0001
        shiftDenom = round(shiftDenom);
        %piecewise cubic hermite interpolation
        tempInterp = cell(nTemplate, shiftDenom);
        for i = 1:nTemplate
            for j = 1:shiftDenom %i = 1 corresonds with integar shift
                %tempInterp{i,j} = pchip(1:lTemplate, templates{i}, (1 : lTemplate) + (j-1)/shiftDenom);
                tempInterp{i,j} = zeros(nElecs, lTemplate);
                for k = 1:nElecs
                    tempInterp{i,j}(k,:) = pchip(1:lTemplate, templates{i}(k,:), (1 : lTemplate) + (j-1)/shiftDenom);
                end
            end
        end
    else
        error('shiftStep must be a fraction of 1 (e.g. 1/2, 1/4)')
    end
end



%% subtraction

offsets = shiftStart:shiftStep:shiftEnd;
offsets = [nan offsets];
subtractionVector = cell(nTemplate, nShifts);

% if shiftStep is a fraction
if exist('shiftDenom','var')
    
    for i = 1:nTemplate
        %% first integer offset
        subtractionVector{i,1} = zeros(nElecs, lTrace);
        if shiftStart >= 0
            if shiftStart + lTemplate <= lTrace
                subtractionVector{i,1}(:, shiftStart + 1 : shiftStart + lTemplate) = tempInterp{i,1};
            else
                subtractionVector{i,1}(:, shiftStart + 1 : lTrace) = tempInterp{i,1}(:, 1 : lTrace - shiftStart);
            end
        else
            if shiftStart + lTemplate <= lTrace
                subtractionVector{i,1}(:, 1 : shiftStart + lTemplate) = tempInterp{i,1}(:, -shiftStart + 1 : lTemplate);
            else
                subtractionVector{i,1} = tempInterp{i,1}(:, -shiftStart + 1 : lTrace - shiftStart);
            end
        end
        
        %% rest of the offsets
        for j = 1:(shiftEnd - shiftStart)
            currentShiftInt = shiftStart + j;

            for k = 1:shiftDenom % k = 1 corresponds to integar shift = currentShiftInt+1
                subtractionVector{i,(j)*shiftDenom - k + 2} = zeros(nElecs, lTrace);
                if currentShiftInt >= 0
                    if currentShiftInt + lTemplate <= lTrace
                        subtractionVector{i,(j)*shiftDenom - k + 2}(:, currentShiftInt + 1 : currentShiftInt + lTemplate) = tempInterp{i,k};
                    else
                        subtractionVector{i,(j)*shiftDenom - k + 2}(:, currentShiftInt + 1 : lTrace) = tempInterp{i,k}(:, 1 : lTrace - currentShiftInt);
                    end
                else
                    if currentShiftInt + lTemplate <= lTrace
                        subtractionVector{i,(j)*shiftDenom - k + 2}(:, 1 : currentShiftInt + lTemplate) = tempInterp{i,k}(:, -currentShiftInt + 1 : lTemplate);
                    else
                        subtractionVector{i,(j)*shiftDenom - k + 2} = tempInterp{i,k}(:, -currentShiftInt + 1 : lTrace - currentShiftInt);
                    end
                end
            end
        end
    end
    
else %shiftStep is an integer
    for i = 1:nTemplate
        for j = 1:nShifts
            subtractionVector{i,j} = zeros(nElecs, lTrace);
            currentShift = offsets(j+1);
            if currentShift >= 0 %no need to truncate beginning of template
                if currentShift + lTemplate <= lTrace %no need to truncate end of template
                    subtractionVector{i,j}(:, currentShift + 1 : currentShift + lTemplate) = templates{i};
                else
                    subtractionVector{i,j}(:, currentShift + 1 : lTrace) = templates{i}(:, 1 : lTrace - currentShift);
                end
            else
                if currentShift + lTemplate <= lTrace
                    subtractionVector{i,j}(:, 1 : currentShift + lTemplate) = templates{i}(:, -currentShift + 1 : lTemplate);
                else
                    subtractionVector{i,j} = templates{i}(:, -currentShift + 1 : lTrace - currentShift);
                end
            end
        end
    end
end


subDimensions = ones(nTemplate,1)'*(nShifts+1);
if nTemplate == 1
    subtractedTraces = cell(nShifts+1, 1);
else
    subtractedTraces = cell(subDimensions);
end

for i1 = 0:nShifts
    if i1 == 0     
        sumSubVector1 = zeros(nElecs, lTrace);
    else
        sumSubVector1 = subtractionVector{1,i1};
    end
    if nTemplate < 2
        subtractedTraces{i1+1} = trace - sumSubVector1;
    else
        for i2 = 0:nShifts
            if i2 == 0
                sumSubVector2 = sumSubVector1;
            else
                sumSubVector2 = subtractionVector{2,i2} + sumSubVector1;
            end
            if nTemplate < 3
                subtractedTraces{i1+1, i2+1} = trace - sumSubVector2;
            else
            
                for i3 = 0:nShifts
                    if i3 == 0
                        sumSubVector3 = sumSubVector2;
                    else
                        sumSubVector3 = subtractionVector{3,i3} + sumSubVector2;
                    end
                    if nTemplate < 4
                        subtractedTraces{i1+1, i2+1, i3+1} = trace - sumSubVector3;
                    else
                        
                        
                        for i4 = 0:nShifts
                            
                            if i4 == 0
                                sumSubVector4 = sumSubVector3;
                            else
                                sumSubVector4 = subtractionVector{4,i4} + sumSubVector3;
                            end
                            subtractedTraces{i1+1, i2+1, i3+1, i4+1} = trace - sumSubVector4;
                        end
                    end
                end             
            end
        end
    end
end

if nargin < 4
    subtractedTraces = subtractedTraces{offsetInd};
    offsets = offsets(offsetInd);
end










function spikeMinWindow = convertTempOffsetWindow(elecResp, movieInd)
%
% elecResp.analysis.details.tempOffsetWindow is specified in terms of
% offset between template (ei trace) and stim trace
%
% this function converts these offset limits to limits in terms of the
% minimum and maximum sample at which the spike minimum is allowed to occur
%
%
%

centInd = find(elecResp.cells.goodElecs == elecResp.cells.recElec);

neuronID = [elecResp.cells.main elecResp.cells.active{movieInd}];
nTemplates = length(neuronID);

templates = cell(nTemplates, 1);
templateMinPos = zeros(nTemplates, 1);


for i = 1:nTemplates
    templates{i} = elecResp.cells.allEIs{elecResp.cells.all == neuronID(i)}(elecResp.cells.goodElecs, :);
    templateMinPos(i) = find(squeeze(templates{i}(centInd,:)) == min(squeeze(templates{i}(centInd,:))));
end

tempOffsetWindow = elecResp.analysis.details.tempOffsetWindow{movieInd};

spikeMinWindow(1) = tempOffsetWindow(1) + max(templateMinPos);
spikeMinWindow(2) = tempOffsetWindow(2) + min(templateMinPos);

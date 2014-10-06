function [cellInfo cellInfoOffArray] = removeSmallSigEdgeCells(cellInfo, pathToEi, varargin)

% separates out cells that are not over 61-elec array based whether they meet one
% of 2 criteria:
% The peak negative signal must EITHER
% (1) not fall on an edge electrode OR
% (2) be >= the threshold value
%
% the threshold is either:
%    the mean - 1 SD of the peak (negative) signal of
%    cells within cell class that meet (1)  (used if cutoffMult parameter is
%    specified as NaN)
% OR
%    cutoffMult*the mean peak (negative) signal of
%    cells within cell class that meet (1) (used if cutoffMult is
%    specified, if not specified it is set to 0.5)
%
% returns cells that don't meet one of the criteria in cellInfoOffArray and
% removes them from cellInfo
%
% also adds/populates cellInfo fields: .maxSigAllElecs, .maxSigRecElec



p = inputParser;

p.addRequired('cellInfo', @isstruct)
p.addRequired('pathToEi', @ischar)

p.addParamValue('cutoffMult', 0.5, @isnumeric)


p.parse(cellInfo, pathToEi, varargin{:})

cutoffMult = p.Results.cutoffMult;


edgeElecs = [2 5 8 10 13 16 19 20 23 27 30 31 34 37 40 42 45 48 51 52 55 59 62 63];



% get ei info and collect max signal values for non-edge cells in each class

onM_NEPS = [];
offM_NEPS = [];
onP_NEPS = [];
offP_NEPS = [];
sbc_NEPS = [];

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
for ii = 1:length(cellInfo)
    if isfield(cellInfo, 'id')
        cellID = cellInfo(ii).id; %for threshold summary script
    else
        cellID = cellInfo(ii).targetCell; %for selectivity summary script
    end
    ei = eiFile.getImage(cellID);
    %calculate negative peak (absolute) signal on each electrode
    eiAmps = zeros(64, 1);
    for j = 1:64
        if ~(j==9||j==25||j==57)
            eiAmps(j) = max(max(-ei(1,j+1,:)));
        end
    end
    cellInfo(ii).maxSigAllElecs = max(eiAmps);
    cellInfo(ii).maxSigRecElec = find(eiAmps == max(eiAmps));
    
    if length(cellInfo(ii).maxSigRecElec) ~= 1 %hopefully this won't happen
        error('more than one electrode has same negative peak signal value')
    end
    if ~any(edgeElecs == cellInfo(ii).maxSigRecElec) %peak sig not on edge elec
        if     strcmpi(cellInfo(ii).type, 'onMidg')
            onM_NEPS  = [onM_NEPS  cellInfo(ii).maxSigAllElecs]; %#ok<*AGROW>
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            offM_NEPS = [offM_NEPS cellInfo(ii).maxSigAllElecs];
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            onP_NEPS  = [onP_NEPS  cellInfo(ii).maxSigAllElecs];
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            offP_NEPS = [offP_NEPS cellInfo(ii).maxSigAllElecs];
        elseif strcmpi(cellInfo(ii).type, 'sbc')
            sbc_NEPS = [sbc_NEPS cellInfo(ii).maxSigAllElecs];
        else
            error(['invalid cell type for neuron ' num2str(cellID)])
        end
    end
end
eiFile.close() %put your tird in the toilet!

% if ~isnan(cutoffMult)
%     disp(['standard deviation of ON midget (N= ' num2str(length(onM_NEPS)) ') non-edge peak signals = ' num2str(std(onM_NEPS)) ' and ' num2str(cutoffMult) '*mean = ' num2str(cutoffMult*mean(onM_NEPS))])
%     disp(['standard deviation of OFF midget (N= ' num2str(length(offM_NEPS)) ') non-edge peak signals = ' num2str(std(offM_NEPS)) ' and ' num2str(cutoffMult) '*mean = ' num2str(cutoffMult*mean(offM_NEPS))])
%     disp(['standard deviation of ON parasol (N= ' num2str(length(onP_NEPS)) ') non-edge peak signals = ' num2str(std(onP_NEPS)) ' and ' num2str(cutoffMult) '*mean = ' num2str(cutoffMult*mean(onP_NEPS))])
%     disp(['standard deviation of OFF parasol (N= ' num2str(length(offP_NEPS)) ') non-edge peak signals = ' num2str(std(offP_NEPS)) ' and ' num2str(cutoffMult) '*mean = ' num2str(cutoffMult*mean(offP_NEPS))])
%     disp(['standard deviation of sbc (N= ' num2str(length(sbc_NEPS)) ') non-edge peak signals = ' num2str(std(sbc_NEPS)) ' and ' num2str(cutoffMult) '*mean = ' num2str(cutoffMult*mean(sbc_NEPS))])
% end

%remove cells that don't meet one of the criteria
count = 0;
for ii = length(cellInfo):-1:1
    if ~isnan(cutoffMult)
        if     strcmpi(cellInfo(ii).type, 'onMidg')
            cutoff = cutoffMult*mean(onM_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            cutoff = cutoffMult*mean(offM_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            cutoff = cutoffMult*mean(onP_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            cutoff = cutoffMult*mean(offP_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'sbc')
            cutoff = cutoffMult*mean(sbc_NEPS);
        else
            error(['invalid cell type for neuron ' num2str(cellID)])
        end
    else
        if     strcmpi(cellInfo(ii).type, 'onMidg')
            cutoff = mean(onM_NEPS)-std(onM_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'offMidg')
            cutoff = mean(offM_NEPS)-std(offM_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'onPar')
            cutoff = mean(onP_NEPS)-std(onP_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'offPar')
            cutoff = mean(offP_NEPS)-std(offP_NEPS);
        elseif strcmpi(cellInfo(ii).type, 'sbc')
            cutoff = mean(sbc_NEPS)-std(sbc_NEPS);
        else
            error(['invalid cell type for neuron ' num2str(cellID)])
        end
    end
    
    if any(edgeElecs == cellInfo(ii).maxSigRecElec) && cellInfo(ii).maxSigAllElecs < cutoff
        count = count + 1;
        cellInfoOffArray(count) = cellInfo(ii);
        cellInfo(ii) = [];
    end
end

if ~exist('cellInfoOffArray', 'var')
    cellInfoOffArray = [];
end

end


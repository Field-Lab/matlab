function [templates recElecs]=makeTemplatesFromElecResp(pathToAnalysisData,patternNo,neuronIds,varargin)

if(nargin>3)
    recElecs=varargin{1};
elseif(nargin==3)
    recElecs=[];
    for n=1:length(neuronIds)
        str=[pathToAnalysisData 'elecResp_n' num2str(neuronIds(n)) '_p' num2str(patternNo) '.mat'];
        load(str)
        recElecs=[recElecs elecResp.cells.goodElecs];
    end
end


for n=1:length(neuronIds)
    str=[pathToAnalysisData 'elecResp_n' num2str(neuronIds(n)) '_p' num2str(patternNo) '.mat'];
    load(str);
    templates{n}=elecResp.cells.mainEI(recElecs,:);
end


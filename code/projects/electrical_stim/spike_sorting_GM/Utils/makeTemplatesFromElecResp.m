function [templates recElecs]=makeTemplatesFromElecResp(pathToAnalysisData,patternNo,neuronIds,Tmin,varargin)
% makeTemplatesFromElecResp Creates templates from elecResp files
% loadData loads data from a movies adds manually an Axonal Bundle Breakpoint to the input structure, defined as the
% last condition j for which there is no observed traveling wave following stimulation.
% 
% inputs:   -pathToAnalysisData: path to movie and elecResp files
%           -patternNo: pattern for spike sorting 
%           -neuronIds: Vector of Ids of neurons that will be part of the analysis         
%           -Tmin: Minimum time of templates that will be considered. Choose one if original templates
%            are aligned to onset of spike at time ~10. (this should be regularly the case for templates of length either 71 or 111)
%            choose Tmin>10 to cut the first part of the template, if there is misalignment with above rule (That seems to be the case when
%            templates have length=81, Tmin=10 in that case
% optional: -recElecs: vector of recording electrodes. If they are not provided, they will be chosen from elecResp.cells.goodElecs file
%          
% Output:   templates: a nNeuron dimensional cell array, such that templates{n} is a E*T_template matrix 
%           recElecs: E dimensional vector of electrodes recElecs (same as input, if they are provided)  
% Gonzalo Mena 6/2015 



if(nargin>4)
    recElecs=varargin{1};
elseif(nargin==4)
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
    templates{n}=elecResp.cells.mainEI(recElecs,Tmin:end);
end


function [templates, recElecs] = makeTemplatesFromEi(pathToEi, neuronIds,  varargin)
% makeTemplatesFromEi creates cell templates from .ei files
% loadData loads data from a movies adds manually an Axonal Bundle Breakpoint to the input structure, defined as the
% last condition j for which there is no observed traveling wave following stimulation.
% 
% inputs:   -pathToEi: path/to/my/data[kkk].ei 
%           -neuronIds: Vector of Ids of neurons that will be part of the analysis  
%           
% optional: -recElecs: vector of recording electrodes. If they are not
% provided, one will be chosen based on the highest recorded amplitude in
% the ei
%          
% Output:   templates: a nNeuron dimensional cell array, such that templates{n} is a E*T_template matrix 
%           recElecs: E dimensional vector of electrodes recElecs (same as input, if they are provided)  
% Lauren Grosberg 7/2015


if(nargin>2)
    recElecs=varargin{1};
elseif(nargin==2)
    recElecs=[];
end

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(pathToEi);
for n=1:length(neuronIds)
    neuronEI = eiFile.getImage(neuronIds(n));
    neuronEI_volt = squeeze(neuronEI(1,2:end,:)); 
    if isempty(recElecs)
        recElecs = find(min(neuronEI_volt,[],2) == min(neuronEI_volt(:)));
    end
    templates{n}=neuronEI_volt(recElecs,:);
end

end %function end
function prefElectrodes = PreferredElectrodes(templates,varargin)
% PreferredElectrodes() defines the preferred or preferreds electrodes for each neuron.
% 
% inputs:   path: cell array templates 
% optional: nNeurons dimensional vector with the number of preferred
% electrodes to be chosen for each neuron. In that case, for each neuron
% the electrodes with the strongest signal will be chosen.
%          
% example: prefElectrodes = PreferredElectrodes(templates,[1 2 1])
% will chose 1 preferred electrode for neuron 1, 2 for neuron 2 and 1 for
% neuron 3
% output:  prefElectrodes{n}(e)
%
% Gonzalo Mena 6/2015 



if(nargin==1)
    nPrefElectrode=ones(1,length(templates));
else
    nPrefElectrode=varargin{1};
end
for n=1:length(templates)
    [StrengthSorted indStrength]  = sort(max(abs(templates{n})'),'descend');
    prefElectrodes{n}             = indStrength(1:nPrefElectrode(n));
end
function prefElectrodes = PreferredElectrodes(templates,varargin)

if(nargin==1)
    nPrefElectrode=ones(1,length(templates));
else
    nPrefElectrode=varargin{1};
end
for n=1:length(templates)
    [StrengthSorted indStrength]  = sort(max(abs(templates{n})'),'descend');
    prefElectrodes{n}             = indStrength(1:nPrefElectrode(n));
end
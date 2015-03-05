function load_forel

global thr fileName pathName rgc spikes getXdata getYdata visibleUnits signUnit fatPoints ifUpdateISI

sep=filesep;
a=find(pathName==sep,2,'last');
forelPath=[pathName(1:a(1)),'forel',sep];
load([forelPath,fileName])
set_thr(3)

spikesSubset=min(spikes)<thr;
tmp=find(spikesSubset);
for i=1:min(length(clusters),9)
    clusters{i}=tmp(clusters{i})';    
end
rgc=clusters;
visibleUnits=ones(1,size(rgc,2));
signUnit=zeros(1,size(rgc,2));
fatPoints=zeros(1,size(rgc,2));
ifUpdateISI=1;
redraw(get(getXdata,'value'),get(getYdata,'value'))

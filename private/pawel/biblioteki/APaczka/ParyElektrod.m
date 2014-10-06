clc
clear all
close all

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1); 
%define the electrode map - must be different than 1 for silicon probes
Radius=1;
CenterChannel = 58;
ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)'

% for i=1:64
%     CenterChannel=i;
%     ChannelsPlot=electrodeMap.getAdjacentsTo(CenterChannel,Radius)'
% end

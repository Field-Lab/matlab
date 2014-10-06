function [ ] = PlotCoordinates( ElectrodeOnFig )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);     %define the electrode map - must be different than 1 for silicon probes

figure (1)
                
   for i=1:ElectrodeOnFig
      CoordinatesX=electrodeMap.getXPosition(i);
      CoordinatesY=electrodeMap.getYPosition(i);
      plot(CoordinatesX,CoordinatesY,'.');
 	hold on
      text(CoordinatesX,CoordinatesY,num2str(i))
      hold on
   end
      axis off
      hold on

end


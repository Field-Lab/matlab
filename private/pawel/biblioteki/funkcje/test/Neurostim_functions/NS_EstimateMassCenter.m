function [x,y]=NS_EstimateMassCenter(Amplitudes,electrodes);

ArrayID=1; %61-channel
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
for i=1:length(electrodes)
    CoordinatesX(i)=electrodeMap.getXPosition(electrodes(i));
    CoordinatesY(i)=electrodeMap.getYPosition(electrodes(i));
end

A=sum(Amplitudes);
X=sum(Amplitudes.*CoordinatesX);
Y=sum(Amplitudes.*CoordinatesY);

x=X/A;
y=Y/A;
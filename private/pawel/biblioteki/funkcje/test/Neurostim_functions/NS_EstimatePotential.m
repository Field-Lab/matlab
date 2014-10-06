function [V]=NS_EstimatePotential(Amplitudes,Electrodes,x,y,z);

ArrayID=1; %61-channel
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
for i=1:length(Electrodes)
    elX(i)=electrodeMap.getXPosition(Electrodes(i));
    elY(i)=electrodeMap.getYPosition(Electrodes(i));
    elZ(i)=0; %for planar array :)
end
Amplitudes
V=0;
for i=1:length(Electrodes)
    V=V+Amplitudes(i)/sqrt((elX(i)-x)^2+(elY(i)-y)^2+(elZ(i)-z)^2);
end
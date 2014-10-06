function [GX,GY,GZ]=NS_EstimatePotentialGradientForPattern(Amplitudes,Electrodes,x0,y0,z0);

ArrayID=1; %61-channel
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
for i=1:length(Electrodes)
    elX(i)=electrodeMap.getXPosition(Electrodes(i));
    elY(i)=electrodeMap.getYPosition(Electrodes(i));
    elZ(i)=0; %for planar array :)
end

[GX,GY,GZ]=NS_EstimatePotentialGradient(Amplitudes,elX,elY,elZ,x0,y0,z0);
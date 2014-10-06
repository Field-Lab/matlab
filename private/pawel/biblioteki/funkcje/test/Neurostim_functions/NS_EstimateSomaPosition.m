function [x,y]=NS_EstimateSomaPosition(EI,electrodes);

ArrayID=1; %61-channel
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
for i=1:length(electrodes)
    CoordinatesX(i)=electrodeMap.getXPosition(electrodes(i));
    CoordinatesY(i)=electrodeMap.getYPosition(electrodes(i));
end

Amplitude=zeros(1,length(electrodes));
for i=1:length(electrodes)
    Amplitude(i)=abs(min(EI(i,:)));
end
A=sum(Amplitude);
X=sum(Amplitude.*CoordinatesX);
Y=sum(Amplitude.*CoordinatesY);

x=X/A;
y=Y/A;
Amplitude;
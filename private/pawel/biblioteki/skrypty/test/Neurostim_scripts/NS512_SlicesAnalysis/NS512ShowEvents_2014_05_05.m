electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
% Jakosc eventu:

fid = fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001_Matlab\figures3\NeuronsPatterns2.bin'],'r');
Data0=fread(fid,'int32');
fclose(fid);
Data1=reshape(Data0,3,length(Data0)/3);

EventsToShow=find(Data1(3,:)==1)

% primary recording electrode:
fid = fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\NeuronsPatternsPrimaryElectrodes.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data2=reshape(Data0,3,length(Data0)/3);
%EventData=Data(:,EventNumber);
figure(5)
clf
hold on

for i=1:length(EventsToShow)
    i
    Event=EventsToShow(i);
    StimEl=Data1(2,Event);
    RecEl=Data2(3,Event);
    
    X1=electrodeMap.getXPosition(StimEl);
    Y1=electrodeMap.getYPosition(StimEl);
    
    X2=electrodeMap.getXPosition(RecEl);
    Y2=electrodeMap.getYPosition(RecEl);
    
    dx=NS512_ConnectivityLine([X1 X2],[Y1 Y2],[0 0 1],[1 0 0],100,[1:5])%[1:5 96:100]);
end

PrimaryElectrodesToShow=unique(Data2(3,EventsToShow))
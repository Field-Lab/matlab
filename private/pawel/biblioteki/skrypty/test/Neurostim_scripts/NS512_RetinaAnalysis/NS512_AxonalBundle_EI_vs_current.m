clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
ArrayID=500;
DataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
Movies=[1:2:63];

ThresholdFilePath1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_1']
fid1=fopen(ThresholdFilePath1,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons1=NeuronInformation(:,1);
Electrodes1=NeuronInformation(:,3);
Amplitudes1=NeuronInformation(:,4);

ThresholdFilePath2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_2']
fid1=fopen(ThresholdFilePath2,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons2=NeuronInformation(:,1);
Electrodes2=NeuronInformation(:,3);
Amplitudes2=NeuronInformation(:,4);

Neurons=[Neurons1' Neurons2']';
Electrodes=[Electrodes1' Electrodes2']';
Amplitudes=[Amplitudes1' Amplitudes2']';

Patterns=sort(Electrodes);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
Channels=[1:512];
Coordinates=zeros(2,length(Channels));
ChannelsNumber=length(Channels);
for i=1:ChannelsNumber
    Coordinates(1,i)=electrodeMap.getXPosition(Channels(i));
    Coordinates(2,i)=electrodeMap.getYPosition(Channels(i));
end

N=length(Patterns)

% znajdz dla kazdej stymulujacej elektrody te elektrody odczytowe, ktora
% pokauja dla najwiekszego pradu sygnal (EI) wiekszy niz 20
for i=[25 58 64]%[1:6 8:11 13:22 24:26 28:32 34:54 56:70]
    OutputFileName1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby5\p' num2str(i) '.gif'];
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Patterns(i),63,0,0);
    DataTraces=DataTraces0(1:50,[1:512],[8:137]);
    d=reshape(mean(DataTraces),512,130);
    
    P=max(d')-min(d'); 
    P(Patterns(i))=0;
    ElectrodesWithSignal=find(P>100);
    length(ElectrodesWithSignal)
end
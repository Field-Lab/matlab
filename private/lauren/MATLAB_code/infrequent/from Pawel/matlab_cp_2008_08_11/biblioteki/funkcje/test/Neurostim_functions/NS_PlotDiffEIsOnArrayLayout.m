function [DiffEI]=NS_PlotDiffEIsOnArrayLayout(StimulatingChannel,MovieNumber1,MovieNumber2,RecElectrodes,BadChannels,ReadPath,FileName,WritePath);
%Calculates PCA and writes two figures to hard drive: one showing two
%dimensions of PCA space, and one showing 7 electrodes with clustered
%waveform - all of that for given channel and movie number.

% 1. Global constants
ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

% 2. File, channel, movie
%cd E:\2008-03-21-0;
cd(ReadPath);
%FileName='003';

% 3. Analysis parameters (for PCA)
TimeStart=0; %-10;
NumberOfSamples=60;
SamplesForPCA=[9:38];
PCA_channels=3; %1 - only stimulating channels; 2 - only neighbors; 3 - both stimulating one and neighbors
Dimensions=3;
NumberOfClusters=4;

% 4. 

% 5. Build data for PCA
[Timings1,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber1,NS_GlobalConstants);
[Timings2,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber2,NS_GlobalConstants);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
neighbors0 = electrodeMap.getAdjacentsTo(StimulatingChannel,1);
%neighbors0=[4 6 7 10:22];
neighbors=[];
for i=1:length(neighbors0)
    active=1;
    for j=1:length(BadChannels)
        if neighbors0(i)==BadChannels(j)
            active=0;
        end
    end
    if active==1
        neighbors=[neighbors neighbors0(i)];
    end
end
if RecElectrodes(1)~=0
    neighbors=RecElectrodes;
end
%neighbors=[16 13];
Offsets=ones(length(neighbors))*(-370);
Traces1=NS_ReadManyTracesFromRaw(FileName,neighbors,Timings1,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
Traces2=NS_ReadManyTracesFromRaw(FileName,neighbors,Timings2,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
DiffEI=mean(Traces2)-mean(Traces1);
Traces=Traces2-Traces1;
%Traces0=Traces([1:73 75:80 82:98],:,:);
%Traces=Traces0;
STraces=size(Traces);
b=mean(Traces);
TracesOrig=Traces;
for i=1:STraces(1)
Traces(i,:,:)=Traces(i,:,:)-b(1,:,:);
end

switch PCA_channels
    case 1,
        dane1=Traces1(:,1,SamplesForPCA);
    case 2,
        dane1=Traces1(:,2:length(neighbors),SamplesForPCA);      
    case 3,
        dane1=Traces1(:,:,SamplesForPCA);
end
Sdane1=size(dane1)
dane11=NS_ConcacenateWaveforms(dane1);

[Types1,PCA_Coeffs1,Inc1]=NS_ClusterSignatures(dane11,Dimensions,NumberOfClusters);
TypesNew1=NS_CleanClustersNew(PCA_Coeffs1,Types1,5);
a=find(TypesNew1~=0);
EI1=mean(Traces1(a,:,:));

switch PCA_channels
    case 1,
        dane2=Traces2(:,1,SamplesForPCA);
    case 2,
        dane2=Traces2(:,2:length(neighbors),SamplesForPCA);      
    case 3,
        dane2=Traces2(:,:,SamplesForPCA);
end
Sdane2=size(dane2)
dane21=NS_ConcacenateWaveforms(dane2);

[Types2,PCA_Coeffs2,Inc2]=NS_ClusterSignatures(dane21,Dimensions,NumberOfClusters);
TypesNew2=NS_CleanClustersNew(PCA_Coeffs2,Types2,5);
a=find(TypesNew2~=0);
EI2=mean(Traces2(a,:,:));

DiffEI=EI2-EI1;

%for i=1:98
%Traces(i,:,:)=Traces(i,:,:)-b(1,:,:);
%end
AmplitudeRange=[-250 250];
TimeRange=[8 50];
TimeRange=[0 50];
Types=ones(STraces(1));
FigureProperties=struct('FigureNumber',15,'Subplot',[2 2 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',12,'Colors',['g' 'r' 'b' 'm' 'k']);
y15=NS_PlotManySignaturesOnArrayLayout(DiffEI,neighbors,Types,ArrayID,FigureProperties,NS_GlobalConstants);

return;

figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
%name=[WritePath '\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\Traces' '_StimChannel' num2str(StimulatingChannel) '_M' num2str(MovieNumber2) '_minus_' 'M' num2str(MovieNumber1)]
print(hj, '-dtiff', name);
AmplitudeRange=[-1000 1000];
FigureProperties=struct('FigureNumber',25,'Subplot',[2 2 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',12,'Colors',['g' 'r' 'b' 'm' 'k']);
y15=NS_PlotManySignaturesOnArrayLayout(Traces1,neighbors,TypesNew1,ArrayID,FigureProperties,NS_GlobalConstants);

FigureProperties=struct('FigureNumber',35,'Subplot',[2 2 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',12,'Colors',['g' 'r' 'b' 'm' 'k']);
y15=NS_PlotManySignaturesOnArrayLayout(Traces2,neighbors,TypesNew2,ArrayID,FigureProperties,NS_GlobalConstants);

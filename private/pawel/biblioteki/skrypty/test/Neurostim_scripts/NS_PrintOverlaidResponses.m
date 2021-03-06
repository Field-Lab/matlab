function [Traces]=NS_PrintOverlaidResponses(StimulatingChannel,MovieNumber,RecElectrodes,BadChannels,ReadPath,FileName,WritePath,FigureProperties);
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
cd(ReadPath);

% 4. 
Fn=2;
%FileName
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
%figure(Fn);
%clf;
%subplot(2,2,1);
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,StimulatingChannel,0,1,'b-',NS_GlobalConstants);

% 5. Build data for PCA
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
Timings=Timings-1;
Timings0=Timings(1:100);    
Timings=Timings0;

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
neighbors0 = electrodeMap.getAdjacentsTo(StimulatingChannel,1);

neighbors=NS_RemoveBadChannels(neighbors0,BadChannels);

if RecElectrodes(1)~=0
    neighbors=RecElectrodes;
end
TimeStart=-10;
NumberOfSamples=60;
Offsets=ones(length(neighbors))*(-370);
Traces=NS_ReadManyTracesFromRaw(FileName,neighbors,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
%Traces0=Traces([1:73 75:80 82:98],:,:);
%Traces=Traces0;
STraces=size(Traces);
b=mean(Traces);
TracesOrig=Traces;
for i=1:STraces(1)
    Traces(i,:,:)=Traces(i,:,:)-b(1,:,:);
end
%Traces=TracesOrig;




FigureProperties.FigureNumber=15;
FigureProperties.LineWidth=1;
y15=NS_PlotManySignaturesOnArrayLayout(Traces,neighbors,ones(1,length(Timings)),ArrayID,FigureProperties,NS_GlobalConstants);

Amp=num2str(Amplitude);
if Amplitude<10
    AmpStr=[Amp(1) Amp(3:length(Amp))];
else
    AmpStr=[Amp(1:2) Amp(4:length(Amp))];
end

figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
%name=[WritePath '\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\Traces' '_Ch' num2str(StimulatingChannel) '_M' num2str(MovieNumber) '_' AmpStr 'uA']
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

return;

SC=size(ClustersForDiff)
for i=1:SC(1)
    EIsDiff(i,:,:)=EIs(ClustersForDiff(i,1),:,:)-EIs(ClustersForDiff(i,2),:,:);
end
size(EIsDiff)
%FigureProperties=struct('FigureNumber',18,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',14,'Colors',['n' 'r' 'b' 'm' 'k' 'g' 'c']);
FigureProperties.FigureNumber=18;
FigureProperties.AmplitudeRange=[-150 100],
EITypes=[1:SC(1)];
y18=NS_PlotManySignaturesOnArrayLayout(EIsDiff,neighbors,EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
%name=['E:\analysis\2008-03-21-0\2008_05_05\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\EI_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%fid = fopen(name,'wb')
%fwrite(fid,EIs,'double');
%fclose(fid);
%clear EIs;

%return;

FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 3 3],'TimeRange',[0 30],'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['g' 'r' 'b' 'm' 'k']);

%[PCA_CoeffsNew,TypesNew]=NS_CleanClusters(PCA_Coeffs,Types,1);

FigureProperties.Subplot=[2 2 4];
h=NS_PlotPCA_Coeffs(NumberOfClusters,TypesNew,PCA_Coeffs,FigureProperties);
%h=gca;
set(h,'XLim',PCA_XLim);
set(h,'YLim',PCA_YLim);
refresh;

FigureProperties.Subplot=[2 3 6];
a=find(TypesNew~=0);
TypesNew2=TypesNew(a);
Waveforms=Waveforms(a,:);
%h=NS_PlotManyWaveformsOnFigure(Waveforms,TypesNew2,ArrayID,FigureProperties,NS_GlobalConstants);
%h=gca;
%set(h,'XLim',Wave_XLim);
%set(h,'YLim',Wave_YLim);

%break;
subplot(2,2,3);
axis([0 100 0 100]);
h=gca;
set(h,'XTick',[]);
set(h,'YTick',[]);
set(h,'Color',[1 1 1]);
set(h,'XColor',[1 1 1]);
set(h,'YColor',[1 1 1]);
%set(h,'Visible','off');
xstart=10;
ystart=84;
ystep=9;
text(xstart,ystart+ystep,['Traces: ' num2str(length(TypesNew2)) '/' num2str(length(Types)) ' - ' num2str(length(TypesNew2)/length(Types)*100,3) '%']);
for i=1:NumberOfClusters    
    text(xstart,ystart-(0.5+i)*ystep,['Cluster ' num2str(i) ': ' num2str(length(find(TypesNew==i))) '/' num2str(length(find(Types==i))) ' - ' num2str(length(find(TypesNew==i))/length(find(Types==i))*100,3) '%']);
end


Amp=num2str(Amplitude);
if Amplitude<10
    AmpStr=[Amp(1) Amp(3:length(Amp))];
else
    AmpStr=[Amp(1:2) Amp(4:length(Amp))];
end

%return;

%figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');

%name=['E:\analysis\2008-03-21-0\2008_05_05\PCA_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\PCA' '_StimChannel' num2str(StimulatingChannel) '_M' num2str(MovieNumber) '_' AmpStr 'uA'];
%print(hj, '-dtiff', name);

figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'rotated');
%name=[WritePath '\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\Traces' '_M' num2str(MovieNumber) '_' AmpStr 'uA'];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(16);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\Traces_FullScale_' '_StimChannel' num2str(StimulatingChannel) '_M' num2str(MovieNumber) '_' AmpStr 'uA'];
%print(hj, '-dtiff', name);

figure(17);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EIs' '_M' num2str(MovieNumber) '_' AmpStr 'uA'];
%name=[WritePath '\Traces_FullScale_' '_StimChannel' num2str(StimulatingChannel) '_M' num2str(MovieNumber) '_' AmpStr 'uA'];
print(hj, '-dtiff', name);

figure(18);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EIdiff' '_M' num2str(MovieNumber) '_' AmpStr 'uA']
%name=[WritePath '\EIdiff' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(19);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EI4clust' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Cleaned_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
%print(hj, '-dtiff', name);
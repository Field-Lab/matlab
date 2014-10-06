function [Types,TypesNew,EIs,Traces]=NS_PCA_test_function_meeting(StimulatingChannel,MovieNumber,RecElectrodes,BadChannels,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);
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
NumberOfSamples=61; %51 !!
SamplesForPCA=[9:38];
PCA_channels=3; %1 - only stimulating channels; 2 - only neighbors; 3 - both stimulating one and neighbors
Dimensions=3; % ??
%NumberOfClusters=4; % ??

% 4. 
Fn=2;
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
figure(Fn);
clf;
subplot(2,2,1);
%[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,StimulatingChannel,0,1,'b-',NS_GlobalConstants);

% 5. Build data for PCA
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
%Timings=Timings-1;
Timings=Timings+7;
Timings0=Timings(TracesNumbers);
Timings=Timings0;
    
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
neighbors=[16 13];
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
Traces=TracesOrig;
%EI0=NS_CalculateEI(Traces);
%Traces=Traces-EI0;
switch PCA_channels
    case 1,
        dane=Traces(:,1,SamplesForPCA);
    case 2,
        dane=Traces(:,2:length(neighbors),SamplesForPCA);      
    case 3,
        dane=Traces(:,:,SamplesForPCA);
end
Sdane=size(dane)
dane1=NS_ConcacenateWaveforms(dane);

[Types,PCA_Coeffs,Inc]=NS_ClusterSignatures(dane1,Dimensions,NumberOfClusters);
%Types=NS_CleanClustersNew(PCA_Coeffs,Types,10);
%Types=ones(1,198);
figure(1)
plot(Types,'b*-')
%FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 2 2],'TimeRange',[0 30],'AmplitudeRange',[-200 200],'FontSize',14,'Colors',['g' 'r' 'b' 'm' 'k']);
FigureProperties.FigureNumber=Fn;
y=NS_PlotPCA_Coeffs(NumberOfClusters,Types,PCA_Coeffs,FigureProperties);
%return;
h=gca;
set(h,'XLim',[-0.2 0.3]);
set(h,'YLim',[-0.3 0.3]);
set(h,'XTick',[-0.2:0.1:0.3]);
set(h,'FontSize',30);
xlabel('1st PCA variable');
ylabel('2nd PCA variable');


PCA_XLim=get(h,'XLim');
PCA_YLim=get(h,'YLim');
Waveforms=reshape(Traces(:,1,SamplesForPCA),Sdane(1),Sdane(3));
FigureNumber=1;
AmplitudeRange=[-300 200];
TimeRange=[8 50];
TimeRange=[0 50];

FigureNumber=15;
FigureProperties.FigureNumber=15;
FigureProperties.LineWidth=1;
%FigureProperties=struct('FigureNumber',15,'Subplot',[2 2 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',12,'Colors',['n' 'r' 'b' 'm' 'k']);
a=find(Types==3);

TypesNew=NS_CleanClustersNew(PCA_Coeffs,Types,5);
y15=NS_PlotManySignaturesOnArrayLayout_SP(Traces,neighbors,TypesNew,ArrayID,FigureProperties,NS_GlobalConstants);

FigureNumber=26;
%FigureProperties=struct('FigureNumber',16,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',[-800 800],'FontSize',12,'Colors',['n' 'r' 'b' 'm' 'k']);
FigureProperties.FigureNumber=16;
y15=NS_PlotManySignaturesOnArrayLayout_SP(TracesOrig,neighbors,TypesNew,ArrayID,FigureProperties,NS_GlobalConstants);

for i=1:NumberOfClusters
    a1=find(Types==i);
    Waves1=Traces(a1,:,:);
    EIs(i,:,:)=NS_CalculateEI(Waves1);
end

%FigureProperties=struct('FigureNumber',19,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k' 'g' 'c']);
FigureProperties.FigureNumber=17;
FigureProperties.LineWidth=2;
size(EIs);
Types=[1:NumberOfClusters];
y18=NS_PlotManySignaturesOnArrayLayout_SP(EIs,neighbors,Types,ArrayID,FigureProperties,NS_GlobalConstants);

%EIs3(1,:,:)=EI2-EI3;
%EIs3(2,:,:)=EI1-EI4;
%EIs3(3,:,:)=EI4-EI3;
%EIs3(4,:,:)=EI1-EI2;

SC=size(ClustersForDiff)
for i=1:SC(1)
    EIsDiff(i,:,:)=EIs(ClustersForDiff(i,1),:,:)-EIs(ClustersForDiff(i,2),:,:);
end
size(EIsDiff)
%FigureProperties=struct('FigureNumber',18,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',14,'Colors',['n' 'r' 'b' 'm' 'k' 'g' 'c']);
FigureProperties.FigureNumber=18;
FigureProperties.AmplitudeRange=[-150 100],
EITypes=[1:SC(1)];
y18=NS_PlotManySignaturesOnArrayLayout_SP(EIsDiff,neighbors,EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
%name=['E:\analysis\2008-03-21-0\2008_05_05\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\EI_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%fid = fopen(name,'wb')
%fwrite(fid,EIs,'double');
%fclose(fid);
%clear EIs;

%return;

FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 3 3],'TimeRange',[0 30],'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['g' 'r' 'b' 'm' 'k']);

%FigureProperties.Subplot=[2 2 4];
%h=NS_PlotPCA_Coeffs(NumberOfClusters,TypesNew,PCA_Coeffs,FigureProperties);
%h=gca;
%set(h,'XLim',PCA_XLim);
%set(h,'YLim',PCA_YLim);
%refresh;
figure(123);
FigureProperties.Subplot=[2 3 6];
a=find(TypesNew~=0);
TypesNew2=TypesNew(a);
Waveforms=Waveforms(a,:);

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

hc = figure(15);
name=[WritePath '\FourClusters_Traces'];
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[10 7]);
set(hc,'PaperPosition',[0 0 10 7]);
print(hc, '-dtiff', '-r120', name);

figure(16);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\FourClusters_ble'];
%print(hc, '-dtiff', '-r120', name);

%return;

hc=figure(17);
name=[WritePath '\FourClusters_AveragedTraces'];
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[10 7]);
set(hc,'PaperPosition',[0 0 10 7]);
print(hc, '-dtiff', '-r120', name);

hc=figure(18);
name=[WritePath '\FourClusters_EIs'];
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[10 7]);
set(hc,'PaperPosition',[0 0 10 7]);
%print(hc, '-dtiff', '-r120', name);

hc=figure(101);
name=[WritePath '\FourClusters_PCA'];
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[7 7]);
set(hc,'PaperPosition',[0 0 7 7]);
print(hc, '-dtiff', '-r120', name);

return;

figure(19);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EI4clust' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Cleaned_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
%print(hj, '-dtiff', name);
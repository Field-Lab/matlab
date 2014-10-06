function [Types,TypesNew,EIs,Traces]=NS_PCA_test_function_meeting(StimulatingChannel,MovieNumber,RecElectrodes,BadChannels,ReadPath,FileName,WritePath,NumberOfClusters,ClustersForDiff,TracesNumbers,FigureProperties);
%Calculates PCA and writes two figures to hard drive: one showing two
%dimensions of PCA space, and one showing 7 electrodes with clustered
%waveform - all of that for given channel and movie number.
FontSize=20;
BoxLineWidth=2;
% 1. Global constants
ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

% 2. File, channel, movie
cd(ReadPath);

% 3. Analysis parameters (for PCA)
TimeStart=0; %-10;
NumberOfSamples=61; %51 !!
SamplesForPCA=[9:38];
PCA_channels=3; %1 - only stimulating channels; 2 - only neighbors; 3 - both stimulating one and neighbors
Dimensions=3; % ??

% 4. 
Fn=2;
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
figure(Fn);
clf;
subplot(2,2,1);

% 5. Build data for PCA
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
Timings=Timings+7;
Timings0=Timings(TracesNumbers);
Timings=Timings0;
    
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
neighbors0 = electrodeMap.getAdjacentsTo(StimulatingChannel,1);
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

STraces=size(Traces);
b=mean(Traces);
TracesOrig=Traces;
for i=1:STraces(1)
Traces(i,:,:)=Traces(i,:,:)-b(1,:,:);
end
Traces=TracesOrig;

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

%plot(Types,'b*-')

x1=0.05;
x2=0.56;
dx=0.18;
dx1=0.24;
y1=0.58;
y2=0.07;
dy=0.4;
figure(201);
clf;
subplot('Position',[0.75 y1 0.23 dy]);
y=NS_PlotPCA_Coeffs(NumberOfClusters,Types,PCA_Coeffs,FigureProperties);

h=gca;
set(h,'XLim',[-0.2 0.3]);
set(h,'YLim',[-0.3 0.3]);
set(h,'XTick',[-0.2:0.1:0.3]);
set(h,'FontSize',FontSize);
xlabel('1st PCA variable');
ylabel('2nd PCA variable');
set(h,'LineWidth',BoxLineWidth);

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
a=find(Types==3);

TypesNew=NS_CleanClustersNew(PCA_Coeffs,Types,5);

subplot('Position',[0.57 0.67 0.01 0.01]);
h=plot(1,1,'r-',1,1,'b-',1,1,'m-',1,1,'k-');
set(h,'LineWidth',2);
h=gca;
set(h,'FontSize',FontSize);
%set(h,'LineWidth',BoxLineWidth);
%set(h,'LineWidth',2);
h1=legend('neuron 1','n. 1 and 2','n. 1 and 3','n. 1, 2 and 3');
set(h1,'LineWidth',2)
p=get(h1,'Position');
p(1)=0.41;
p(2)=0.49;
set(h1,'Position',p);

subplot('Position',[x1 y1 dx dy]);
y15=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(Traces(:,1,:),neighbors(1),TypesNew,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
subplot('Position',[x1+dx1 y1 dx dy]);
y15=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(Traces(:,2,:),neighbors(2),TypesNew,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
h=ylabel('');

dfgdfgdfh=size(Traces)

for i=1:NumberOfClusters
    a1=find(Types==i);
    Waves1=Traces(a1,:,:);
    EIs(i,:,:)=NS_CalculateEI(Waves1);
end

FigureProperties.FigureNumber=17;
FigureProperties.LineWidth=2;
size(EIs)
Types=[1:NumberOfClusters];
subplot('Position',[x1 y2 dx dy]);
y18=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(EIs(:,1,:),neighbors(:,1,:),Types,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
subplot('Position',[x1+dx1 y2 dx dy]);
y18=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(EIs(:,2,:),neighbors(:,2,:),Types,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
h=ylabel('');

SC=size(ClustersForDiff)
for i=1:SC(1)
    EIsDiff(i,:,:)=EIs(ClustersForDiff(i,1),:,:)-EIs(ClustersForDiff(i,2),:,:);
end
size(EIsDiff)
FigureProperties.FigureNumber=18;
FigureProperties.AmplitudeRange=[-150 100],
EITypes=[1:SC(1)];
subplot('Position',[x2 y2 dx dy]);
y18=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(EIsDiff(:,1,:),neighbors(:,1,:),EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
h1=legend('cl. 2 - cl. 1','cl. 3 - cl. 1','cl. 4 - cl. 2','cl. 3 - cl. 2');
subplot('Position',[x2+dx1 y2 dx dy]);
y18=NS_PlotManySignaturesOnArrayLayout_SP_2011_03_01(EIsDiff(:,2,:),neighbors(:,2,:),EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
h=gca;
set(h,'LineWidth',BoxLineWidth);
h=ylabel('');

name=[WritePath '\EI_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];

return;
FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 3 3],'TimeRange',[0 30],'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['g' 'r' 'b' 'm' 'k']);

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
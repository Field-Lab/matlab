function [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(StimulatingChannel,RecordingChannels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
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

% 3. Analysis parameters (for PCA)
TimeStart=-10;
NumberOfSamples=71;
SamplesForPCA=[18:48];
PCA_channelsMode=3; %1 - only stimulating channels; 2 - only neighbors; 3 - both stimulating one and neighbors
Dimensions=3;
NumberOfClusters=4;

% 4. 
Fn=2;
[Pulse,Status]=NS_FindPulseShapeForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
figure(Fn);
clf;
subplot(2,2,1);
[Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,StimulatingChannel,0,1,'b-',NS_GlobalConstants);

% 5. Build data for PCA
[Timings,PDChunkIndex]=NS_FindPulsesTimingsForMovie(FileName,StimulatingChannel,MovieNumber,NS_GlobalConstants);
GoodChannels=NS_RemoveBadChannels(RecordingChannels,BadChannels);
Offsets=ones(length(GoodChannels))*(-370);
Traces=NS_ReadManyTracesFromRaw(FileName,GoodChannels,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
b=mean(Traces);
sb=size(b);
for i=1:sb(1)
    Traces(i,:,:)=Traces(i,:,:)-b(1,:,:);
end

%EI0=NS_CalculateEI(Traces);
%Traces=Traces-EI0;
switch PCA_channelsMode
    case 1,
        dane=Traces(:,1,SamplesForPCA);
    case 2,
        dane=Traces(:,2:length(GoodChannels),SamplesForPCA);      
    case 3,
        dane=Traces(:,:,SamplesForPCA);
end
Sdane=size(dane);
dane1=NS_ConcacenateWaveforms(dane);

% 6. Calculate PCA
[Types,PCA_Coeffs,Inc]=NS_ClusterSignatures(dane1,Dimensions,NumberOfClusters);

% 7. Plot PCA coeeficients
figure(1)
plot(Types,'b*-')
FigureProperties=struct('FigureNumber',Fn,'Subplot',[2 2 2],'TimeRange',[0 30],'AmplitudeRange',[-200 200],'FontSize',10,'Colors',['g' 'r' 'b' 'm' 'k']);
y=NS_PlotPCA_Coeffs(NumberOfClusters,Types,PCA_Coeffs,FigureProperties);
h=gca;
PCA_XLim=get(h,'XLim');
PCA_YLim=get(h,'YLim');

hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
Amp=num2str(Amplitude);
if Amplitude<10
    AmpStr=[Amp(1) Amp(3:length(Amp))];
else
    AmpStr=[Amp(1:2) Amp(4:length(Amp))];
end
name=[WritePath '\PCA_' AmpStr 'uA' '_Channel' num2str(StimulatingChannel)];
print(hj, '-dtiff', name);

% 8. Plot traces after clustering
Waveforms=reshape(Traces(:,1,SamplesForPCA),Sdane(1),Sdane(3));
FigureNumber=1;
AmplitudeRange=[-400 100];
TimeRange=[18 70];
FigureNumber=15;
FigureProperties=struct('FigureNumber',15,'Subplot',[2 2 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['g' 'r' 'b' 'm' 'k']);
y15=NS_PlotManySignaturesOnArrayLayout(Traces,GoodChannels,Types,ArrayID,FigureProperties,NS_GlobalConstants);

%figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\Traces' AmpStr 'uA' '_Channel' num2str(StimulatingChannel)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

TypesNew=Types;
EIs=0;
return;

% 9. Clean clusters and plot waveforms after cleaning
TypesNew=NS_CleanClustersNew(PCA_Coeffs,Types,5);
FigureNumber=16;
FigureProperties=struct('FigureNumber',16,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k']);
y15=NS_PlotManySignaturesOnArrayLayout(Traces,GoodChannels,TypesNew,ArrayID,FigureProperties,NS_GlobalConstants);

EIs=zeros(NumberOfClusters*(NumberOfClusters-1)/2,Sdane(2),NumberOfSamples);
EIs2=zeros(NumberOfClusters,Sdane(2),NumberOfSamples);
WaveformTypes=ones(Sdane(2));
FigureProperties=struct('FigureNumber',17,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',[-200 200],'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k']);


for i=1:NumberOfClusters
    for j=i+1:NumberOfClusters
        if i~=j            
            %cluster 1:
            a2=find(Types==i);
            Waves1=Traces(a2,:,:);
            EI1=NS_CalculateEI(Waves1);
            %cluster 2:
            a2=find(Types==j);
            Waves2=Traces(a2,:,:);
            EI2=NS_CalculateEI(Waves2);
            EI=EI2-EI1;
            
            EImax=max(max(EI));
            indexEImax=find(EI==EImax);
            EImin=min(min(EI));
            indexEImin=find(EI==EImin);
            
            if indexEImin>indexEImax
                %EI=-EI;
                nr1=i;
                nr2=j;
            else
                nr1=i; %j;
                nr2=i;
            end            
            EIs2(i,:,:)=EI1;
            
            index=(i-1)*NumberOfClusters-sum([1:1:i])+j;
            EIs(index,:,:)=EI;
            y=NS_PlotSignatureOnArrayLayout(EI,GoodChannels,WaveformTypes,ArrayID,FigureProperties,NS_GlobalConstants);
            %hj=gcf;
            %set(hj, 'PaperOrientation', 'portrait');
            %name=['D:\analysis\2008-03-21-0\PCA_EI\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber) 'Clusters' num2str(nr1) '-' num2str(nr2)];
            %print(hj, '-dtiff', name);
                   
            %figure(18);
            %FigureProperties=struct('FigureNumber',17,'Subplot',[2 3 3],'TimeRange',[0 30],'AmplitudeRange',[-200 200],'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k']);
            %EITypes=[1:NumberOfClusters*(NumberOfClusters-1)/2];
            %y15=NS_PlotManySignaturesOnArrayLayout(EIs,neighbors,EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
        end
    end
end


return;



figure(18);
FigureProperties=struct('FigureNumber',18,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',[-200 200],'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k' 'g' 'c']);
EITypes=[1:NumberOfClusters*(NumberOfClusters-1)/2];
y18=NS_PlotManySignaturesOnArrayLayout(EIs,GoodChannels,EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
%name=['E:\analysis\2008-03-21-0\2008_05_05\EI_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\EI_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
fid = fopen(name,'wb')
fwrite(fid,EIs,'double');
fclose(fid);
%clear EIs;

figure(19);
FigureProperties=struct('FigureNumber',19,'Subplot',[2 3 3],'TimeRange',TimeRange,'AmplitudeRange',AmplitudeRange,'FontSize',10,'Colors',['n' 'r' 'b' 'm' 'k' 'g' 'c']);
EITypes=[1:NumberOfClusters*(NumberOfClusters-1)/2];
y18=NS_PlotManySignaturesOnArrayLayout(EIs2,GoodChannels,EITypes,ArrayID,FigureProperties,NS_GlobalConstants);
%clear EIs2;

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

%figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');

%name=['E:\analysis\2008-03-21-0\2008_05_05\PCA_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
name=[WritePath '\PCA_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(15);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\Array_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(16);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\Array_Cleaned_Channel' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Cleaned_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(18);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EIdiff' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Cleaned_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);

figure(19);
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
name=[WritePath '\EI4clust' num2str(StimulatingChannel) 'Movie' num2str(MovieNumber)];
%name=['E:\analysis\2008-03-21-0\2008_05_05\Array_Cleaned_Channel' num2str(Channel) 'Movie' num2str(MovieNumber)];
print(hj, '-dtiff', name);
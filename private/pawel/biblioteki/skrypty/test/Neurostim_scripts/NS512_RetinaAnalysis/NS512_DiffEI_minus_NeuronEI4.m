NS_GlobalConstants=NS_GenerateGlobalConstants(500);
%UpsamplingFactor=5;

Movies=[1:2:63];

RawDataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
CrosstalkFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\CrosstalkInfo';
ClusterFilePathCore='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\ClusterFiles_protected\ClusterFile_data000_ID'; %sciezka musi byc dopelniona numerem ID neuronu
EIfilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\EIs';

%EventNumber=82;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes

CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,EventNumber)
Pattern=CrosstalkInfo(1)
Movie=CrosstalkInfo(2)
NeuronID=CrosstalkInfo(6)
PrimaryRecEl=CrosstalkInfo(7)
ChannelsPlot=electrodeMap.getAdjacentsTo(PrimaryRecEl,1)';

EI0=NS512Read_EI_File(EIfilePath,NeuronID,ChannelsPlot);
SEI=size(EI0);
EI=zeros(SEI(1),(SEI(2)-1)*UpsamplingFactor+1);
for i=1:SEI(1)    
    EI0(i,:)=EI0(i,:)-EI0(i,1);
    EI(i,:)=Upsample(EI0(i,:),UpsamplingFactor);
end

[x1,x2]=min(EI,[],2);
MostNegativeSampleIndex_NeuronEI=x2(1)

ClusterFilePath=[ClusterFilePathCore num2str(NeuronID)]
ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath);
eff=sum(ClusterIndexes-1,3);
EffVsAmp=eff(Movies,Pattern);
MoviesToAnalyze=Movies(find(EffVsAmp));
L=120*UpsamplingFactor;
t=[1:L];

%MoviesToAnalyze=[57 59];

for i=1:length(MoviesToAnalyze)-1
    %1 * * * * * Plot EIs on seven electrodes
    Movie1=MoviesToAnalyze(i);
    Movie2=MoviesToAnalyze(i+1);
    
    Movie1Indexes=find(ClusterIndexes(Movie1,Pattern,:)==1);
    Movie2Indexes=find(ClusterIndexes(Movie2,Pattern,:)==1);
    
    [DataTraces1,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie1,1500,Movie1Indexes);
    [DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie2,1500,Movie2Indexes);   
    
    % 1. Extract data for plotting differential EI on seven electrodes
    SDT1=size(DataTraces1);
    p1=DataTraces1(:,ChannelsPlot,:);
    SDT2=size(DataTraces2);
    p2=DataTraces2(:,ChannelsPlot,:);
    DiffEI0=mean(p2)-mean(p1);
    SDEI=size(DiffEI0);
    DiffEI=zeros(SDEI(1),SDEI(2),(SDEI(3)-1)*UpsamplingFactor+1);
    for kl=1:SDEI(2);
        DiffEI(1,kl,:)=Upsample(DiffEI0(1,kl,:),UpsamplingFactor);
    end
    
    % 2. Extract data for animated EI
    EI0=reshape(mean(DataTraces2)-mean(DataTraces1),SDT1(2),SDT1(3));
    SEI=size(EI0);
    for j=1:SEI(1)
         EI0(j,:)= EI0(j,:)-mean([ EI0(j,1)  EI0(j,SEI(2))]);
    end
    
    Norm=max(max(EI0')-min(EI0'));
    EI_animation=EI0/Norm*150;
    
    % 2. Find best match of neuron EI and differential EI on the primary
    % electrode
    [x3,x4]=min(DiffEI,[],3);
    MostNegativeSampleIndex_DiffEI=x4(1);
    Delay0=MostNegativeSampleIndex_NeuronEI-MostNegativeSampleIndex_DiffEI-1;
    % 3. Generate plots on seven electrodes and 51 time shits
    DataTraces=zeros(2,length(ChannelsPlot),L);    
    HalfAmplitudeWidth=[FindNegativeHalfAmplitudeWidth(EI(1,:)) FindNegativeHalfAmplitudeWidth(reshape(DiffEI(1,1,1:L),1,L))];
    NegAmplitude=[FindNegativeAmplitude(EI(1,:)) FindNegativeAmplitude(reshape(DiffEI(1,1,1:L),1,L))];
    Delay=Delay0;
    for Frame=1:51%Delay=Delay0-5*UpsamplingFactor:Delay0+5*UpsamplingFactor    
        Delay=Delay0-5*UpsamplingFactor+(Frame-1)
        figure(120);
        clf
        DataTraces(1,:,:)=DiffEI(:,:,1:L);
        if Delay>0
            DataTraces(2,:,1:min(Delay+L,length(EI))-Delay-1)=EI(:,Delay+2:min(Delay+L,length(EI))); %
        else
            DataTraces(2,:,1-Delay:L)=EI(:,1:L+Delay);
        end
                
        WaveformTypes=[1 0];
        FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',[0 L],'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        y=NS_PlotClustersOfSignaturesOnArrayLayout2(DataTraces,ChannelsPlot,WaveformTypes,500,FigureProperties,NS_GlobalConstants,1,[1:L]/20/UpsamplingFactor);    
        FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',[0 L]/20/UpsamplingFactor,'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\DiffEIminusNeuronEI\e' num2str(EventNumber) 'delay' num2str(Delay+5) 'EIs'];            
        DiffDiffEI=DataTraces(1,:,:)-DataTraces(2,:,:);
        HalfAmplitudeWidth=[HalfAmplitudeWidth FindNegativeHalfAmplitudeWidth(reshape(DiffDiffEI(1,1,:),1,L))];
        NegAmplitude=[NegAmplitude FindNegativeAmplitude(reshape(DiffDiffEI(1,1,:),1,L))];
        
        y=NS_PlotClustersOfSignaturesOnArrayLayout2(DataTraces(1,:,:)-DataTraces(2,:,:),ChannelsPlot,1,500,FigureProperties,NS_GlobalConstants,2,[1:L]/20/UpsamplingFactor);    
        
        subplot('position',[0.02 0.02 0.6 0.43]);
        h=NS512_ShowEIFrameAsCircles(EI_animation(:,Frame),500,[1:512],PrimaryRecEl,Pattern,Stim,[-1005 1005],[-505 505]);
        set(h,'Visible','off');
        h15=NS_AddFrameForArrayLayout2(500,2);
        
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\Animations\e' num2str(EventNumber) 'delay' num2str(Delay+5) 'Neuron'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r120', FigureName);                    
    end
    %DiffEI-reshape(EI,1,7,140);
end
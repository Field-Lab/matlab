NS_GlobalConstants=NS_GenerateGlobalConstants(500);
if DelaySpike==0
    UpsamplingFactor=5;
    NumberOfFrames=51;
    TR=[0 60];
else
    UpsamplingFactor=10;
    NumberOfFrames=101;
    TR=[0 120];
end

Movies=[1:2:63];
%DelaySpike
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
    EI3=reshape(mean(DataTraces2)-mean(DataTraces1),SDT1(2),SDT1(3));
    SEI=size(EI3);
    for j=1:SEI(1)
         EI3(j,:)= EI3(j,:)-mean([ EI3(j,1)  EI3(j,SEI(2))]);
    end
    
    Norm=max(max(EI3')-min(EI3'));
    EI_animation=EI3/Norm*150;
    
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
    
    for Frame=1:101
        Delay=Delay0-5*UpsamplingFactor+(Frame-1);
        DataTraces(1,:,:)=DiffEI(:,:,1:L);
        if Delay>0
            DataTraces(2,:,1:min(Delay+L,length(EI))-Delay-1)=EI(:,Delay+2:min(Delay+L,length(EI))); %
        else
            DataTraces(2,:,1-Delay:L)=EI(:,1:L+Delay);
        end
        DiffDiffEI=DataTraces(1,:,:)-DataTraces(2,:,:);
        HalfAmplitudeWidth=[HalfAmplitudeWidth FindNegativeHalfAmplitudeWidth(reshape(DiffDiffEI(1,1,:),1,L))];
        NegAmplitude=[NegAmplitude FindNegativeAmplitude(reshape(DiffDiffEI(1,1,:),1,L))];        
    end
    
    for Frame=1:101
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
        %if DelaySpike==1
        %    TR=[0 60]%[0 120]
        %else
        %    TR=[0 60];
        %end
            
        FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',TR,'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',12,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        y=NS_PlotClustersOfSignaturesOnArrayLayout2(DataTraces,ChannelsPlot,WaveformTypes,500,FigureProperties,NS_GlobalConstants,1,[1:L]/20/UpsamplingFactor);    
        FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',TR,'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',12,'Colors',['b' 'b' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        y=NS_PlotClustersOfSignaturesOnArrayLayout2(DataTraces(1,:,:)-DataTraces(2,:,:),ChannelsPlot,1,500,FigureProperties,NS_GlobalConstants,2,[1:L]/20/UpsamplingFactor);    
        
        % * * * * * Negative Amplitude vs Half-Amplitude Width plot:
        subplot('position',[0.72 0.53 0.24 0.41]);
        plot(HalfAmplitudeWidth(3:length(HalfAmplitudeWidth))/UpsamplingFactor/20,NegAmplitude(3:length(NegAmplitude)),'b-');
        hold on
        %h=plot(HalfAmplitudeWidth(i,length(HalfAmplitudeWidth))/UpsamplingFactor/20,NegAmplitude(i,length(HalfAmplitudeWidth)),'bd');
        h=plot(HalfAmplitudeWidth(Frame+2)/UpsamplingFactor/20,NegAmplitude(Frame+2),'bd');
        set(h,'MarkerSize',10);
        set(h,'MarkerFaceColor','b');    
        h=plot(HalfAmplitudeWidth(i,1)/UpsamplingFactor/20,NegAmplitude(i,1),'kd');
        set(h,'MarkerSize',10);
        set(h,'MarkerFaceColor','k');
        h=plot(HalfAmplitudeWidth(i,2)/UpsamplingFactor/20,NegAmplitude(i,2),'rd');
        set(h,'MarkerSize',10);
        set(h,'MarkerFaceColor','r');    
        axis([0 0.6 -300 0]);
        grid on;
        h=gca;
        set(h,'FontSize',12);
        xlabel('Half-amplitude Width [ms]');
        ylabel('Negative Amplitude [mV]');
        
        % * * * * Plot differential EI animation on array layout
        subplot('position',[0.02 0.02 0.45 0.4]);
         if Frame==1
            Stim=Pattern;
        else
            Stim=[];
         end
        EI_animation=EI_animation/max(max(EI_animation))*50;
        h=NS512_ShowEIFrameAsCircles(EI_animation(:,Frame),500,[1:512],PrimaryRecEl,Pattern,Stim,[-1005 1005],[-505 505]); %tutaj popatrzec
        set(h,'Visible','off');
        h15=NS_AddFrameForArrayLayout2(500,2);
        
        % * * * * Plot Neuron EI animation on array layout
        EI0=NS512Read_EI_File(EIfilePath,NeuronID,[1:512]);
        SEI=size(EI0);
        for j=1:SEI(1)
            EI0(j,:)= EI0(j,:)-mean([ EI0(j,1)  EI0(j,SEI(2))]);
        end
        EI0=EI0/max(max(EI0))*50;
        subplot('position',[0.52 0.02 0.45 0.4]);
        h=NS512_ShowEIFrameAsCircles(EI0(:,Frame+14),500,[1:512],PrimaryRecEl,Pattern,[],[-1005 1005],[-505 505]);
        set(h,'Visible','off');
        h15=NS_AddFrameForArrayLayout2(500,2);
        
        % * * * * Print figure
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\proby2\e' num2str(EventNumber) 'delay' num2str(Delay+5+1000) 'Neuron'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r120', FigureName);                    
    end
    %DiffEI-reshape(EI,1,7,140);
end
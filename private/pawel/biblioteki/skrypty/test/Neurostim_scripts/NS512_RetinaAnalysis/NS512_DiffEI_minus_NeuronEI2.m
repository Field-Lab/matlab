NS_GlobalConstants=NS_GenerateGlobalConstants(500);

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

EI=NS512Read_EI_File(EIfilePath,NeuronID,ChannelsPlot);
SEI=size(EI);
for i=1:SEI(1)
    EI(i,:)=EI(i,:)-EI(i,1);
end
[x1,x2]=min(EI,[],2);
MostNegativeSampleIndex_NeuronEI=x2(1)

ClusterFilePath=[ClusterFilePathCore num2str(NeuronID)]
ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath);
eff=sum(ClusterIndexes-1,3);
EffVsAmp=eff(Movies,Pattern);
MoviesToAnalyze=Movies(find(EffVsAmp));
L=120;
t=[1:L];

%MoviesToAnalyze=[57 59];
for i=1:length(MoviesToAnalyze)-1
    Movie1=MoviesToAnalyze(i);
    Movie2=MoviesToAnalyze(i+1);
    
    Movie1Indexes=find(ClusterIndexes(Movie1,Pattern,:)==1);
    Movie2Indexes=find(ClusterIndexes(Movie2,Pattern,:)==1);
    
    [DataTraces1,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie1,1500,Movie1Indexes);
    [DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie2,1500,Movie2Indexes);   
    
    SDT1=size(DataTraces1);
    p1=DataTraces1(:,ChannelsPlot,:);
    SDT2=size(DataTraces2);
    p2=DataTraces2(:,ChannelsPlot,:);
    DiffEI=mean(p2)-mean(p1);
    [x3,x4]=min(DiffEI,[],3);
    MostNegativeSampleIndex_DiffEI=x4(1)
    Delay0=MostNegativeSampleIndex_NeuronEI-MostNegativeSampleIndex_DiffEI-1
    
    DataTraces=zeros(2,length(ChannelsPlot),L);
    
    HalfAmplitudeWidth=[FindNegativeHalfAmplitudeWidth(EI(1,:)) FindNegativeHalfAmplitudeWidth(reshape(DiffEI(1,1,1:L),1,L))];
    NegAmplitude=[FindNegativeAmplitude(EI(1,:)) FindNegativeAmplitude(reshape(DiffEI(1,1,1:L),1,L))];
    Delay=Delay0;
    for Delay=Delay0-5:Delay0+5    
        DataTraces(1,:,:)=DiffEI(:,:,1:L);
        if Delay>0
            DataTraces(2,:,1:min(Delay+L,length(EI))-Delay-1)=EI(:,Delay+2:min(Delay+L,length(EI))); %
        else
            DataTraces(2,:,1-Delay:L)=EI(:,1:L+Delay);
        end
                
        WaveformTypes=[1 0];
        FigureProperties=struct('FigureNumber',120,'Subplot',[2 3 3],'TimeRange',[0 L],'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        y=NS_PlotClustersOfSignaturesOnArrayLayout(DataTraces,ChannelsPlot,WaveformTypes,500,FigureProperties,NS_GlobalConstants);    
        FigureProperties=struct('FigureNumber',121,'Subplot',[2 3 3],'TimeRange',[0 L],'AmplitudeRange',[min(EI(1,:))*1.5 max(EI(1,:))*1.5],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',2,'YLabel','Signal [mV]');
        
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\DiffEIminusNeuronEI\e' num2str(EventNumber) 'delay' num2str(Delay+5) 'EIs'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        %print(h, '-dtiff', '-r120', FigureName);   
    
        %y=NS_PlotClustersOfSignaturesOnArrayLayout(DiffEI(:,:,[1:L])-reshape(EI(:,[Delay0+1:Delay0+L]),1,length(ChannelsPlot),L),ChannelsPlot,1,500,FigureProperties,NS_GlobalConstants);
        DiffDiffEI=DataTraces(1,:,:)-DataTraces(2,:,:);
        HalfAmplitudeWidth=[HalfAmplitudeWidth FindNegativeHalfAmplitudeWidth(reshape(DiffDiffEI(1,1,:),1,L))];
        NegAmplitude=[NegAmplitude FindNegativeAmplitude(reshape(DiffDiffEI(1,1,:),1,L))];
        
        y=NS_PlotClustersOfSignaturesOnArrayLayout(DataTraces(1,:,:)-DataTraces(2,:,:),ChannelsPlot,1,500,FigureProperties,NS_GlobalConstants);    
        
        FigureName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\DiffEIminusNeuronEI\e' num2str(EventNumber) 'delay' num2str(Delay+5) 'Neuron'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        %print(h, '-dtiff', '-r120', FigureName);                    
    end
    %DiffEI-reshape(EI,1,7,140);
end
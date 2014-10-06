NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
% 1) Find primary electrode for each cell
NeuronIDs=[6 110 139 320 349 367 636 726 873 889 1026 1190 1283 1428 1445 1518 1732 1804 1880 1908 2017 2044 2089 2269 2299 2483 2677 2722 2900 3037 3218 3245 3425 3605 3622 3917 4055 4158 4310 4401 4578 4609 4865 4895 4907 5102 5162 5255 5389 5602 5633 5798 5884 5977 6064 6110 6263 6320 6425 6473 6503 7011 7025 7069 7249 7278 7324 7427 7461 7638];

CrosstalkFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\CrosstalkInfo'; 

Movies=[1:2:63];
Events=[150];
UpsamplingFactor=5;

RawDataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
EIsFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\EIs';
EIsFigureFolderPath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\NeuronsEIs';
ClusterFilePathCore='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\ClusterFiles\ClusterFile_data000_ID'; %sciezka musi byc dopelniona numerem ID neuronu

figure(11)
for i=1:length(Events)           
    CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,Events(i));
    Pattern=CrosstalkInfo(1);
    Movie=CrosstalkInfo(2);
    NeuronID=CrosstalkInfo(6);
    PrimaryRecEl=CrosstalkInfo(7);
    ChannelsPlot=electrodeMap.getAdjacentsTo(PrimaryRecEl,1)';
    ChannelsPlot=[1:512];
            
    EI0=NS512Read_EI_File(EIsFilePath,NeuronID,[1:512]);   % * ** * * * * * * * *
    
    
    
    ClusterFilePath=[ClusterFilePathCore num2str(NeuronID)];
    ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath);    
    eff=sum(ClusterIndexes-1,3);
    EffVsAmp=eff(Movies,Pattern);
    MoviesToAnalyze=Movies(find(EffVsAmp))
    Movie1=MoviesToAnalyze(i);
    Movie2=MoviesToAnalyze(i+1);
    
    
    Movie1Indexes=find(ClusterIndexes(Movie1,Pattern,:)==1);
    Movie1IndexesFalse=find(ClusterIndexes(Movie1,Pattern,:)~=1);
    Movie2Indexes=find(ClusterIndexes(Movie2,Pattern,:)==1);
    
    [DataTraces1,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie1,1500,Movie1Indexes);
    [DataTraces1False,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie1,1500,Movie1IndexesFalse);
    [DataTraces2,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(RawDataPath,RawDataPath,0,Pattern,Movie2,1500,Movie2Indexes);   
       
    SDT1=size(DataTraces1);
    p1=DataTraces1(:,ChannelsPlot,:);
    SDT2=size(DataTraces2);
    p2=DataTraces2(:,ChannelsPlot,:);
    
    SEI0=size(p1);
    EI0=reshape(mean(p2)-mean(p1),SEI0(2),SEI0(3));
                    
    SEI=size(EI0);
    for j=1:SEI(1)
        EI0(j,:)=EI0(j,:)-mean([EI0(j,1) EI0(j,SEI(2))]);
    end
    
    
    DiffEI=zeros(512,(SEI(2)-1)*UpsamplingFactor+1);
    for kl=1:SEI(1);
        DiffEI(kl,:)=Upsample(EI0(kl,:),UpsamplingFactor);
    end
    
    
    
    
    Norm=max(max(DiffEI')-min(DiffEI'));
    EI=EI0/Norm*250;
    %SEI=size(EI);
    for j=1:51%SEI(2)        
        clf    
        subplot('position',[0.01 0.01 0.98 0.98]);
        if j==1
            Stim=Pattern;
        else
            Stim=[];
        end
        h=NS512_ShowEIFrameAsCircles(EI(:,j),500,[1:512],PrimaryRecEl,Pattern,Stim,[-1005 1005],[-505 505]);
        set(h,'Visible','off');
    
        h15=NS_AddFrameForArrayLayout2(500,2);
        pause(1)
    
        FigureName=[EIsFigureFolderPath '\e' num2str(Events(i))  'ID' num2str(NeuronID) '_p' num2str(Pattern) '_f' num2str(j)]            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r60', FigureName);        
    end
end
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
% 1) Find primary electrode for each cell
NeuronIDs=[6 110 139 320 349 367 636 726 873 889 1026 1190 1283 1428 1445 1518 1732 1804 1880 1908 2017 2044 2089 2269 2299 2483 2677 2722 2900 3037 3218 3245 3425 3605 3622 3917 4055 4158 4310 4401 4578 4609 4865 4895 4907 5102 5162 5255 5389 5602 5633 5798 5884 5977 6064 6110 6263 6320 6425 6473 6503 7011 7025 7069 7249 7278 7324 7427 7461 7638];

CrosstalkFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\CrosstalkInfo'; 

%Events=[87 111 161 147 150 121 94 74 207];

Events=[161 147 150 151 121 94 74 87 111 207];
Events=[161];

EIsFilePath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\EIs';
EIsFigureFolderPath='C:\pawel\nauka\analiza\retina\2012-09-27-4\files\NeuronsEIs';

figure(11)
for i=1:length(Events)        
   
    CrosstalkInfo=NS512_ReadCrosstalkFile(CrosstalkFilePath,Events(i));
    Pattern=CrosstalkInfo(1);
    Movie=CrosstalkInfo(2);
    NeuronID=CrosstalkInfo(6);
    PrimaryRecEl=CrosstalkInfo(7);
            
    EI0=NS512Read_EI_File(EIsFilePath,NeuronID,[1:512]);   % * ** * * * * * * * *
    SEI=size(EI0);
    for j=1:SEI(1)
        EI0(j,:)=EI0(j,:)-mean([EI0(j,1) EI0(j,SEI(2))]);
    end
    Norm=max(max(EI0')-min(EI0'));
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
    
        FigureName=[EIsFigureFolderPath '\ID' num2str(NeuronID) '_p' num2str(Pattern) '_f' num2str(j)]            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[16 9]);
        set(h,'PaperPosition',[0 0 16 9]);  
        print(h, '-dtiff', '-r60', FigureName);        
    end
end
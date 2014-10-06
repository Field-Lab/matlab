NS_GlobalConstants=NS_GenerateGlobalConstants(61);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
for i=1:64
    ElectrodesCoordinates(i,1)=electrodeMap.getXPosition(i);
    ElectrodesCoordinates(i,2)=electrodeMap.getYPosition(i);
end

ClusterFilePath='E:\pawel\analysis\retina\2010-09-21-0\data003\ClusterFile_003_id';
NeuronsIDs=[76 227 256 271 391 406 541 616 691 736 856 901];
Electrodes=[1:8 10:24 26:56 58:64];
for i=1:26    
    [StimChannels,Amplitudes]=NS_StimulatedChannels('E:\pawel\analysis\retina\2010-09-21-0\data003',45,i,[1:64],NS_GlobalConstants);
    Amps(i)=Amplitudes;
end
StEl=[6 16 18 3 27 28 37 45 54 51 60 61];
%Colors=['b' 'k' 'g' 'r' 'c' 'm' 'y'];

ColorsRed = [0  1 0 0   0.5 0   0     0.5 0.5 0     0.5  1 1 0]
ColorsGreen=[0  0 1 0   0   0.5 0     0.5 0   0.5   0.5  1 0 1]
ColorsBlue =[0  0 0 1   0   0   0.5   0   0.5 0.5   0.5  0 1 1]

figure(101)
clf
subplot('Position',[0.12 0.15 0.1 0.1]);
for i=1:length(NeuronsIDs)        
    h=plot(1,i,'bd-');
    set(h,'Color',[ColorsRed(i) ColorsGreen(i) ColorsBlue(i)]);
    %axis([1.1 2 0 1]);
    hold on;
end

EffEl=zeros(13,64);
Thresholds100=ones(1,length(NeuronsIDs))*30;
Thresholds0=ones(length(NeuronsIDs),length(NeuronsIDs))*30;

for i=[2 3 5:12]
    ClusterFileName=[ClusterFilePath num2str(NeuronsIDs(i))]
    for R=[2 3 5:12]%1:length(StEl)
        Pattern=StEl(R);
        for Movie=1:26
            WaveformTypes=NS_ReadClusterFile(ClusterFileName,Movie,Pattern,50);
            Eff(Movie)=length(find(WaveformTypes==2))/length(WaveformTypes)*100;
        end
        
        if i==R
            il=find(Eff==100);
            Thresholds100(1,i)=il(1);
        else
            il=find(Eff>0);
            if length(il)>0
                Thresholds0(i,R)=il(1);
            end
        end
    end
            
    
    
    
    
    for P=1:0%length(Electrodes)    !!!!!!!!
        Pattern=Electrodes(P);
        for Movie=1:26
            WaveformTypes=NS_ReadClusterFile(ClusterFileName,Movie,Pattern,50);
            Eff(Movie)=length(find(WaveformTypes==2))/length(WaveformTypes)*100;
        end
        %tutaj mamy - dla danego neuronu i danego patternu - wartosc
        %efektywnosci stymulacji w funkcji movie number (Eff)
        il=find(Eff>0);
                                
        if length(il>0)
            EffEl(i,Electrodes(P))=il(1);
            subplot('Position',[ElectrodesCoordinates(Pattern,1)/600+0.45 ElectrodesCoordinates(Pattern,2)/600+0.45 0.08 0.08]);
            h=semilogx(Amps(il),Eff(il),'bd-');
            axis([0.15 1.6 0 100]);
            text(0.16,90,['el. ' num2str(Pattern)]);
            set(h,'Color',[ColorsRed(i) ColorsGreen(i) ColorsBlue(i)]);
            hold on;
            h=gca;
            set(h,'YTick',[0:20:100]);
            set(h,'XTick',[0.15 0.3 0.6 1.2]);
            if Pattern~=42
                %set(h,'XTickLabel',[]);            
                %set(h,'YTickLabel',[]);
            end
            grid on;            
        end
        
        if StEl(i)==Pattern
            subplot('Position',[ElectrodesCoordinates(Pattern,1)/600+0.45 ElectrodesCoordinates(Pattern,2)/600+0.45 0.08 0.08]);
            axis([0.15 1.6 0 100]);
            h=gca;
            set(h,'XScale','log');
            set(h,'YTick',[0:20:100]);
            set(h,'XTick',[0.15 0.3 0.6 1.2]);
                [StimChannels,Amp1]=NS_StimulatedChannels('E:\pawel\analysis\retina\2010-09-21-0\data017',Pattern,2,[1:64],NS_GlobalConstants);
                [StimChannels,Amp2]=NS_StimulatedChannels('E:\pawel\analysis\retina\2010-09-21-0\data017',Pattern,47,[1:64],NS_GlobalConstants);                
                text(0.16,70,['ID=' num2str(NeuronsIDs(i))]);
                text(0.16,50,['Amp=' num2str(Amp1,'%10.2f') '-' num2str(Amp2,'%10.2f')]);
            end
    end
    stim_el=find(EffEl(i,:));
    leg{i}=['ID=' num2str(NeuronsIDs(i)) ', el=' num2str(stim_el)]
end     

subplot('Position',[0.12 0.8 0.01 0.1]);
for i=1:length(NeuronsIDs)        
    h=plot(1,i,'bd-');
    set(h,'Color',[ColorsRed(i) ColorsGreen(i) ColorsBlue(i)]);
    axis([1.1 2 0 1]);
    hold on;
end
legend(leg)
break;
FullName=['C:\home\pawel\nauka\Stimchip_paper\obrazki\map.tif'];     
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[18.7 13]);
set(h,'PaperPosition',[0 0 18.7 13]); 
print(h, '-dtiff', '-r200', FullName);
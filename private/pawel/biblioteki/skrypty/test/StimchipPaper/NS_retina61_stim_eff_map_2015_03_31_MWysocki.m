% Kluczowa jest wywolanie funkcji NS_ReadClusterFile (linijka 107). 

clear;
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
for i=1:64
    ElectrodesCoordinates(i,1)=electrodeMap.getXPosition(i);
    ElectrodesCoordinates(i,2)=electrodeMap.getYPosition(i);
end

BoxLineWidth=2;
ElectrodeMap([1:6 16 7 8 10:11 18 12:15 17 27 19:24 26 28 29:36 37 38:44 45 46:50 51:54 55:56 60 58:59 61:64])=[1:61];

ClusterFilePath='C:\home\Pawel\nauka\StimchipPaper\2012grudzien\dyskusje\data003\ClusterFile_003_id';

NeuronsIDs=[76 227 256 271 391 406 541 616 691 736 856 901];
NeuronsIDs2=[1:12];
Electrodes=[1:8 10:24 26:56 58:64];

for i=1:26    
    [StimChannels,Amplitudes]=NS_StimulatedChannels('C:\home\Pawel\nauka\StimchipPaper\2012grudzien\dyskusje\data003',45,i,[1:64],NS_GlobalConstants);
    Amps(i)=Amplitudes;
end
StEl=[1 16 18 10 27 28 37 45 54 51 60 61];

ColorsRed = [0  1 0 0   0.5 0   0     0.5 0.5 0     0.5  1 1 0]
ColorsGreen=[0  0 1 0   0   0.5 0     0.5 0   0.5   0.5  1 0 1]
ColorsBlue =[0  0 0 1   0   0   0.5   0   0.5 0.5   0.5  0 1 1]

ColorsRed = [0 1 0 0   0 1 0 0  0 1 0 0];
ColorsGreen=[0 0 1 0   0 0 1 0  0 0 1 0];
ColorsBlue =[0 0 0 1   0 0 0 1  0 0 0 1];

ColorsRed = [0 1 0 0   0 1 0 0  0 1 0 0];
ColorsGreen=[0 0 1 0   0 0 1 0  0 0 1 0];
ColorsBlue =[0 0 0 1   0 0 0 1  0 0 0 1];

figure(101)
clf
subplot('Position',[0.12 0.1 0.1 0.1]);
for i=1:length(NeuronsIDs)        
    h=plot(1,i,'bd-');
    set(h,'Color',[ColorsRed(i) ColorsGreen(i) ColorsBlue(i)]);
    %axis([1.1 2 0 1]);
    hold on;
end

EffEl=zeros(13,64);
duzo=30;
Thresholds100=ones(1,length(NeuronsIDs))*duzo;
Thresholds0=ones(length(NeuronsIDs),64)*duzo;

Colors={'k' 'r' 'b' 'g' 'k' 'r' 'b' 'g' 'k' 'r' 'b' 'g'};
Ms={'pentagram' 'pentagram' 'pentagram' 'pentagram' 'diamond' 'diamond' 'diamond' 'diamond' 'square' 'square' 'square' 'square'};
Ms={'diamond' 'diamond' 'diamond' 'diamond' 'diamond' 'diamond' 'diamond' 'square' 'square' 'square' 'square' 'square'};

Ms={'square' 'square' 'square' 'diamond' 'diamond' 'diamond' 'pentagram' 'pentagram' 'pentagram' 'hexagram' 'hexagram' 'hexagram'};
Ms={'square' 'square' 'square' 'diamond' 'diamond' 'diamond' 'pentagram' 'pentagram' 'pentagram' 'hexagram' 'hexagram' 'hexagram'};
ColorsRed = [1 0 0   1 0 0   1 0 0  1 0 0];
ColorsGreen=[0 1 0   0 1 0   0 1 0  0 1 0];
ColorsBlue =[0 0 1   0 0 1   0 0 1  0 0 1];

Colors={'k' 'r' 'b' 'g' 'b' 'g' 'k' 'r' 'k' 'b' 'g' 'k'};

figure(200)
clf
subplot('Position',[0.06 0.18 0.8 0.75]);
%ghj=[2 5:8 10:12];
ghj=[1:12];
for k=1:length(ghj)    
    i=ghj(k);
    h2=plot([10:20],'ro-');
    set(h2,'LineStyle','none');
    hold on;
    set(h2,'Color',[ColorsRed(k) ColorsGreen(k) ColorsBlue(k)]);
    if k<10
        set(h2,'Marker',Ms{i});
    end
    set(h2,'MarkerFaceColor',[ColorsRed(k) ColorsGreen(k) ColorsBlue(k)]);
    set(h2,'MarkerSize',14);
    leg{k}=['neuron ' num2str(k)];
end
hl=legend(leg);
p=get(hl,'Position');
p(1)=0.83;
p(4)=0.001;
p(4)=0.25;
set(hl,'Position',p);
set(hl,'Box','off');
set(hl,'FontSize',18);

h=gca;
set(h,'LineWidth',BoxLineWidth);

figure(300);
kolorki=jet(1500);
clf;

[76 227 256 271 391 406 541 616 691 736 856 901];
[227 391 541 616 691 736 856];
ghj=[2 5 7 8 9 10 1 3 4 6 11 12];
for k=1:12 % dla kazdego neuronu...
    i=ghj(k) %efektywnie i=k, to przypisanei to zasz?o??z czas�w gdy i oraz k by?y r�?ne
    ClusterFileName=[ClusterFilePath num2str(NeuronsIDs(i))] % to jest nazwa Clusterfile dla danego neuronu
    for Pattern=Electrodes % dla ka?dej z 61 elektro matrycy...        
        %Pattern=StEl(R);
        for Movie=1:26 %dla ka?dego movie - czyli dal ka?dej amplitudy (jedno movie to stymulacja wszystkich elektro dla jednej amplitudy)
            WaveformTypes=NS_ReadClusterFile(ClusterFileName,Movie,Pattern,50); %kazdemu z 50 przebieg�w przypisany jest typ: 1 - czysty artefakt, 2 - spike
            Eff(Movie)=length(find(WaveformTypes==2))/length(WaveformTypes)*100; %efektywnosc stymulacji danego neuronu dana elektroda przy danej amplitudzie
        end
        
        if StEl(i)==Pattern %jesli ten pattern odpowiada primary electrode dla tego neuronu...
            il=find(Eff==100); %... znajdz amplitude, przy ktorej efektywnosc stymulacji jest 100%
            if length(il)>0
                Thresholds100(1,i)=il(1);
            end            
        else
            il=find(Eff>0);
            if length(il)>0
                Thresholds0(i,Pattern)=il(1);
            end
        end
    end
            
    figure(200);
    hold on;
    T0=Thresholds0(i,:);
    S1=find(T0<duzo)
    h1=plot(ElectrodeMap(S1),Amps(T0(S1)),'ro');
    set(h1,'Color',[ColorsRed(k) ColorsGreen(k) ColorsBlue(k)]);
    if (k<10)
        set(h1,'Marker',Ms{k});
    end
    set(h1,'MarkerSize',14);
    set(h1,'LineStyle','none');
    set(h1,'LineWidth',2)
    
    if Thresholds100(i)<duzo
        h2=plot(ElectrodeMap(StEl(i)),Amps(Thresholds100(i)),'ro');
        set(h2,'Color',[ColorsRed(k) ColorsGreen(k) ColorsBlue(k)]);
        if k<10
            set(h2,'Marker',Ms{k});
        end
        set(h2,'MarkerFaceColor',[ColorsRed(k) ColorsGreen(k) ColorsBlue(k)]);
        set(h2,'MarkerSize',14);   
    end
    
    h3=gca;
    axis([0 62 0 1.7])
    grid on;
    %set(h3,'XTick',sort(ElectrodeMap(StEl([2 3 5:12]))))
    set(h3,'XTick',sort(ElectrodeMap(StEl(ghj))))
    set(h3,'FontSize',20);
    xlabel('Electrode ID');
    ylabel('Threshold current [\muA]');        
        
    leg{i}=['ID=' num2str(NeuronsIDs2(i))];
end     

%W porownaniu do wersji 2013+07+18, teraz rysujemy na wykresach dane dla
%data002 (bez synaptic blockers) oraz data005 (z synaptic blockers), nadale
%dla tkanki 2010-09-14-0.
clear
NS_GlobalConstants=NS_GenerateGlobalConstants(512);

data=2;

if data==5
    NeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\data005\2010-09-14-0\data005\data005.neurons');
    PathForImages='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\data005\other_plots\';
    FullName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\data005\histograms4\GaussParameters.bin']; % pokoj 109
    AllPatternsFinal=NS512_AllPatternsInExperiment('G:\data\2010-09-14-0\movie005',NS_GlobalConstants)
else
    if data==2
        NeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.neurons');
        PathForImages='C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2013_10_02\data002\other_plots\';
        FullName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\histograms5\GaussParameters.bin']; % pokoj 109
        AllPatternsFinal=NS512_AllPatternsInExperiment('G:\data\2010-09-14-0\movie002',NS_GlobalConstants)
    end
end

AllPatternsFinal=NS512_AllPatternsInExperiment('G:\data\2010-09-14-0\movie005',NS_GlobalConstants);

fid=fopen(FullName,'r','ieee-le');
a=fread(fid,'double');
GaussParameters=reshape(a,length(a)/5,5); %NeuronID, electrode, Amplitude, tau, sigma

neurons=GaussParameters(:,1);
uniqueneurons=unique(neurons);
NeuronID=4666
Events=find(neurons==NeuronID);
%sigmas(Events)
break
delays=GaussParameters(:,4);
sigmas=GaussParameters(:,5);
figure(100)
hist(sigmas/20,[1:2:130]/20);
h=xlabel('\sigma [ms]');
set(h,'FontSize',10)
h=ylabel('N');
set(h,'FontSize',10)
grid on
h=gca
set(h,'FontSize',10)
FullImageName=[PathForImages 'sigmas_histogram.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);

figure(101);
h=plot(delays/20,sigmas/20,'bd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','b')
h=xlabel('Latency [ms]');
set(h,'FontSize',10)
ylabel('\sigma [ms]');
set(h,'FontSize',10)
grid on
h=gca
set(h,'FontSize',10)
set(h,'XLim',[0 30])
FullImageName=[PathForImages 'latency_vs_sigma.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);

Neurons=GaussParameters(:,1);
Electrodes=GaussParameters(:,2);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
figure(100);
clf
hold on;
for i=1:length(Neurons)
    SeedEl(i) = NeuronFile.getNeuronIDElectrode(Neurons(i));
    XSeed(i)=electrodeMap.getXPosition(SeedEl(i));
    YSeed(i)=electrodeMap.getYPosition(SeedEl(i));
    XStim(i)=electrodeMap.getXPosition(Electrodes(i));
    YStim(i)=electrodeMap.getYPosition(Electrodes(i));
    %h15=plot([XStim(i) XSeed(i)]-XStim(i),[YStim(i) YSeed(i)]-YStim(i),'bd-');
    h15=plot([XStim(i) XSeed(i)],[YStim(i) YSeed(i)],'b-');
    if sigmas(i)<7
        set(h15,'Color','r');
    end
    h16=plot([XSeed(i)],[YSeed(i)],'gd');
    set(h16,'MarkerFaceColor','g')
    set(h16,'MarkerSize',20)
    Distances(i)=sqrt((XSeed(i)-XStim(i))^2+((YSeed(i)-YStim(i))^2));
end
grid on
axis([-1000 1000 -500 500]);

figure(102);
h=plot(Distances,delays/20,'bd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','b');
hold on;
DirectStimulation=find(sigmas<7);
InDirectStimulation=find(sigmas>6);

h=plot(Distances(DirectStimulation),delays(DirectStimulation)/20,'rd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','r');
h=xlabel('Distance [\mum]');
set(h,'FontSize',10)
h=ylabel('Latency [ms]');
set(h,'FontSize',10)
grid on
h=gca
set(h,'FontSize',10)
set(h,'YLim',[0 30])
FullImageName=[PathForImages 'distance_vs_latency.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);
   
figure(103);
h=plot(Distances,sigmas/20,'bd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','b');
hold on;
h=plot(Distances(DirectStimulation),sigmas(DirectStimulation)/20,'rd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','r');
h=xlabel('Distance [\mum]');
set(h,'FontSize',10)
h=ylabel('\sigma [ms]');
set(h,'FontSize',10)
grid on
h=gca
set(h,'FontSize',10)
%set(h,'YLim',[0 30])
FullImageName=[PathForImages 'distance_vs_sigma.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);

ShortLines=1;

figure(108)
clf
hold on
for i=1:length(uniqueneurons)
    EventsForNeuron=find(Neurons==uniqueneurons(i));
    DirectStimOfNeuron=0;
    InDirectStimOfNeuron=0;
    for j=1:length(EventsForNeuron);
        if find(DirectStimulation==EventsForNeuron(j))
            DirectStimOfNeuron=DirectStimOfNeuron+1;
        else
            if find(InDirectStimulation==EventsForNeuron(j))
            	InDirectStimOfNeuron=InDirectStimOfNeuron+1;
            end
        end
    end
    h=plot(DirectStimOfNeuron,InDirectStimOfNeuron,'bd');
    set(h,'MarkerFaceColor','b');
end
h=gca;
set(h,'XLim',[-1 15]);
set(h,'YLim',[-1 30]);
xlabel('No. of direct stim.');
ylabel('No. of indirect stim.');

FullImageName=[PathForImages 'DirectVsIndirect_neurons.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[3 6]);
set(h,'PaperPosition',[0 0 3 6]); 
print(h, '-dtiff', '-r240', FullImageName);

break


%now, the connectivity plot. Zmienna summary: jesli rowna 0, to wtedy
%rysujemy normlane mapy po??cze? na?o?óne na layout matrycy. Jesli równa 1,
%to wszystki linie po??cze? s?przesuwane tak, aby mia?y wspólne pocz?tki.
summary=0;

if summary==1 
    n1=0;               %tutaj okreslamy, czy na rysunku b?d? zaznaczone elektrody
    string_ext='_sum';   %tutaj okreslamy, czy w nazwie pliku b?dzie rozszerzenie
else
    n1=512;
    string_ext='';
end

NumberOfSteps=100;
figure(104);
clf
hold on
for i=1:length(DirectStimulation)
    x1=XStim(DirectStimulation(i));
    x2=XSeed(DirectStimulation(i));
    y1=YStim(DirectStimulation(i));
    y2=YSeed(DirectStimulation(i));
    
    %here, for the summary plot - all the lines having common starting
    %point:
    if summary==1
        x2=x2-x1;
        y2=y2-y1;
        x1=0;
        y1=0;
    end
    
    distance=sqrt((x1-x2)^2+(y1-y2)^2);
    if distance~=0
        StepsForShortLine=ceil(40/distance*NumberOfSteps);
        %h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[1 0 0],[0 0 1],NumberOfSteps,[1:NumberOfSteps]);        
        %h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[1 0 0],[0 0 1],NumberOfSteps,[1:StepsForShortLine]);
        h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[1 0 0],[1 0 0],NumberOfSteps,[NumberOfSteps-StepsForShortLine:NumberOfSteps]);
    end
end

for i=1:n1
    x=electrodeMap.getXPosition(i);
    y=electrodeMap.getYPosition(i);
    h18=plot(x,y,'ko');
    set(h18,'MarkerSize',3);
    if find(AllPatternsFinal==i)
        set(h18,'color','r');
        set(h18,'MarkerSize',8);
        set(h18,'LineWidth',2);
    end
end
h=gca;
set(h,'Visible','off');
        
%grid on
axis([-1000 1000 -500 500]);
FullImageName=[PathForImages 'connectivity_direct' string_ext '.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 5]);
set(h,'PaperPosition',[0 0 10 5]); 
print(h, '-dtiff', '-r240', FullImageName);


figure(105);
%clf
figure(104)
hold on
for i=length(InDirectStimulation):-1:1
    x1=XStim(InDirectStimulation(i));
    x2=XSeed(InDirectStimulation(i));
    y1=YStim(InDirectStimulation(i));
    y2=YSeed(InDirectStimulation(i));
    
    if summary==1
        x2=x2-x1;
        y2=y2-y1;
        x1=0;
        y1=0;
    end
    
    distance=sqrt((x1-x2)^2+(y1-y2)^2);
    if distance~=0
        StepsForShortLine=ceil(40/distance*NumberOfSteps);
        %h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[1 0 0],[0 0 1],NumberOfSteps,[1:NumberOfSteps]);
        %h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[1 0 0],[0 0 1],NumberOfSteps,[1:StepsForShortLine]);
        h17=NS512_ConnectivityLine([x1 x2],[y1 y2],[0 0 1],[0 0 1],NumberOfSteps,[NumberOfSteps-StepsForShortLine:NumberOfSteps]);
    end
end
for i=1:n1
    x=electrodeMap.getXPosition(i);
    y=electrodeMap.getYPosition(i);
    h18=plot(x,y,'ko');
    set(h18,'MarkerSize',3);
    if find(AllPatternsFinal==i)
        set(h18,'color','r');
        set(h18,'MarkerSize',8);
        set(h18,'LineWidth',2);
    end
end
h=gca;
set(h,'Visible','off')
        
%grid on
axis([-1000 1000 -500 500]);
FullImageName=[PathForImages 'connectivity_indirect_exp' string_ext '.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 5]);
set(h,'PaperPosition',[0 0 10 5]); 
print(h, '-dtiff', '-r240', FullImageName);
        

ElectrodesToAnalyze=setdiff(Electrodes,SeedEl);
figure(106)
clf
hold on
for i=1:length(Neurons)
    if find(ElectrodesToAnalyze==Electrodes(i))
        h15=plot([XStim(i) XSeed(i)],[YStim(i) YSeed(i)],'b-');
        if sigmas(i)<7
            set(h15,'Color','r');
        end
        h16=plot([XSeed(i)],[YSeed(i)],'gd');
        set(h16,'MarkerFaceColor','g')
        set(h16,'MarkerSize',20)
        text([XSeed(i)],[YSeed(i)],num2str(Neurons(i)));
    end
end

figure(107);

g1=hist(neurons,[1:8000]);
NeuronWithManyStimulations=find(g1>6);
for i=1:length(NeuronWithManyStimulations)
    Events=find(neurons==NeuronWithManyStimulations(i));
    S1=sigmas(Events);
    DirectivityIndex(i)=length(find(S1<7))/length(S1);    
    
    clf;
    hold on;
    for j=1:length(Events)
        h15=plot([XStim(Events(j)) XSeed(Events(j))],[YStim(Events(j)) YSeed(Events(j))],'b-');
        text(XStim(Events(j)),YStim(Events(j)),num2str(j));
        if sigmas(Events(j))<40 %7
            set(h15,'Color','r');
        end        
    end
    axis([-1000 1000 -500 500]);
end

figure(108)
clf
hold on
for i=1:length(uniqueneurons)
    EventsForNeuron=find(Neurons==uniqueneurons(i));
    DirectStimOfNeuron=0;
    InDirectStimOfNeuron=0;
    for j=1:length(EventsForNeuron);
        if find(DirectStimulation==EventsForNeuron(j))
            DirectStimOfNeuron=DirectStimOfNeuron+1;
        else
            InDirectStimOfNeuron=InDirectStimOfNeuron+1;
        end
    end
    plot(DirectStimOfNeuron,InDirectStimOfNeuron,'bd');
end

figure(109)
hist(DirectivityIndex,100)

NeuronID=5566
Events=find(neurons==NeuronID);
sigmas(Events)
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
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\other_plots\sigmas_histogram.tif'];
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
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\other_plots\latency_vs_sigma.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);

Neurons=GaussParameters(:,1);
Electrodes=GaussParameters(:,2);
NeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile('C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\Vision_output\2010-09-14-0\data002\data002.neurons');
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
for i=1:length(Neurons)
    SeedEl = NeuronFile.getNeuronIDElectrode(Neurons(i));
    XSeed=electrodeMap.getXPosition(SeedEl);
    YSeed=electrodeMap.getYPosition(SeedEl);
    XStim=electrodeMap.getXPosition(Electrodes(i));
    YStim=electrodeMap.getYPosition(Electrodes(i));
    Distances(i)=sqrt((XSeed-XStim)^2+((YSeed-YStim)^2));
end

figure(102);
h=plot(Distances,delays/20,'bd');
set(h,'MarkerSize',3)
set(h,'MarkerFaceColor','b');
hold on;
DirectStimulation=find(sigmas<7);
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
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\other_plots\distance_vs_latency.tif'];
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
FullImageName=['C:\home\Pawel\nauka\analiza\SlicesTTX\2010-09-14-0\analysis_2012_04_10\other_plots\distance_vs_sigma.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[4 3]);
set(h,'PaperPosition',[0 0 4 3]); 
print(h, '-dtiff', '-r240', FullImageName);


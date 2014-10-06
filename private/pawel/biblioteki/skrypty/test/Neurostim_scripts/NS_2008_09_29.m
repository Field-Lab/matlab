DataPath='D:\analysis\2008-08-26-0\data006_proba3';
ArtifactDataPath='D:\analysis\2008-08-26-0\data011_proba3';
ClusterFilePath='D:\analysis\2008-08-26-0\data006_proba3\clusters006_3';
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
ArtifactSubtraction=1;
FontSize=24;

Thresholds(1,:)=[1.2 0.6 0.31];
Thresholds(2,:)=[1.1 0.51 0.28];
Thresholds(3,:)=[1.01 0.45 0.28];
Thresholds(4,:)=[1.01 0.45 0.23];
Thresholds(5,:)=[1.8 0.81 0.45];
Thresholds(6,:)=[0.76 0.38 0.23];
Thresholds(7,:)=[1.5 0.68 0.38];
Thresholds(8,:)=[1.2 0.55 0.3];
Thresholds(9,:)=[1.01 0.5 0.28];
Durations=[50 100 200];
figure(101);
clf
figure(201)
clf
for i=1:9
    figure(101);
    plot(Durations,Thresholds(i,:),'bd-');
    hold on;
    figure(201);
    plot(Durations,Thresholds(i,:).*Durations,'bd-');
    hold on;
end
figure(101)
h=xlabel('Pulse durarion [\mus]');
set(h,'FontSize',FontSize);
h=ylabel('Threshold [\muA]');
set(h,'FontSize',FontSize);
h=gca;
set(h,'FontSize',FontSize);
grid on

figure(201)
h=xlabel('Pulse durarion [\mus]');
set(h,'FontSize',FontSize);
h=ylabel('Threshold [pC]');
set(h,'FontSize',FontSize);
h=gca;
set(h,'FontSize',FontSize);
grid on

Movies=[43:3:76];
Currents=[0.76 0.83 0.91 1.01 1.10 1.20 1.30 1.51 1.61 1.81 2.01 2.21];
PatternNumbers=[41 22 54 37 40 39 8 7 5];
%Channels=[6 7 8 38 39 40];
electrodes=[7 13 42 11 6 5 39 38 40];
Channels=[1:64];
%Amplitudes=zeros(length(PatternNumbers),length(Movies));

for i=1:length(PatternNumbers)
    Amp=NS_2008_10_01(DataPath,ArtifactDataPath,ClusterFilePath,NS_GlobalConstants,ArtifactSubtraction,Movies,PatternNumbers(i),Channels,electrodes(i));
    Amplitudes(i,:)=Amp/0.44;
    %plot(Currents,Amp,'bd-');
    %hold on;
end
figure(102);
clf
for i=1:length(PatternNumbers)
    figure(102);
    plot(Currents,Amplitudes(i,:),'bd-');
    hold on;
    h=text(2.25,Amplitudes(i,12),[num2str(PatternNumbers(i)) '/' num2str(electrodes(i))]);
    set(h,'FontSize',FontSize);
end
grid on;

h=xlabel('Stimulation current amplitude [\muA]');
set(h,'FontSize',FontSize);
h=ylabel('Signal amplitude on distant electrode [\muV]');
set(h,'FontSize',FontSize);
h=gca;
set(h,'FontSize',FontSize);

PatternNumbers=[37 37 37 37 37];
%Channels=[6 7 8 38 39 40];
electrodes=[10 11 12 15 26];
Channels=[1:64];
%Amplitudes=zeros(length(PatternNumbers),length(Movies));

for i=1:length(PatternNumbers)
    Amp=NS_2008_10_01(DataPath,ArtifactDataPath,ClusterFilePath,NS_GlobalConstants,ArtifactSubtraction,Movies,PatternNumbers(i),Channels,electrodes(i));
    Amplitudes2(i,:)=Amp/0.44;
    %plot(Currents,Amp,'bd-');
    %hold on;
end
figure(103);
clf
for i=1:length(PatternNumbers)
    figure(103);
    plot(Currents,Amplitudes2(i,:)/Amplitudes2(i,12),'bd-');
    hold on;
    %h=text(2.25,Amplitudes(i,12),[num2str(PatternNumbers(i)) '/' num2str(electrodes(i))]);
    %set(h,'FontSize',FontSize);
end
grid on;

h=xlabel('Stimulation current amplitude [\muA]');
set(h,'FontSize',FontSize);
h=ylabel('Signal amplitude on distant electrode (normalized)');
set(h,'FontSize',FontSize);
h=gca;
set(h,'FontSize',FontSize);

FigureProperties=struct('FigureNumber',103,'Subplot',[2 3 3],'TimeRange',[5 25],'AmplitudeRange',[-400 400],'FontSize',13,'Colors',['k' 'r' 'b' 'k' 'g' 'm' 'c'],'LineWidth',2,'YLabel','input signal [\muV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(EIsDiff/0.44,Channels,[1:length(Movies)],1,FigureProperties,NS_GlobalConstants);
%Amp=NS_2008_10_01(DataPath,ArtifactDataPath,ClusterFilePath,NS_GlobalConstants,ArtifactSubtraction,Movies,PatternNumber,Channels,electrode);
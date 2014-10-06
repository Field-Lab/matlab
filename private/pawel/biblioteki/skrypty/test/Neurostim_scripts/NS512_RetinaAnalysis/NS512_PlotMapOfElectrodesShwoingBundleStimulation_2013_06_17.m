clear
NS_GlobalConstants=NS_GenerateGlobalConstants(500);
ArrayID=500;
Movies=[1:2:63];

cd C:\pawel\nauka\analiza\retina\Sept2012
day=24;
piece=0;

% Read in the somatic thresholds
ThresholdFilePath1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\threshold_files\2012-09-' num2str(day) '-' num2str(piece) '\stim_scan\thresholds_1']
ThresholdFilePath2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\threshold_files\2012-09-' num2str(day) '-' num2str(piece) '\stim_scan\thresholds_2']
[Neurons Electrodes Amplitudes]=NS512_ReadThresholdFiles(ThresholdFilePath1,ThresholdFilePath2);
ElectrodesForHistogram_Somas=find(Amplitudes>0)

%ThresholdFileName=['Thresholds' num2str(day) '_' num2str(piece)]
load(['Thresholds' num2str(day) '_' num2str(piece)])
%load Thresholds27_4;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

MinimumElectrodeAboveThreshold=4;
ThresholdForSignalAboveThresholdOnEnoughElectrodes=zeros(1,512);

figure(1)

for i=1:512
    i
    clf
    s1=Thresholds(i,:);
    Electrodes=find(s1>0);
    s1(Electrodes)=(40-s1(Electrodes));
    s1(i)=15;
    
    NS512_ShowEIAsCircles(s1'*4,500,[1:512],i,[-1005 1005],[-505 505]);
    FullName=['C:\pawel\nauka\analiza\retina\2012-09-27-4\summary4\figures_maps\map_' num2str(day) '_' num2str(piece) '_p' num2str(i)];            ;
    %h=gcf;
    %set(h,'PaperUnits','inches');
    %set(h,'PaperSize',[8 4]);
    %set(h,'PaperPosition',[0 0 8 4]); 
    %print(h, '-dtiff', '-r200', FullName);
end
function [Neurons Electrodes Amplitudes]=NS512_ReadThresholdFiles(ThresholdFilePath1,ThresholdFilePath2)
%ThresholdFilePath1=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_1']
fid1=fopen(ThresholdFilePath1,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons1=NeuronInformation(:,1);
Electrodes1=NeuronInformation(:,3);
Amplitudes1=NeuronInformation(:,4);

%ThresholdFilePath2=['C:\pawel\nauka\analiza\retina\2012-09-27-4\files\stim_scan\thresholds_2']
fid1=fopen(ThresholdFilePath2,'r');
b=fread(fid1,inf,'double');
fclose(fid1);
NeuronInformation=reshape(b,length(b)/4,4);
Neurons2=NeuronInformation(:,1);
Electrodes2=NeuronInformation(:,3);
Amplitudes2=NeuronInformation(:,4);

Neurons=[Neurons1' Neurons2']';
Electrodes=[Electrodes1' Electrodes2']';
Amplitudes=[Amplitudes1' Amplitudes2']';

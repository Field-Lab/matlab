function Execution = ExecuteDetection(ChannelSignal, ThresholdStandard)

clc

tic

ChannelSignal= ChannelSignal - mean(ChannelSignal);


% this line removes offset
ThresholdSpike = (-1)*ThresholdStandard * std(ChannelSignal);
% determines threshold according to standard deviation
ThresholdHisteresis= (-1)*(ThresholdStandard-1) * std(ChannelSignal); %histereza w przyk³adowych danych nie ma sensu!
% determines histeresis according to standard deviation

%Execution
[Spikes, Histeresis] = SpikesDetection(ThresholdSpike,ThresholdHisteresis, ChannelSignal);

ArrayToMexFile = Spikes+Histeresis;

Execution=(histereza(ArrayToMexFile));

PlotDataPartTopValue =6000;
PlotDataPartBottomValue = 5500;

PlotData(ChannelSignal, Execution*ThresholdSpike, ThresholdSpike, ThresholdHisteresis); 
%PlotDataPart(ChannelSignal, Execution*ThresholdSpike, PlotDataPartBottomValue, PlotDataPartTopValue, ThresholdSpike, ThresholdHisteresis);
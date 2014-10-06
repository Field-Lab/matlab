function t = PlotData(Signal ,SpikesDetected, ThresholdSignal, ThresholdHisteresis)
x = [1:length(Signal)];
% determines vector x, its values will be compared to values of input
% vectors in order to plot fig2
%ThresholdSignal = ThresholdSignal * ones(1,length(Signal));
% vector with constant value = threshold, and length of input vectors 
figure;
plot(x,Signal);
% plots signal
hold all
plot(x,SpikesDetected,'rs');
% plots spikes
hold all
plot(x,ThresholdSignal,'g');
% plots threshold
hold all
plot(x, ThresholdHisteresis,'g');
%plots histeresis
t = 91;
function t = PlotDataPart(Signal ,SpikesDetected, BottomValue, TopValue, ThresholdSignal, ThresholdHisteresis)
tic



figure(2);
x = [1:length(Signal)];
% determines vector x, its values will be compared to values of input
% vectors in order to plot fig2
ThresholdSignal = ThresholdSignal * ones(1,length(Signal));
% vector with constant value = threshold, and length of input vectors 
x = x(BottomValue:TopValue);
ThresholdSignal = ThresholdSignal(BottomValue:TopValue);

ThresholdHisteresis = ThresholdHisteresis * ones(1,length(Signal));
ThresholdHisteresis = ThresholdHisteresis(BottomValue:TopValue);

Signal = Signal(BottomValue:TopValue);
SpikesDetected = SpikesDetected(BottomValue:TopValue);
% this is part where section to be ploted is chosen
plot(x,Signal);
% plots signal
hold all
plot(x,SpikesDetected,'rs');
% plots spikes
hold all
plot(x,ThresholdSignal,'g');
% plots threshold
hold all
plot(x,ThresholdHisteresis,'g');
toc
t = 0;
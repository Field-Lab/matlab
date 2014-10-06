function [SignalSpike, HisteresisSpike] = SpikesDetection(ThresholdSpike,ThresholdHisteresis,ChannelData)
tic
ChannelDataSpike = ChannelData - ThresholdSpike;
ChannelDataHisteresis = ChannelData - ThresholdHisteresis;
SpikeSign = sign(ChannelDataSpike);
HisteresisSign = sign(ChannelDataHisteresis);
% sign - decides wheather signal is above or below threshold
SpikeDerivative = diff(SpikeSign);
HisteresisDerivative = diff(HisteresisSign);
% deriv shows the moment when sign is changed <=> threshold is passed
SignalSpike = abs(SpikeDerivative) - SpikeDerivative;
HisteresisSpike =abs(HisteresisDerivative) +HisteresisDerivative;
% spikes shows the moment when threshold is passed from value above to
% value below (note you are working on values below 0 so function detects
% fist change of sign, if those where valus above it wolud detect the
% second change! => if you want to work on abs of signal then exchange
% minus with plus
SignalSpike =(1/4)*[0; SignalSpike];
HisteresisSpike =(1/2)*[0; HisteresisSpike];
toc
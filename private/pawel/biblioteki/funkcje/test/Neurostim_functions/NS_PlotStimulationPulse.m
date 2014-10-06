function [Amplitude,PlotPointer]=NS_PlotStimulationPulse(Pulse,Status,Channel,StartTime,Style,NS_GlobalConstants);
%Amplitude - the highest abolute value of the current in microamps

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

t=[StartTime:StartTime+length(Pulse)-1]/Fs*1000;
CurrentStep=CurrentRanges(Status.ChannelsStatus(Channel).range+1)/127;
[t0,s0,h]=plotPH(t,Pulse(1,:).*Pulse(3,:)*CurrentStep,20,Style);
grid on;

xlabel('t [ms]');
ylabel('\mu A');

Amplitude=max(abs(Pulse(1,:).*Pulse(3,:)*CurrentStep));
PlotPointer=h;
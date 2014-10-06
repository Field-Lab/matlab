function Amplitude=NS_PlotStimulationPulse(Pulse,Status,Channel,StartTime,TimeStep,PulsesNumber,Multiplier,Style,NS_GlobalConstants);
%Amplitude - the highest abolute value of the current in microamps

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

CurrentStep=CurrentRanges(Status.ChannelsStatus(Channel).range+1)/127;
for i=1:PulsesNumber
    StartTime=StartTime+(i-1)*TimeStep;
    t=[StartTime:StartTime+length(Pulse)-1]/Fs*1000;
    plotPH(t,Pulse(1,:).*Pulse(3,:)*CurrentStep,20,Style);
    hold on;
end
grid on;

xlabel('t [ms]');
ylabel('\mu A');

Amplitude=max(abs(Pulse(1,:).*Pulse(3,:)*CurrentStep));
function [Name,Channels,AmpOut]=NS_PatternAmplitudes(Patterns,PatternsIndexes,Status,PatternNumber,NS_GlobalConstants)
%Amplitude - the highest absolute value of the current in microamps

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

IndexEnd=PatternsIndexes(PatternNumber);
if PatternNumber==1
    IndexStart=1;
else
    IndexStart=PatternsIndexes(PatternNumber-1)+1;
end

Name=[];
Channels=[];
for i=IndexStart+1:IndexEnd
    Channel=Patterns(i).channel;
    CurrentStep=CurrentRanges(Status.ChannelsStatus(Channel).range+1)/127;
    %CurrentStep=1;
    
    Data=Patterns(i).data;
        
    Amin=min(Data(1,:).*Data(3,:)*CurrentStep);
    Amax=max(Data(1,:).*Data(3,:)*CurrentStep);
    Amplitude=max(abs(Amax),abs(Amin));
    
    Amp=num2str(Amplitude);
    if Amplitude<10
        AmpStr=[Amp(1) Amp(3:length(Amp))];
    else
        AmpStr=[Amp(1:2) Amp(4:length(Amp))];
    end
    
    if abs(Amin)>abs(Amax)
        pol='n';
    else
        pol='p';
    end   
        
    %Amplitude=max(abs(Data(1,:).*Data(3,:)*CurrentStep));
    if Amplitude>0
        Name=[Name '_e' num2str(Channel) '_r' num2str(Status.ChannelsStatus(Channel).range) '_' pol AmpStr 'uA'];
        Channels=[Channels Channel];
        AmpOut=Amplitude;
    end
end


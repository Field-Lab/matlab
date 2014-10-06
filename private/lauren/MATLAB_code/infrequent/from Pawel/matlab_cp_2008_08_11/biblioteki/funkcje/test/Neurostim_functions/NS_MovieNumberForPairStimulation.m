function [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);

%Filenames(1)=struct('Filename','001',Movies,[12:34]);
[PatternNumber,number]=NS_PatternNumber(Pattern,Indexes);

SI=size(Indexes);
%MovieNumber=(AmplitudeNumber-1)*SI(1)+PatternNumber
for i=1:length(Filenames)
    a=Filenames(i).Amplitudes
    i;
    b=find(a==AmplitudeNumber)
    if b
        FileNumber=i
    end
end
M=Filenames(FileNumber).Amplitudes(1)
AmplitudeNumber=AmplitudeNumber-M+1
MovieNumber=(AmplitudeNumber-1)*SI(1)+PatternNumber
FileName=Filenames(FileNumber).Filename

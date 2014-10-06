function p=NS_PCAAndPlotsForPairStimulation(ReadPath,ParentWriteDir,Channels,ChannelsPairIndexes,Amplitudes,Indexes,Filenames,BadChannels);

Channels
ChannelsPairIndexes
Channel1=Channels(ChannelsPairIndexes(1))
Channel2=Channels(ChannelsPairIndexes(2))
p=[num2str(Channel1) 'and' num2str(Channel2)];
mkdir(ParentWriteDir,p);
PD=[ParentWriteDir '\' p];

Pattern=zeros(1,length(Channels));

y=mkdir(PD,'np');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=-1;
Pattern(ChannelsPairIndexes(2))=1;
WritePath=[PD '\np']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end

y=mkdir(PD,'nn');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=-1;
Pattern(ChannelsPairIndexes(2))=-1;
WritePath=[PD '\nn']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end

y=mkdir(PD,'pn');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=1;
Pattern(ChannelsPairIndexes(2))=-1;
WritePath=[PD '\pn']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end

y=mkdir(PD,'pp');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=1;
Pattern(ChannelsPairIndexes(2))=1;
WritePath=[PD '\pp']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end


y=mkdir(PD,'p0');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=1;
%Pattern(ChannelsPairIndexes(2))=-1;
WritePath=[PD '\p0']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end

y=mkdir(PD,'n0');
Pattern=zeros(1,length(Channels));
Pattern(ChannelsPairIndexes(1))=-1;
%Pattern(ChannelsPairIndexes(2))=1;
WritePath=[PD '\n0']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel1,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end


y=mkdir(PD,'0p');
Pattern=zeros(1,length(Channels));
%Pattern(ChannelsPairIndexes(1))=1;
Pattern(ChannelsPairIndexes(2))=1;
WritePath=[PD '\0p']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel2,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end

y=mkdir(PD,'0n');
Pattern=zeros(1,length(Channels));
%Pattern(ChannelsPairIndexes(1))=-1;
Pattern(ChannelsPairIndexes(2))=-1;
WritePath=[PD '\0n']
for i=Amplitudes
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,i,Indexes,Filenames);       
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(Channel2,Channels,MovieNumber,BadChannels,ReadPath,FileName,WritePath);
end
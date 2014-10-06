% dane z 2008-05-20, poszatkowane w kilku plikach (niestety)
BadChannels=[1 9 25 28 31 33 37 41 57 64];
BadChannels=[9 25 28 31 33 37 41 57 64];
Channels=[50 46 49 53 52 51 48];
%Channels=[6 1 4 12 11 7 3];
GoodChannels=[];
for i=Channels
    active=1;
    for j=BadChannels
        if i==j
            active=0;
        end
    end
    if active==1
        GoodChannels=[GoodChannels i];
    end
end
Channels=GoodChannels;

Filenames(1)=struct('Filename','000','Amplitudes',[1:11]);
Filenames(2)=struct('Filename','001','Amplitudes',[12:34]); % dane z 2008-05-20!!! poszatkowane w kilku plikach (niestety)

%Filenames(1)=struct('Filename','003','Amplitudes',[1:11]);
%Filenames(2)=struct('Filename','006','Amplitudes',[12:34])

Indexes=NS_IndexesForPairStimulation(5);
Pattern=[1 0 0 0 -1 0 0];
AmplitudeNumber=14;
'dfhdfh'
[FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames)

ReadPath='D:\2008-05-20-0\stimulation';
ParentDir='D:\analysis\2008-05-20-0\2008_05_25';
%y=NS_PCAAndPlotsForPairStimulation(ReadPath,ParentDir,Channels,[1 5],[8:14],Indexes,Filenames,BadChannels)

Amplitudes=[10:16];
for i=6
    ChannelsPair=[2 1]
    y=NS_PCAAndPlotsForPairStimulation(ReadPath,ParentDir,Channels,ChannelsPair,Amplitudes,Indexes,Filenames,BadChannels);
end
    
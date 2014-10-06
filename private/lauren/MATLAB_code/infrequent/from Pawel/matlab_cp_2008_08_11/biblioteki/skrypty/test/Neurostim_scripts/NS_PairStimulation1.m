% dane z 2008-05-20, poszatkowane w kilku plikach (niestety)
BadChannels=[1 9 25 28 31 33 37 41 57 64];
Channels=[50 46 49 53 52 51 48];
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
Channels=GoodChannels

Filenames(1)=struct('Filename','000','Amplitudes',[1:11]);
Filenames(2)=struct('Filename','001','Amplitudes',[12:34]); % dane z 2008-05-20!!! poszatkowane w kilku plikach (niestety)

Indexes=NS_IndexesForPairStimulation(5);
Pattern=[1 0 0 0 -1 0 0];
AmplitudeNumber=14;
'dfhdfh'
[FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames)

ReadPath='D:\2008-05-20-0\stimulation';
ParentDir='D:\analysis\2008-05-20-0\proba';
y=NS_PCAAndPlotsForPairStimulation(ReadPath,ParentDir,Channels,[1 4],[9:13],Indexes,Filenames,BadChannels)

Pattern=[1 0 0 0 -1 0 0];
%WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\np';
for AmplitudeNumber=9:8
    Pattern=[0 0 0 0 -1 0 0];
    %Pattern=[-1 0 0 1 0 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\n0';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(52,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    break;
    
    Pattern=[0 0 0 0 1 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\p0';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(52,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[1 0 0 0 0 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\0p';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[-1 0 0 0 0 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\0n';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[1 0 0 0 1 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\pp';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[1 0 0 0 -1 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\np';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[-1 0 0 0 1 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\pn';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
    Pattern=[-1 0 0 0 -1 0 0];
    WritePath='D:\analysis\2008-05-20-0\proba\stim52and50rec52\nn';
    [FileName,MovieNumber]=NS_MovieNumberForPairStimulation(Pattern,AmplitudeNumber,Indexes,Filenames);
    [Types,TypesNew,EIs,Traces]=NS_PCA_Pattern_Stimulation(50,[52 51 50 53 55],MovieNumber,BadChannels,ReadPath,FileName,WritePath);
    
end
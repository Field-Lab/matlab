function [SuccIfSucc,SuccIfFail,PrevSucc,PrevFail]=NS_FindConditionalEfficacy(RawDataPath,RawFileNumber,ClusterFilePath,Patterns,Movie,NS_GlobalConstants);
%The function finds the efficacy of stimulation by given pattern under
%codition that the previous pattern initiated spike or not. For given
%pattern P1, the previous pattern P0 (the one that was applied during stimulation
%30 ms before pattern P1) is identified, the succesfull stimulation events
%for P0 are found, and the succesful events for P! following successful
%events of P) are found. The same way, the efficacy of P! under condition
%that P0 did NOT elicit spike is found.
%NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

ClusterIndexes=NS_ReadClusterFileAll(ClusterFilePath);
ClusterIndexesNeg=ClusterIndexes; %when the previous pattern did not elicit spike
ClusterIndexesPos=ClusterIndexes; % when previous pattern did elicit spike

Timings1=Patterns;
Timings2=Patterns;
cd (RawDataPath);
ChunkData=NS_MovieData(RawFileNumber,Movie,NS_GlobalConstants);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
MovieData
for i=1:length(Patterns)
    a=find(MovieData==Patterns(i))
    Timings1(i)=MovieData(a(1)-1);
    Timings2(i)=MovieData(a(2)-1);
end
Timings1
SuccIfSucc=zeros(1,length(Patterns));
SuccIfFail=zeros(1,length(Patterns));
PrevSucc=zeros(1,length(Patterns));
PrevFail=zeros(1,length(Patterns));
for i=1:length(Patterns)
    T=Timings1(i)
    PreviousPattern=Patterns(find(Timings1==T-600))
    if PreviousPattern
        a=find(Patterns==PreviousPattern);
        if length(a)==0
            SuccIfSucc(i)=-1;
            SuccIfFail(i)=-1;
        else
            PreviousSucc=find(ClusterIndexes(Movie,PreviousPattern,1:50)~=1);
            CurrentSuccIfPrevSucc=find(ClusterIndexes(Movie,Patterns(i),PreviousSucc)~=1);
            PrevSucc(i)=length(PreviousSucc);
            SuccIfSucc(i)=length(CurrentSuccIfPrevSucc)/length(PreviousSucc);
            PreviousFail=find(ClusterIndexes(Movie,PreviousPattern,1:50)==1);
            PrevFail(i)=length(PreviousFail);
            CurrentSuccIfPrevFail=find(ClusterIndexes(Movie,Patterns(i),PreviousFail)~=1);
            SuccIfFail(i)=length(CurrentSuccIfPrevFail)/length(PreviousFail);    
        end    
    end
end
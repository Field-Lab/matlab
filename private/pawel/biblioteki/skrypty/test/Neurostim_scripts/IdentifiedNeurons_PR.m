sid1=size(id1);
IfOnes=zeros(1,sid1(1));

for i=1:64
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,i,1,[1:512],NS_GlobalConstants);
    sss=find(StimChannels==393);
    if length(sss)>0
        i
    end
end

for i=150
    electrode=id1(i,1);
    pattern=id1(i,2);
    movie=id1(i,3);    
    [StimChannels,Amplitudes]=NS_StimulatedChannels(DataPath,pattern,movie,[1:512],NS_GlobalConstants);
    StimChannels;
    sss=find (StimChannels==electrode);
    if length(sss)>0
        IfOnes(i)=1;
    end
end

soma=find(IfOnes==1);
l=id1(soma,3);

axon=find(IfOnes==0);
l1=id1(axon,3);
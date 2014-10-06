NS_GlobalConstants=NS_GenerateGlobalConstants(61);
Filename='D:\analysis\2008-08-26-0\data008_proba4\clusters008';

PrimaryEl=44;
P1Start=293;
Patterns1=[P1Start:P1Start+5];

Movies=[69:5:104];
Curves=zeros(6,length(Movies);

for i=1:6
    PatternNumber=Patterns1(i);
    for j=1:length(Movies)
        WaveformTypes=NS_ReadClusterFile(FileName,Movies(j),PatternNumber);
        a=find(WaveformTypes==1); % it is assumed that in the cluster files, for each clustering the cluster number 1 cirrespond to "artifact only"
        Curves(i,j)=100-length(a);
    end
end
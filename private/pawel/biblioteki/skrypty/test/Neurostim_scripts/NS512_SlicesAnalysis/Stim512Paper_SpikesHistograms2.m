AllTraces=zeros(400,600);
AllDelays=[];

Patterns=[13];

Movies=[8:8:136];

% 1. Znajdz ciekawe elektrody
InterestingElectrodes=NS512_SpikesAnalysis_FindElectrodesWithTimeLockedResponses('C:\pawel\nauka\512paper\SpikesAnalysis',13,[8:8:136]);

% 2. 
for m=10:length(Movies)
    Movie=Movies(m);
    fid=fopen(['C:\pawel\nauka\512paper\SpikesAnalysis\sp_p13m' num2str(Movie)],'r','ieee-le');
    a=fread(fid,'int32');
    b=reshape(a,length(a)/3,3);
    fclose(fid);    
    
    SpikesIDs=find(b(:,1)==electrode);

    Delays=b(SpikesIDs,3);
    AllDelays=[AllDelays' Delays']';
    Histograms=hist(Delays,[1:600]);
    CFs=NS512_FitWithMultiGauss([1:600],Histograms*50);

    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData('C:\pawel\nauka\512paper\ProprocessedData','C:\pawel\nauka\512paper\ProprocessedData',0,13,Movie,0,0);    
    s=reshape(DataTraces(:,electrode,:),50,600);
    AllTraces((m-10)*50+1:(m-10)*50+50,:)=s;
    figure(1)
    subplot(3,3,m-9);
    plot(s');
    figure(2)
    subplot(3,3,m-9);
    hist(Delays,[1:600]);
end

HistogramsAll=hist(AllDelays,[10:10:600]);
CFs=NS512_FitWithMultiGauss([1:60],HistogramsAll);
figure(1)
subplot(3,3,9);
plot(AllTraces');
figure(2)
subplot(3,3,9);
hist(AllDelays,[10:10:600]);
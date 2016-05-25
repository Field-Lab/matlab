electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;

X=zeros(1,512);
Y=X;
for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end

figure(1)
clf
%hold on
N=0;
for Pattern=AllPatternsUsed(8:64)
    fid=fopen(['C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2010-09-14-0\data005\GaussesFiles600points\AllGausses_p' num2str(Pattern)],'r');
    a=fread(fid,'double');
    fclose(fid);
    
    b=reshape(a,length(a)/5,5);  
    sb=size(b)
    ChannelsToExclude=electrodeMap.getAdjacentsTo(Pattern,Radius)';
    %g=zeros(1,1);
    g=[];
    for c=ChannelsToExclude
        l=find(b(:,2)==c);
        g=[g' l']'
    end            
    HistogramsToPlot=setdiff([1:sb(1)],g);            
    
    Distances=sqrt((X(b(HistogramsToPlot,1))- X(b(HistogramsToPlot,2))).^2+(Y(b(HistogramsToPlot,1))- Y(b(HistogramsToPlot,2))).^2);    N=N+length(HistogramsToPlot);
    %HistogramsToExclude=find(b(:,2)==ChannelsToExclude)
    subplot(2,2,1)
    hold on
    axis([0 30 0 6])
    plot(b(HistogramsToPlot,4)/20,b(HistogramsToPlot,5)/20,'bd');
    subplot(2,2,2)
    hold on
    axis([0 2000 0 30])
    plot(Distances,b(HistogramsToPlot,4)/20,'bd');
    subplot(2,2,3)
    hold on
    axis([0 2000 0 6])
    plot(Distances,b(HistogramsToPlot,5)/20,'bd')
end

subplot(2,2,1)
axis([0 30 0 6])
grid on

subplot(2,2,2)
axis([0 2000 0 30])
grid on

subplot(2,2,3)
axis([0 2000 0 6])
grid on
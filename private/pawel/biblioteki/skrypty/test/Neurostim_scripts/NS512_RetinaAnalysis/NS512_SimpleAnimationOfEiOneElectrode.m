NS_GlobalConstants=NS_GenerateGlobalConstants(500);
DataPath='F:\analiza\retina\2012-09-27-4\files\scan_new';
Movies=[1:2:63];

Pattern=182;
RecElectrodes=[122 123 124];
RecElectrodes=[150 214];
FigurePath=['C:\pawel\nauka\analiza\retina\2012-09-27-4\proby\p' num2str(Pattern) '.gif'];

%figure(10);
%clf
t=[1:130]/20;
for i=1:length(Movies)
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movies(i),0,0);
    DataTraces=DataTraces0(1:50,[1:512],SamplesToAnalyze);
    d=reshape(mean(DataTraces),512,130);   
    means=mean(d1(:,111:130),2);
    for k=1:512
        d(k,:)=d(k,:)-means(k);
    end
    plot(t,d(RecElectrodes,:)');   
    axis([0 4 -500 500]);
    
    h=gcf;
    frame = getframe(h);
    im = frame2im(frame);
    [imind,map] = rgb2ind(im,256);
    if i == 1
        imwrite(imind,map,FigurePath,'gif', 'Loopcount',1);
    else
        imwrite(imind,map,FigurePath,'gif','WriteMode','append','DelayTime',0.2);
    end
    %pause(1);
end
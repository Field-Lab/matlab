Chns=[233 234 241 242 243 249 250];

%T(1,:,:)=Data(Chns+1,7:106)/N;
%g=mean(Spike);
%T(2,:,:)=Spike(1,Chns,:);

FigureProperties=struct('FigureNumber',3,'TimeRange',[0 80],'AmplitudeRange',[-150 60],'FontSize',16,'Colors',['k' 'r' 'b' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','signal [\muV]');
%y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(T,[Chns],[1 2],500,FigureProperties,NS_GlobalConstants,[242]);

g1=DataTraces(artifacts,242,:);
g2=DataTraces(spikes,242,:);

g=zeros(90,1,100);
g(1:47,1,:)=g1;
g(48:90,1,:)=g2;

FigureProperties=struct('FigureNumber',3,'TimeRange',[0 80],'AmplitudeRange',[-750 -200],'FontSize',16,'Colors',['k' 'r' 'b' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','signal [\muV]');

types=ones(1,90);
types(1:47)=1;
types(48:90)=2;
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(g/0.84,[242],types,500,FigureProperties,NS_GlobalConstants,[242]);

break;
h=gca;
set(h,'XTickLabel',{'0' '' '1' '' '2' '' '3' '' '4'});
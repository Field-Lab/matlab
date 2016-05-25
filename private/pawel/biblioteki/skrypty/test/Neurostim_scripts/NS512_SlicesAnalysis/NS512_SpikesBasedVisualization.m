electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

DataPath='E:\analysis\2010-09-14-0\SpikeFiles';

Pattern=13;
Movie=112;

FigurePath=['C:\pawel\nauka\analiza\slices\2010-09-14-0\analysis_2014_08_19\SpikesPropagation' num2str(Pattern) 'm' num2str(Movie)];

fid=fopen([DataPath '\sp_p' num2str(Pattern) 'm' num2str(Movie)],'r','ieee-le');
a=fread(fid,'int32');
SpikesData=reshape(a,length(a)/3,3);
fclose(fid);

for i=[1:512]%[1:64 321:512]
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end

Channels=[1:512]
TimePoints=[2:4:598];
SpikeTimesHistograms=zeros(512,length(TimePoints));

for Channel=Channels
    SpikesForChannel=find(SpikesData(:,1)==Channel);
    SpikeTimes=SpikesData(SpikesForChannel,3);
    SpikeTimesHistograms(Channel,:)=hist(SpikeTimes,TimePoints);
end

%colormap(jet(50))
cmp=colormap(jet(50));
figure(1)
clf
hold on
for Frame=1:length(TimePoints)
    Frame
    clf
    hold on
    for Channel=Channels
        SpikeNumber=SpikeTimesHistograms(Channel,Frame);
        if SpikeNumber>0
            if SpikeNumber>50
                SpikeNumber=50;
            end
            h=plot(X(Channel),Y(Channel),'bo');
            set(h,'Markersize',12)
            set(h,'Color',cmp(SpikeNumber,:));
            set(h,'MarkerFaceColor',cmp(SpikeNumber,:));
        else
            h=plot(X(Channel),Y(Channel),'bo');
            set(h,'Markersize',2);
        end
    end
    text(-950,450,num2str(Frame*0.2))
    axis([-1000 1000 -500 500]);
    
    h=gcf;
    set(h,'PaperUnits','inches');
    set(h,'PaperSize',[9 16]);
    set(h,'PaperPosition',[0 0 9 16]); 
    print(h, '-dtiff', '-r120', [FigurePath 'f' num2str(Frame)); 
    %h=gcf;
    %frame = getframe(h);
    %im = frame2im(frame);
    %[imind,map] = rgb2ind(im,256);
    %if Frame == 1
    %    imwrite(imind,map,FigurePath,'gif', 'Loopcount',1);
    %else
    %   imwrite(imind,map,FigurePath,'gif','WriteMode','append','DelayTime',0);
    %end
end
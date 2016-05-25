load C:\home\Pawel\nauka\512stim_paper\SpikesAnalysis\2013-12-12-3\data001\AllSpikes2 AllSpikes2
SAS=size(AllSpikes2);
sp10=zeros(SAS(1),SAS(2));
sp25=sp10;
sp40=sp10;
AllSpikes=zeros(1,SAS(1));
for Amplitude=1:SAS(1)
    Spikes=reshape(AllSpikes2(Amplitude,:,:),SAS(2),SAS(3));   
    AllSpikes(Amplitude)=sum(sum(Spikes));
    %[P1,E1]=find(Spikes>25);
    for p=1:SAS(2)
        Spikes=reshape(AllSpikes2(Amplitude,p,:),1,SAS(3));   
        e10=find(Spikes>10);
        sp10(Amplitude,p)=length(e10);
        e25=find(Spikes>25);
        sp25(Amplitude,p)=length(e25);
        e40=find(Spikes>40);
        sp40(Amplitude,p)=length(e40);
    end
end

DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc';
movies=[16:18:448]
for i=1:length(movies)
    amp1(i)=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],13,movies(i),NS_GlobalConstants);
end

subplot('position',[0.77 0.57 0.19 0.38]);
%plot(amp1,sum(sp25,2)/SAS(2),'gd-',amp1,sum(sp40,2)/SAS(2),'rd-');
hold on
[ax,p1,p2]=plotyy(amp1,AllSpikes/(128*512),amp1,sum(sp10,2)/SAS(2));
set(p1,'Marker','o');
set(p1,'Color','k');
set(p2,'Marker','d');
set(p2,'Color','b');
h10=gca
set(h10,'FontSize',FontSize);
set(ax(1),'XScale','log')
set(ax(2),'XScale','log')
set(ax(1),'YColor','k')
set(ax(2),'YColor','b')
%set(ax(2),'HitTest','on')
h10=ylabel(ax(1),'Average number of spikes per electrode')
set(h10,'FontSize',FontSize);
h10=ylabel(ax(2),'Average number of activated electrodes')
set(h10,'FontSize',FontSize);

%set(ax(2),'YLabel','dfg')

set(ax(1),'YLim',[0 25]);
set(ax(1),'YTick',[0:5:25])
set(ax(2),'YLim',[0 300]);
set(ax(2),'YTick',[0:60:300])
set(ax(1),'XTick',[0.4:0.1:1 1.4 2 3 4])
set(ax(1),'XTickLabel',{'' '0.5' '' '0.7' '' '' '1' '1.4' '2' '3' '4'})
set(ax(2),'XTickLabel',{''})
set(ax(2),'FontSize',FontSize)

%set(ax(1),'YTick',[0:2000:120000])
s1=get(ax(1),'YLim')
s2=get(ax(2),'YLim')


plot(amp1,s1/s2*sum(sp25,2)/SAS(2),'gd-',amp1,s1/s2*sum(sp40,2)/SAS(2),'rd-');
xlabel('Stimulation amplitude [\muA]');

legend('Number od spikes','Threshold=20%','Threshold=50%','Threshold=80%','Location','NorthWest');
grid on
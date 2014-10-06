ID=96;
MovieID=5;
StimEl=127;

fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\RawDataPieces\ID' num2str(ID) '_stim' num2str(StimEl)],'r','ieee-le'); 
a=fread(fid,'int32');
fclose(fid)

t=[0:9999]/20;
s=reshape(a,16,50,10000);
s1=reshape(s(MovieID,:,:),50,10000);
figure(16)
clf
h=plot(t,s1');
set(h,'Color','b');
hold on

%trzeba znalezc te spiki, ktore zostaly znalezione przez Vision i narysowac
%je na czerwono:

IDs=[96 216 217 289];

for i=1:4

    ID=IDs(i);
fid=fopen(['D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data005\Matlab\dane\ID=' num2str(ID)],'r','ieee-le'); 
a=fread(fid,'int32');
l=length(a);
b=reshape(a,5,l/5);
fclose(fid);

SpikesForMovie=find(b(1,:)==MovieID);
SpikesForElectrode=find(b(3,:)==StimEl);

for rep=1:0%50
    SpikesForRepetition=intersect(find(b(2,:)==rep),SpikesForMovie);
    if SpikesForRepetition
        SignalForRepetition=s1(rep,:);
        for spike=1:length(SpikesForRepetition)
            Latency=b(5,SpikesForRepetition(spike));
            SpikeWaveform=SignalForRepetition(Latency-20:Latency+40);
            plot([Latency-21:Latency+39]/20,SpikeWaveform-i*150,'r');
            %text(230,-50-i*150,['ID=' num2str(IDs(i))])
            text(Latency/20,-100-i*150,num2str(rep));
        end
    end
end

end
axis([4 100 -200 50])
FontSize=14
h=gca
set(h,'FontSize',FontSize);
h=xlabel('Time [ms]')
set(h,'FontSize',FontSize);
h=ylabel('Signal [\muV]')
set(h,'FontSize',FontSize);

grid on

FullImageName=['C:\home\Pawel\nauka\analiza\report_2014_01_23\Traces4_5.tif'];
h=gcf;
set(h,'PaperUnits','inches');
set(h,'PaperSize',[10 3]);
set(h,'PaperPosition',[0 0 10 3]); 
print(h, '-dtiff', '-r240', FullImageName);
    
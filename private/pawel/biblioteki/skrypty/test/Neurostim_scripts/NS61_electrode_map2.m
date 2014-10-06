electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
figure(300)
clf;

MarkedElectrodes=[1 10 11 16 17 18 20 24 27 28 30 34 43 45 51 53 54 60]
MovieNumbers=[14 15 17 9 12 23 12 15 15 16 11 18 10 6 18 13 9 17]

MarkedElectrodesNew=[16 17 20 30 32 37 45 51 54 41 64 10 24];
for i=1:length(MarkedElectrodesNew)
    index=find(MarkedElectrodes==MarkedElectrodesNew(i))
    MovieNumbersNew(i)=MovieNumbers(index);
end
MT=20;
for i=1:64    
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
    h=text(X(i),Y(i),num2str(i));
    j=find(MarkedElectrodesNew==i & MovieNumbersNew<MT);
    if j
        set(h,'Color','r');
    end
    set(h,'FontSize',18);
end
axis([-300 300 -300 300])

break;
%plot(X,Y,'bd')

figure(1)
plot([0:0.1:100])
axis([0 1000 0 1000])
break;
FontSize=26;
h=text(100,900,'A)');
set(h,'FontSize',FontSize);

h=text(200,900,'B)');
set(h,'FontSize',FontSize);

h=text(300,900,'C)');
set(h,'FontSize',FontSize);

h=text(400,900,'D)');
set(h,'FontSize',FontSize);


FontSize=30;
h=text(100,700,'A)');
set(h,'FontSize',FontSize);

h=text(300,700,'B)');
set(h,'FontSize',FontSize);

h=text(500,700,'C)');
set(h,'FontSize',FontSize);

h=text(700,700,'D)');
set(h,'FontSize',FontSize);

FontSize=36;
h=text(100,500,'A)');
set(h,'FontSize',FontSize);

h=text(300,500,'B)');
set(h,'FontSize',FontSize);

h=text(500,500,'C)');
set(h,'FontSize',FontSize);

h=text(700,500,'D)');
set(h,'FontSize',FontSize);
break;

name='C:\home\pawel\2010\Stimchip_paper\literki';
hc=figure(1);
set(hc,'PaperUnits','inches');
set(hc,'PaperSize',[10 7]);
set(hc,'PaperPosition',[0 0 10 7]);
print(hc, '-dtiff', '-r120', name);
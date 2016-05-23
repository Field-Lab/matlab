electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
figure(102)
clf;

%MarkedElectrodes=PatternsToUse
MarkedElectrodes=[193:448];
MarkedElectrodes=[61 334  331 328 ]
%MarkedElectrodes=AllPatterns%[63 64  79 127 128 129 131 213 273 447 ];
%MarkedElectrodes=ElectrodesWithSpikes
%MovieNumbers=[14 15 17 9 12 23 12 15 15 16 11 18 10 6 18 13 9 17]
clear X;
clear Y;
for i=[1:512]%[1:64 321:512]
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
    h=text(X(i),Y(i),num2str(i));
    j=find(MarkedElectrodes==i);
    set(h,'FontSize',12);
    if j
        set(h,'Color','r');
        set(h,'FontSize',14);
    end
    %set(h,'FontSize',14);
end
axis([-950 950 -460 460])

break
h=gca
set(h,'Visible','off');
FullName=['D:\Home\Pawel\analysis\2010-09-11-0\ForAlan\ElectrodeMap512'];            
        h=gcf;
        set(h,'PaperUnits','inches');
        set(h,'PaperSize',[20 10]);
        set(h,'PaperPosition',[0 0 12 10]); 
        print(h, '-dtiff', '-r200', FullName);
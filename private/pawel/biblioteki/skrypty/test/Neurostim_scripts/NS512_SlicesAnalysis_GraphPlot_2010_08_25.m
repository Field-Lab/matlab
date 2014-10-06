Pattern=197;
%Pattern=126;
Movie=153;
%Movie=113;
RecElectrodes=[81 82 102 122 166 193 209 242 296];
Pairs=zeros(length(RecElectrodes),2);
Pairs(:,2)=RecElectrodes';
Pairs(:,1)=Pattern;

figure(101)
Dst=NS512_ElectrodesDistance(197,RecElectrodes,500)';
Delays=[2.2 6.5 2.0 1.7 4.5 2.5 2.5 3.0 2.0]-0.3;
y=NS512_PlotGraphWithVelocities(Pairs,ArrayID,Dst./Delays'/100);    
h=gcf;
FullName='G:\analysis\2010-07-29-0\figures2\Graf';
set(h,'PaperUnits','inches');
set(h,'PaperSize',[19 9]);
set(h,'PaperPosition',[0 0 19 9]); 
print(h, '-dtiff', '-r80', FullName);
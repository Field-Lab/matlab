figure(1);
clf;
plot([0:0.01:1]);
axis([0 1 0 1]);

h=text(0.1,0.9,'A  B  C  D');
set(h,'FontSize',28);
h=text(0.45,0.905,'(A)  (B)  (C)  (D)');
set(h,'FontSize',28);
h=text(0.5,0.8,'(   )  (   )  (  )  (   )');
set(h,'FontSize',28);


h=text(0.1,0.75,'A)  B)  C)  D)');
set(h,'FontSize',30);

h=text(0.1,0.55,'A)  B)  C)  D)');
set(h,'FontSize',40);

h=gcf;
FullName=['C:\home\pawel\nauka\Stimchip_paper\obrazki\letters'];
set(h,'PaperUnits','inches');
set(h,'PaperSize',[7 7]);
set(h,'PaperPosition',[0 0 7 7]);  
print(h, '-dtiff', '-r400', FullName);
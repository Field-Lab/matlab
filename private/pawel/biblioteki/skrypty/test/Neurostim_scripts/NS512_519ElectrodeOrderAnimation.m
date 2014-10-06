ElectrodeOrder=AllPatterns;

X0=electrode_positions(519);
X=X0(ElectrodeOrder,:);

a=X(:,1);
b=X(:,2);

figure(3);
axis([0 33 0 17]);
clf;
for i=1:519
    h=plot(a(i,1),b(i,1),'bd');
    set(h,'MarkerSize',15);
    set(h,'MarkerEdgeColor','r');
    set(h,'MarkerFaceColor','r');
    axis([-400 400 -400 400]);
    pause(0.2);
    set(h,'MarkerSize',6);
    set(h,'MarkerEdgeColor','b');
    set(h,'MarkerFaceColor','none');

    hold on;
    %axis([0 33 0 17]);
end
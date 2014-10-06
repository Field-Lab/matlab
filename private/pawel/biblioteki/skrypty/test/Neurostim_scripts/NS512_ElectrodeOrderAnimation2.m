X=NS512_OptimalElectrodeOrder();
a=X(:,1);
b=X(:,2);

figure(3);
axis([0 33 0 17]);
clf;
hold on;
for i=1:512
    if find(ElectrodesCombinedOrder==i)
        plot(a(i,1),b(i,1),'bd');     
        pause(0.25); 
    end       
    axis([0 33 0 17]);
end
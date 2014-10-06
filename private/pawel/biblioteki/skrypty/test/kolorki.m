figure(400)
c=jet(1000);
for i=1:1000
    h=plot([i,i],[5 10],'r-');
    set(h,'Color',c(i,:));
    hold on;
end
    
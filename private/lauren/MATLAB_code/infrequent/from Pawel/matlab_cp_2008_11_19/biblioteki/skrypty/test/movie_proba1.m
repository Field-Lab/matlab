a=rand(1,10);
t=[1:length(a)];
figure(11);
%h=plot(a,'bo');

for i=1:20
    clf;
    b=a*sin(i/5);
    %h=plot(b,'bo');
    for j=1:length(b) 
        signal=b(j);
        h=plot(t(j),signal,'bo');
        hold on;
        if signal<0
            cl=[1 0 0];
        else
            cl=[0 0 1];
        end
        cl
        set(h,'MarkerEdgeColor',cl);
        set(h,'MarkerFaceColor',cl);
        set(h,'MarkerSize',10); %abs(signal)*i);
    end
    axis([1 10 0 1]);
    M(i)=getframe;
    refresh;
    pause(0.5);
end 
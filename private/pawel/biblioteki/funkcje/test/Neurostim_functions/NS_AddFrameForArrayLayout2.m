function h=NS_AddFrameForArrayLayout2(ArrayID,LineWidth);

if ArrayID==1
    X=[-300 -150 150 300 150 -150];
    Y=[0 300 300 0 -300 -300];
end

if ArrayID==500
    X=[-1000 -1000 1000 1000 -1000];
    Y=[-500 500 500 -500 -500];
end
h=plot(X,Y,'k-');
set(h,'LineWidth',LineWidth);
    %for i=1:length(X)-1
    %    h=plot([X(i) X(i+1)],[Y(i) Y(i+1)],'k-');
    %    set(h,'LineWidth',LineWidth);
    %nd
    %h=plot([X(i+1) X(1)],[Y(i+1) Y(1)],'k-');
    %set(h,'LineWidth',LineWidth);
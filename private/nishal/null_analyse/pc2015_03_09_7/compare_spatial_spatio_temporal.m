
on=logical([1,0,1,0,1,0,1]);
off=logical([0,1,0,1,0,1,0]);

xx = firing_rate_log%firing_rate_log%PSTH_var_log;
figure('Color','w');
plot(xx(on,1),xx(on,2),'r+');
hold on
plot(xx(off,1),xx(off,2),'r*');
hold on
plot(xx(on,1),xx(on,3),'b+');
hold on
plot(xx(off,1),xx(off,3),'b*');
hold on
plot([0,max(xx(:,1))],[0,max(xx(:,1))],'g');
xlabel('Original');
ylabel('Null');
h=legend('Spatial, ON cells','Spatial, OFF cells','Spatio-Temporal, ON cells','Spatio-Temporal, OFF cells');
set(h,'FontSize',10)
axis equal

figure('Color','w');
plot((xx(on,2).^2)./(xx(on,1).^2),(xx(on,3).^2)./(xx(on,1).^2),'r+');
hold on;
plot((xx(off,2).^2)./(xx(off,1).^2),(xx(off,3).^2)./(xx(off,1).^2),'r*');
hold on;
plot([0,1],[0,1],'g')
axis equal
xlim([0,1]);
ylim([0,1]);
xlabel('Spatial');
ylabel('Spatio-Temporal')
h=legend('ON Cells','Off Cells','equality line');
set(h,'FontSize',10)
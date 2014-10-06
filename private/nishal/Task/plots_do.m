f=figure;
load('2_0pt4.mat')
plot(DD_log,II_log,'r','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')
% hold on


load('3_0pt4.mat')
plot(DD_log,II_log,'g','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')
% hold on


load('4_0pt4.mat')
plot(DD_log,II_log,'b','Marker','o')
% hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')

legend('2','3','4');

xlabel('Distortion');
ylabel('Rate');
title('discrimination 0.4')
print(f,'Discrimination p0=0pt4.png','-dpng')
%% 


f=figure;
load('5_2_recog.mat')
plot(DD_log,II_log,'r','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r','Marker','o')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g','Marker','o')
% hold on

load('3_2_recog.mat')
plot(DD_log,II_log,'g','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r','Marker','o')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g','Marker','o')
% hold on

load('5_3_recog.mat')
plot(DD_log,II_log,'b','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r','Marker','o')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g','Marker','o')
% hold on

legend('Q=5,P=2, p_0 = 0.1667','Q=3,P=2,p_0 = 0.4286','Q=5,P=3,p_0 =0.1364')
ylim([0,4]);

xlabel('Distortion');
ylabel('Rate');
title('Recognition with : (p_0)/(1-p_0)=1/3')

print(f,'Recognition.png','-dpng')
%%
f=figure;
load('2_0pt4.mat')
plot(DD_log,II_log,'r','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')
% hold on


load('2_0pt3.mat')
plot(DD_log,II_log,'g','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')
% hold on

load('2_0pt2.mat')
plot(DD_log,II_log,'b','Marker','o')
hold on
% plot(DD_log,UB_log+s_log.*DD_log,'r')
% hold on
% plot(DD_log,LB_log+s_log.*DD_log,'g')
% hold on

legend('0.4','0.3','0.2')

xlabel('Distortion');
ylabel('Rate');
title('discrimination 2')

print(f,'Discrimination 2.png','-dpng')


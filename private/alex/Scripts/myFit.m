function res=myFit(m,br,pol,bound)

m=m*pol;

[s,k]=max(m);

lower_a=0;
upper_a=0.01;
lower_b=40;
upper_b=200;
lower_c=0;
upper_c=100;


fit_res=fit((br:br+bound)',m,'gauss1',...
    'Lower',[lower_a, lower_b, lower_c],...
    'Upper',[upper_a, upper_b, upper_c],...
    'Startpoint',[s, k, 30]);
% resFit(1)=pol*fit_res.a1;
% resFit(2)=fit_res.b1;
% resFit(3)=fit_res.c1;
% 
% x=(br:br+bound)';
% y=resFit(1)*exp(-((x-resFit(2))/resFit(3)).^2);
% plot(x,m,'b','LineWidth',2)
% hold on
% plot(x,y,'r','LineWidth',2)
res=fit_res.b1;
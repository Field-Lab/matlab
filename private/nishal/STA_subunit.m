clear all

gam=[-10:0.01:10]';
std=1;
pgam = (1/(2*pi*std))*(exp(-(gam.^2)/(2*std^2)));
plot(gam,pgam)

%f= @(x) double(x>0).*(0.2*x).^2;%exp(0.5*x);
f= @(x) exp(0.2*x);
N= @(x) exp(x);

orig_rat=[];
calc_rat=[];
figure
for fact_n2=[0:0.05:1]

n2=5;
n1=n2*fact_n2;


eta1 = pgam'*N(n1*f(gam))
eta2 = pgam'*N(n2*f(gam))
sig1 = pgam'*(gam.*N(n1*f(gam)))
sig2 = pgam'*(gam.*N(n2*f(gam)))

fact1=(eta2*sig1)
fact2=(eta1*sig2)

fact1/fact2

orig_rat=[orig_rat;n1/n2];
calc_rat=[calc_rat;fact1/fact2];

extimated=vpa(fact1*sym('k1')+fact2*sym('k2'))

plot(N(n1*f(gam)),pgam,'r')
hold on
plot(N(n2*f(gam)),pgam,'b')
hold off
pause(1/120)
end


figure;
plot(calc_rat,orig_rat,'LineWidth',1.5);
hold on
plot(orig_rat,orig_rat,'r','LineWidth',1.5);
xlim([0,1]);
ylim([0,1]);
ylabel('O');
xlabel('$$\hat{O}$$','Interpreter','Latex');
title('O v/s $$\hat{O}$$','Interpreter','Latex');
legend('Plot','y=x line')


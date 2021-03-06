gam=[-10:0.01:10]';
std=1;
pgam = (1/(2*pi*std))*(exp(-(gam.^2)/(2*std^2)));
plot(gam,pgam)

%f= @(x) double(x>0).*(0.2*x).^1;%double(x>0).*(0.1*x).^2;%exp(0.5*x);
f= @(x) exp(0.05*x);
N= @(x) exp(x);

orig_rat=[];
calc_rat=[];
figure
for fact_n2=1:0.5:50
n1=0.1;
n2=n1*fact_n2;

eta1 = pgam'*N(n1*f(gam))
eta2 = pgam'*N(n2*f(gam))
sig1 = pgam'*(gam.*N(n1*f(gam)))
sig2 = pgam'*(gam.*N(n2*f(gam)))

fact1=(eta2*sig1)
fact2=(eta1*sig2)

fact1/fact2

orig_rat=[orig_rat;n2/n1];
calc_rat=[calc_rat;fact2/fact1];

extimated=vpa(fact1*sym('k1')+fact2*sym('k2'))

plot(N(n1*f(gam)),pgam,'r')
hold on
plot(N(n2*f(gam)),pgam,'b')
hold off
pause(1/120)
end

figure;
plot(orig_rat,calc_rat,'LineWidth',1.5);
hold on
plot(orig_rat,orig_rat,'r','LineWidth',1.5);
xlim([1,60]);
ylim([1,60]);
xlabel('O');
ylabel('$$\hat{O}$$','Interpreter','Latex');
title('$$\hat{O}$$ v/s O','Interpreter','Latex');
legend('Plot','y=x line')


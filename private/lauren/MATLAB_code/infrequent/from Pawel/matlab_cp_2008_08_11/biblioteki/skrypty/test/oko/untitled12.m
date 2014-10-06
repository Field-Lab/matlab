clear;

t=[1:1000];

A=0.7;
f=0.004;
fi=3;
offset=1.5;

szum=0.6*rand(1,1000)-0.3;
x=A*sin(2*pi*f*t+fi)+offset+szum;

[b,p]=nielfit(t,x,'sinus_fit',[0.8 0.0037 3 0],[0.001 0.00001 0.002 0.01]);

b
figure(1);
for i=1:4
    subplot(2,2,i);
    plot(p(:,i));
end

figure(2);
plot(t,x,t,feval('sinus_fit',b,t));
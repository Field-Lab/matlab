t=[1:32];

A=0.4;
f=0.037;
fi=10;
offset=-0.2;

b=[A f fi offset];
%szum=0.6*rand(1,1000)-0.3;
%x=A*sin(2*pi*f*t+fi)+offset+szum;
x=w(27,:);

[b,p]=nielfit(t,x,'sinus_fit',[A f fi offset],[0.0005 0.0002 0.0005 0.0002]);

b
figure(1);
for i=1:4
    subplot(2,2,i);
    plot(p(:,i));
end

figure(2);
s=feval('sinus_fit',b,t);
plot(t,x,t,feval('sinus_fit',b,t));
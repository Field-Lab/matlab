clear;
n=0.7;
b=n-1;

x=[0.1:0.1:8];
%amplituda, przy ktorej ekstremum wypadnie po czasie x: 
A=((x+2).^b-2*(x+1).^b+x.^b)./((x+3).^b-(x+2).^b-(x+1).^b+x.^b);
%wartosc napiecia w ekstremum w zaleznosci od czasu, po jakim pojawia sie
%ekstremum:
umax=A.*(x+3).^n-(1+A).*(x+2).^n+(2-A).*(x+1).^n-(1-A).*x.^n;
%dla danej wartosci opoznienia ekstremum, kiedy artefakt wchodzi w
%przedzial okreslony granicami "wartosc w ekstremum i minus wartosc
%ekstremum"
t=[0.01:0.01:5];
for i=1:length(A)
    a=A(i);
    u=a*(t+3).^n-(1+a).*(t+2).^n+(2-a).*(t+1).^n-(1-a).*t.^n;    
    
    b=find(abs(u)<=abs(umax(i)));    
    t0(i)=t(b(1));
end
figure(1)
plot(u)
figure(2)
subplot(2,2,1);
plot(x,A);
grid on;
subplot(2,2,2);
plot(x,umax);
grid on;
subplot(2,2,3);
plot(abs(umax),t0);
grid on;

clear;
n=0.8;
b=n-1;
x=[0.1:0.1:15];
%amplituda, przy ktorej ekstremum wypadnie po czasie x: 
A=((x+2).^b-(x+1).^b)./((x+1).^b-x.^b);
%wartosc napiecia w ekstremum w zaleznosci od czasu, po jakim pojawia sie
%ekstremum:
umax=-(x+2).^n+(1+A).*(x+1).^n-A.*x.^n;
%dla danej wartosci opoznienia ekstremum, kiedy artefakt wchodzi w
%przedzial okreslony granicami "wartosc w ekstremum i minus wartosc
%ekstremum"
t=[0.01:0.01:15];
for i=1:length(A)
    a=A(i);
    u=-(x+2).^n+(1+a).*((x+1).^n)-a.*x.^n;
    %figure(5)
    %plot(u);
    %grid on;
    %pause(0.3);

    b=find(abs(u)<=abs(umax(i)));    
    t0(i)=t(b(1));
end
figure(3)
subplot(2,2,1);
plot(x,A);
grid on;
subplot(2,2,2);
plot(x,umax);
grid on;
subplot(2,2,3);
plot(abs(umax),t0);
grid on;

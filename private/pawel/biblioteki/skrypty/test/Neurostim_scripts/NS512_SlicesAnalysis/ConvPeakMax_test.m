x=[0:0.05:30];

A=0.5;
tau=10;
sigma=1;
s1=A*(exp(-(x-tau).^2/(2*sigma^2)));

A=10;
tau=20;
sigma=0.05;
s2=A*(exp(-(x-tau).^2/(2*sigma^2)));

A=12;
tau=20.5;
sigma=0.05;
s3=A*(exp(-(x-tau).^2/(2*sigma^2)));

w=ones(1,10);

s=s1+s2+rand(1,length(s1))/3;
s=s1.*(0.4+rand(1,length(x))*1.2)+s2.*(0.4+rand(1,length(x))*1.2)+s3+rand(1,length(x))/3;

[c,smax]=ConvPeakMax(s,w);
smax
figure(11)
clf

subplot(2,1,1)
plot(s,'bd-')
grid on
%hold on
subplot(2,1,2)
plot(c/length(w))
grid on
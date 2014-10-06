fs=100000;

I0=10;
t0=0.001;

time=t0*fs*18;
t=t0*fs;
marg=t*2;

v_stim=zeros(1,time);
v_stim(1,marg+1:marg+t)=-I0;
v_stim(1,marg+t+1:marg+11*t)=I0/10;

figure(1);
clf;
plot(v_stim);
h=gca;
h=axis([0 1700 -12 2]);
set(h,'XTickLabel','   ');
h=text(100,-10,'-I_0');
set(h,'FontSize',14);
h=text(150,0.7,'I_0/10');
set(h,'FontSize',14);

fid=fopen('D:\Adapt_experiment\Stimuli_files\Tom_st.txt','r') 
a=fscanf(fid,'%f');
fclose(fid);
figure(4)
plot(0:0.0125:32-0.0125,a)
title('Full stimulus Tom')
 
%chirp
b=[(0:0.0125:32-0.0125)',a];
b=b(b(:,1)>=10.05&b(:,1)<=18,:);
b=b(:,2);
tt=0:0.0125:7.95;
f0=0.05;
k=1;
x=0:1/60:8;
ampl=30;
mean_shift=30;
y=(sin(2*pi*(f0+k/2*x).*x))*ampl+mean_shift;
figure(2)
hold off
plot(x,round(y))
hold on
plot(tt,b,'r')
title('Chirp')
legend('Function','Tom')
 
%resonance
b=[(0:0.0125:32-0.0125)',a];
b=b(b(:,1)>=20.12&b(:,1)<=28,:);
b=b(:,2);
tt=0:0.0125:7.875;
x=0:1/60:8;
ampl=30;
mean_shift=30;
y=(0.127*x.*cos(4*pi*x))*ampl+mean_shift;
yi = interp1(0:1/60:8,y,tt);%interpolated
figure(3)
hold off
plot(x,round(y))
hold on
plot(tt,b,'r')
title('Resonance')
legend('Function','Tom')
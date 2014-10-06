clear;

v0=1;
v1=0;

clock=5e6; %clock dla data_bus
Ts=50e-6; %us
N=clock*Ts; %ilosc stanow na kazdej linii data_bus dla jednego kompletu probek

%1. Okreslenie sekwencji bitow DLA JEDNEGO KANALU
td=1; %us

time_neg=2; %w jednostkach Ts
time_pos=2;
time_zero=0;

time_disch=6;

time_start=1;

time_full=15;
tf=time_full;

dac_value=[1 1 1 0 1 1 1]';

connect=ones(1,tf).*v0;
sw1=connect;
sw2=connect;
sw3=connect;
polarity=connect;
dac=ones(7,tf).*v0;

connect(1,(time_start+2):(time_start+1+time_neg+time_pos+time_zero))=v1;
sw1(1,(time_start+1):(time_start+1+time_neg+time_pos+time_zero+time_disch))=v1;
sw2(1,(time_start+time_neg+time_pos+time_zero+2):(time_start+time_neg+time_pos+time_zero+time_disch+1))=v1;
sw3=sw2;
polarity((time_start+1):(time_start+1+time_neg))=v1;
for i=(time_start+1):(time_start+1+time_neg+time_pos+time_zero+1)
	dac(:,i)=dac_value;    
end

figure(1);
%t=[1:time_full];
subplot(7,1,1);
plot(connect);
ylabel('connect');
subplot(7,1,2);
plot(sw1);
ylabel('sw1');
subplot(7,1,3);
plot(sw2);
ylabel('sw2');
subplot(7,1,4);
plot(sw3);
ylabel('sw3');
subplot(7,1,5);
plot(polarity);
ylabel('polarity');
subplot(7,1,6);
plot(dac(6,:));


%2. Wpisanie bitow na linie data_bus:
data_bus=ones(4,tf*N); %calkowita ilosc stanow dla calej stymulacji
channel=5;
for i=1:tf
    data_bus(4,(i-1)*N+(channel-1)*3+1)=connect(1,i);
    data_bus(3,(i-1)*N+(channel-1)*3+1)=sw1(1,i);
    data_bus(2,(i-1)*N+(channel-1)*3+1)=sw2(1,i);
    data_bus(1,(i-1)*N+(channel-1)*3+1)=sw3(1,i);
    
    data_bus(4,(i-1)*N+(channel-1)*3+2)=polarity(1,i);
    data_bus(3,(i-1)*N+(channel-1)*3+2)=dac(7,i);
    data_bus(2,(i-1)*N+(channel-1)*3+2)=dac(6,i);
    data_bus(1,(i-1)*N+(channel-1)*3+2)=dac(5,i);
    
    data_bus(4,(i-1)*N+(channel-1)*3+3)=dac(4,i);
    data_bus(3,(i-1)*N+(channel-1)*3+3)=dac(3,i);
    data_bus(2,(i-1)*N+(channel-1)*3+3)=dac(2,i);
    data_bus(1,(i-1)*N+(channel-1)*3+3)=dac(1,i);
end
size(data_bus)
t0=[0:tf*N-1]/(tf*N)*tf*Ts;
t=[0:2e-8:max(t0)];
data_bus=round(interp1(t0,data_bus',t,'linear'))';
%t=[1:td/Ts:tf];

figure(2)
for i=1:4
subplot(4,1,i)
plot(t,data_bus(i,:),'bd-')
end
figure(3)
plot(t,data_bus(1,:),'bd-',t,data_bus(2,:)+0.01,'gd-',t,data_bus(3,:)+0.02,'kd-',t,data_bus(4,:)+0.03,'rd-')
legend('db0','db1','db2','db3');

cd H:\pliki\nauka\stymulacja\chip\symulacje\pliki_dla_spice;
fid = fopen('data_bus1.txt','w');
fprintf(fid,'%12.10f %12.8f\n',[t; data_bus(1,:)]);
fclose(fid);

fid = fopen('data_bus2.txt','w');
fprintf(fid,'%12.10f %12.8f\n',[t; data_bus(2,:)]);
fclose(fid);

fid = fopen('data_bus3.txt','w');
fprintf(fid,'%12.10f %12.8f\n',[t; data_bus(3,:)]);
fclose(fid);

fid = fopen('data_bus4.txt','w');
fprintf(fid,'%12.10f %12.8f\n',[t; data_bus(4,:)]);
fclose(fid);
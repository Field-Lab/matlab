R=9200;
Vss=-2.5
Vdd=2.5

%cd /home/pawel/pliki/dane/;
cd C:\pliki\189e\pliki\dane;

s1c1=importdata('stim1_chip1.dat')/R;
s2c1=importdata('stim2_chip1.dat')/R;
s1c2=importdata('stim1_chip2.dat')/R;
s2c2=importdata('stim2_chip2.dat')/R;

a=size(s1c1);

vin=[-1:0.005:1];

figure(1);
for i=1:32
    subplot(4,8,i);
    %plot(vin,s1c1(i,:),vin,s2c1(i,a(2):-1:1));
    plot(vin,s1c1(i,:));
    axis([-1 1 -2e-4 2e-4]);
    grid on;
end

figure(2);
for i=1:32
    subplot(4,8,i);
    plot(vin,s1c2(i,:),vin,s2c2(i,a(2):-1:1));
    axis([-1 1 -2e-4 2e-4]);
    grid on;
end

figure(3);
for i=1:32
    subplot(4,8,i);
    plot(vin,s1c1(i,:),vin,s2c1(i,a(2):-1:1));
    axis([-.1 .1 -2e-5 2e-5]);
    grid on;
end

figure(4);
for i=1:32
    subplot(4,8,i);
    plot(vin,s1c2(i,:),vin,s2c2(i,a(2):-1:1));
    axis([-0.1 0.1 -2e-5 2e-5]);
    grid on;
end

dac=1;
figure(5);
subplot(2,1,1);
plot(vin,s1c1(dac+1,:),vin,s2c1(dac+1,a(2):-1:1));
grid on;
subplot(2,1,2);
plot(vin,s1c2(dac+1,:),vin,s2c2(dac+1,a(2):-1:1));
grid on;

figure(6);
for i=1:32
    subplot(4,8,i);
    plot(vin,s2c2(i,a(2):-1:1));
    %axis([-0.1 0.1 -2e-5 2e-5]);
    grid on;
end
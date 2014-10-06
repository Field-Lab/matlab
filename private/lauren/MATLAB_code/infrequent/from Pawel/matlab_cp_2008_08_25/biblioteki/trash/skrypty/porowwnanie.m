cd /mnt/win3/oko/2000-12-11/;

channels=[40 50 57 64];

figure(5);
t=[1: 15996];
for i=1:length(channels)
    subplot(2,2,i);
    r0=readconv('Data000cnst',header,nrchns,channels(i),samples);
    r1=readconv('Data001cnst',header,nrchns,channels(i),samples);
    r2=readconv('Data002cnst',header,nrchns,channels(i),samples);
    plot(t,r0,t,r1,t,r2);
    grid on;
end

    
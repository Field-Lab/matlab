cd /home/pawel/pliki/nauka/neuro64/listopad2002/3-11-2002;
cd /mnt/e-fat16/pawel/neuro64/listopad2002/4-11-2002;

%name='300.dat';
N=20000;
n=10;
fp=20000;
nchns=64; % ilosc kanalow
frqs=[20 30 40 50 100 300 500 1000 1500 2000 3000 4000 6000];

names(1)={'20.dat'};
names(2)={'30.dat'};
names(3)={'40.dat'};
names(4)={'50.dat'};
names(5)={'100.dat'};
names(6)={'300.dat'};
names(7)={'500.dat'};
names(9)={'1000.dat'};
names(1)={'1500.dat'};
names(10)={'2000.dat'};
names(11)={'3000.dat'};
names(12)={'4000.dat'};
names(12)={'6000.dat'};

tlumik=800; %orientacyjnie!
gen_ampl=200*0.001; %mV

a=importdata(name);
size(a)
data=zeros(n,nchns,N);
single=zeros(nchns,N);

for i=1:n % ze wzgledu na specyfike pliku
    a1=a(1,(i*nchns*N+1):(i+1)*nchns*N);
    for j=1:N
        single(1:nchns,j)=a1(1,((j-1)*nchns+1):j*nchns)';
    end    
    data(i,:,:)=single;
end

'data converted'
amplitudes=zeros(nchns,n);
characteristics=zeros(nchns,n,lengrh(frqs));

figure(1)
for i=1:nchns
    subplot(8,8,i)
    plot(single(i,:))
    grid on;
end

figure(2)
for i=1:nchns
    subplot(8,8,i)
    plot(abs(fft(single(i,:)))*2/N);
    axis([0 10000 0 1]);
    grid on;
end



for i=1:nchns
    for j=1:n
        ampl=max(abs(fft(data(j,i,:)))*2/N)/gen_ampl*tlumik;
        amplitudes(i,j)=ampl;
    end
end

figure(3)
plot(amplitudes);
grid on;

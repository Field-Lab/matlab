cd /home/pawel/pliki/nauka/neuro64/listopad2002/6-11-2002/ust4;
%cd /mnt/e-fat16/pawel/neuro64/listopad2002/4-11-2002/ust1;

%name='300.dat';
N=2000;
n=6;
fp=20000;

nchns=64; % ilosc kanalow
%frqs=[20 30 40 50 100 300 500 1000 1500 2000 3000 4000 6000];
frqs=[20 30 40 70 120 230 420 780 1020 1530 2020 4020 6020];

names(1)={'20.dat'};
names(2)={'30.dat'};
names(3)={'40.dat'};
names(4)={'70.dat'};
names(5)={'120.dat'};
names(6)={'230.dat'};
names(7)={'420.dat'};
names(8)={'780.dat'};
names(9)={'1020.dat'};
names(10)={'1530.dat'};
names(11)={'2020.dat'};
names(12)={'4020.dat'};
names(13)={'6020.dat'};


tlumik=800; %orientacyjnie!
%gen_ampl=200*0.001; %mV

data=zeros(n,nchns,N);
single=zeros(nchns,N);

amplitudes=zeros(nchns,n);
%characteristics=zeros(nchns,n,length(frqs));


for fr=1:length(frqs)
    
    name=names{fr}
    a=importdata(name);

for i=1:n % ze wzgledu na specyfike pliku
    a1=a(1,(i*nchns*N+1):(i+1)*nchns*N);
    for j=1:N
        single(1:nchns,j)=a1(1,((j-1)*nchns+1):j*nchns)';
    end    
    data(i,:,:)=single;
end

'data converted'

for i=1:nchns
    for j=1:n
        %ampl=max(abs(fft(data(j,i,:)))*2/N)/gen_ampl*tlumik;
        ampl=abs(fft(data(j,i,:)))*2/N;
        indeks=round(frqs(fr)*N/fp)+1;
        amplitudes(i,j)=ampl(1,indeks);
    end
end

characteristics(:,:,fr)=amplitudes;

end



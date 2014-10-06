cd /home/pawel/pliki/nauka/neuro64/listopad2002/4-11-2002/ust1;
name='300.dat';
czest=300;
fp=20000;
N=2000;

amplitudes=zeros(nchns,n);

a=importdata(name);

for i=1:n % ze wzgledu na specyfike pliku
    a1=a(1,(i*nchns*N+1):(i+1)*nchns*N);
    for j=1:N
        single(1:nchns,j)=a1(1,((j-1)*nchns+1):j*nchns)';
    end    
    data(i,:,:)=single;
end

figure(5);

for i=1:64
    signal=single(i,:);
    f=fft(signal);
    subplot(8,8,i);
    plot(abs(f));
end
    


for i=1:nchns
    for j=1:n
        %ampl=max(abs(fft(data(j,i,:)))*2/N)/gen_ampl*tlumik;
        ampl=abs(fft(data(j,i,:)))*2/N/gen_ampl*tlumik;
        indeks=round(czest*N/fp)+1
        amplitudes(i,j)=ampl(1,indeks);
    end
end

amplitudes;
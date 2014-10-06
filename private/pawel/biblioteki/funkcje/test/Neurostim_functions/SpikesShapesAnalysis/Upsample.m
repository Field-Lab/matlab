function y=Upsample(Waveform,N);
%Waveform = defined at x0=1,2,3,....
for i=1:length(Waveform)-1
    y(1+(i-1)*N)=Waveform(i);
    y([(i-1)*N+2:i*N+1])=Waveform(i)+(Waveform(i+1)-Waveform(i))*[1:N]/N;
end
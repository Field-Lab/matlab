function y=Quantization(x,qTime,qAmplitude,Dither);

Window=qTime;
for i=1:floor(length(x)/Window)
    x1(i)=mean(x([1:Window]+(i-1)*Window));
    %x([1:Window]+(i-1)*Window)=mean(x([1:Window]+(i-1)*Window));
end

if Dither
    x1=x1+(rand(1,length(x1))-0.5)*qAmplitude;
end
y0=round(x1/qAmplitude)*qAmplitude;

for i=1:length(y0)
    y([1:Window]+(i-1)*Window)=y0(i);
end
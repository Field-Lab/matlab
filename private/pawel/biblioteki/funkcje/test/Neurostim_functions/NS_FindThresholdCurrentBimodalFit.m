function ThresholdCurrent=NS_FindThresholdCurrentBimodalFit(Amplitudes,Efficacies,ThresholdEff);
%The function gives the values of current for which the effiacy is 50%. It
%uses linear interpolation. If the Efficacy vs Amplitude curve is not
%monotonic, there may be more than one output value!

f = @(p,x) 100 ./ (1 + exp(-(x-p(1))/p(2)));
p=nlinfit(Amplitudes,Efficacies,f,[0.5 0.1]);
N=100;
Amin=round(min(Amplitudes)*N)/N;
Amax=round(max(Amplitudes)*N)/N;

x=[Amin-0.1:1/N:Amax+0.1];
y=f(p,x);
figure(103)
plot(x,y,Amplitudes,Efficacies,'rd-');
grid on;
SF=sign(y-ThresholdEff);
SFN=find(SF==-1);
SFP=find(SF==1);
ThresholdCurrent=(x(max(SFN))+x(min(SFP)))/2;
Amp=0.26;
step=0.024;

for i=1:14
    Amp=Amp*1.1;
    A(i)=Amp;
    
    Amp_q=round(Amp/step);
    Aq(i)=Amp_q;
    
    Amp2=Amp*3/4;
    Amp2_q=round(Amp2/step);
    A2q(i)=Amp2_q;
    
    factor(i)=Amp2_q/Amp_q;
end
Aq=Aq*step;
A2q=A2q*step;
plot(A,factor,'bd-');
grid on;
function dydt = myIAF(t,y)
%this function solves system of ODE for neuron model
%parameters for ode
global Vth
tauM=30; %ms
eL=-65; %mV
Vreset=-65; %mV
Rm=90; %MOhm


%external stimulation
if t<50
    Iapp=0;
elseif t<150
    Iapp=0.2;
else Iapp=0;
end

%ODE
if y<Vth
    dydt = (eL-y+Rm*Iapp)/tauM;
else display(y)
    dydt = Vreset;
    display('Ok')
    display(t)
end
end
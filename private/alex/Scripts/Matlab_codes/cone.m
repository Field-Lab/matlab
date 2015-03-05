function cone
%this function defines conditions of differential equations, calls
%computational function and creates voltage plots
tic
tspan = [0 5000];
y0=[-45.92; 0.069975; 0.390990; 0.98679; 0.0162183; 0.600613; 0; 0; 21.9556; 21.9556];
% Solve the problem using ode45 (based on Runge-Kutta, medium accuracy)
solution=ode45(@f,tspan,y0);
plot(solution.x,solution.y(1,:)')
% fprintf('%f\n',solution.y(1,:))
a=toc
end


function dydt = f(t,y)
%this function solves system of ODE for neuron model
%parameters for ode

cm = 3;
hbar = 0.217;
eh = -32;
eKv = -80;
Kvbar = 2.0;
eCa = 38.8; 
Cabar = 7.99;
aoCa = 0.0031;
VhalfCa = -16.6;
Sca = 5.7;
eCl = -45;
Clbar = 6.77;
SCl = 0.09;
Clh = 0.37;
FactorCaI = 0.45;
eKca = -80;
Kcabar = 0.4;
lbar = 0.011;
el = 0;
tauR = 2.97;
tauE = 39.84;
Cbeta = 0.015;
Kbeta = 0.0004;
ac = 0.065;
nc = 2;
tauC = 26.16;
nx = 1;
scale = 1;
current = 10000000;

if t<500
    Ihat=0;
elseif t<600
    Ihat=current;
else Ihat=0;
end


%functions

%h channel
Ih=hbar*(1-(1+3*y(2))*(1-y(2)^3))*(y(1)+eh);
alphah=0.001*18/(exp((y(1)+88)/12)+1);
betah=0.001*18/(exp(-(y(1)+18)/19)+1);

% Ikv channel
Ikv=Kvbar*y(3)^3*y(4)*(y(1)-eKv);
alphamKv=0.001*5*(100-y(1))/(exp(100-y(1)/42)-1);
betamKv=0.0001*9*exp((20-y(1))/40);
alphahKv=0.001*0.15*exp(-y(1)/22);
betahKv=0.001*0.4125/(exp((10-y(1))/7)+1);
taumKv=1/(alphamKv+betamKv);
infmKv=alphamKv/(alphamKv+betamKv);
tauhKv=1/(alphahKv+betahKv);
infhKv=alphahKv/(alphahKv+betahKv);

% Ca channel
Ica=Cabar*y(5)*(y(1)-eCa);
alphaCa=aoCa*exp((y(1)-VhalfCa)/(2*Sca));
betaCa=aoCa*exp(-(y(1)-VhalfCa)/(2*Sca));
infCa=alphaCa/(alphaCa+betaCa);
tauCa=1/(alphaCa+betaCa);

% Cl channel
Cas=-0.2+FactorCaI*(-Ica)*1.8821875e50;
mCl=1/(1+exp((Clh-Cas)/SCl));
Icl=Clbar*mCl*(y(1)-eCl);

% Kca channel
Ikca=Kcabar*y(6)^2*(y(1)-eKca);
alphamKca=0.001*15*(80-y(1))/(exp((80-y(1))/40)-1);
betamKca=0.001*20*exp(-y(1)/35);
infmKca=alphamKca/(alphamKca+betamKca);
taumKca=1/(alphamKca+betamKca);

% leak
Il=lbar*(y(1)-el);

% Outer current
beta=Cbeta+Kbeta*y(8);
alpha=1/(1+(ac*y(9))^nc);
Ios=y(10)^nx;
Iouter=(21.9615-Ios)*scale;

%ODE system - mostly gating variables
dydt = [
    % potential
    -1/cm*(Ih+Ikv+Ica+Icl+Ikca+Il+Iouter); % (1)
    
    %h channel
    (alphah/(alphah+betah)-y(2))/1/(alphah+betah); % (2)
    
    % Ikv channel
    (infmKv-y(3))/taumKv; %(3)
    (infhKv-y(4))/tauhKv; %(4)    
    
    % Ca channel
    (infCa-y(5))/tauCa; % (5)
    
    % Kca channel
    (infmKca-y(6))/taumKca; % (6)
    
    % Outer current
    (Ihat-y(7))/tauR; %(7)
    (y(7)-y(8))/tauE; % (8)
    (Ios-y(9))/tauC; % (9)
    alpha-y(10)*beta; % (10)
    ];

end
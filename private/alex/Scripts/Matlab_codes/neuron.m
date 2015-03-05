function neuron
%this function defines conditions of differential equations, calls
%computational function and creates voltage plots

tspan = [0 200];
y0 = [-70.038; 0.8522; 0.000208; 0.2686; 0.5016;-70.038; 0.8522; 0.000208; -70.038; 0.8522; 0.000208; -70.038; 0.8522; 0.000208];
% Solve the problem using ode45 (based on Runge-Kutta, medium accuracy)
[t,solution]=ode45(@f,tspan,y0);
plot(t,solution)
end

function y=gammaf(VV,theta,sigma)
%this function serves computation of the gamma function, called multiply
y=1.0/(1.0+exp(-(VV-theta)/sigma));
end

function dydt = f(t,y)
%this function solves system of ODE for neuron model
%parameters for ode
g=0.34;
g1=0.34;
gA=0.19;
theta_m=-24.0;
gNa=112.5;
gK=225.0;
gL=0.25;
sigma_m=11.5;
theta_h=-58.3;
sigma_h=-6.7;
theta_n=-12.4;
sigma_n=6.8;
theta_t_h=-60;
sigma_t_h=-12.0;
theta_tna=-14.6;
sigma_tna=-8.6;
theta_tnb=1.3;
sigma_tnb=18.7;
theta_a=-50;
sigma_a=20;
theta_b=-70;
sigma_b=-6;
tau_b=150;
tau_a=2;
V_Na=50.0;
V_K=-90.0;
V_L=-70.0;

%functions
%soma
m_inf=gammaf(y(1),theta_m,sigma_m);
h_inf=gammaf(y(1),theta_h,sigma_h);
n_inf=gammaf(y(1),theta_n,sigma_n);
a_inf=gammaf(y(1),theta_a,sigma_a);
b_inf=gammaf(y(1),theta_b,sigma_b);
tau_h=0.5+14.0*gammaf(y(1),theta_t_h,sigma_t_h);
tau_n=(0.087+11.4*gammaf(y(1),theta_tna,sigma_tna))*(0.087+11.4*gammaf(y(1),theta_tnb,sigma_tnb));
i_na=gNa*m_inf^3*y(2)*(y(1)-V_Na);
i_k=gK*(y(3)^2)*(y(1)-V_K);
i_l=gL*(y(1)-V_L);
i_a=gA*y(4)^3*y(5)*(y(1)-V_K);

%basal dendrite
m_infd=gammaf(y(6),theta_m,sigma_m);
h_infd=gammaf(y(6),theta_h,sigma_h);
n_infd=gammaf(y(6),theta_n,sigma_n);
tau_hd=0.5+14.0*gammaf(y(6),theta_t_h,sigma_t_h);
tau_nd=(0.087+11.4*gammaf(y(6),theta_tna,sigma_tna))*(0.087+11.4*gammaf(y(6),theta_tnb,sigma_tnb));
i_nad=gNa*m_infd^3*y(7)*(y(6)-V_Na);
i_kd=gK*(y(8)^2)*(y(6)-V_K);
i_ld=gL*(y(6)-V_L);

%apical dendrite 1
m_infd1=gammaf(y(9),theta_m,sigma_m);
h_infd1=gammaf(y(9),theta_h,sigma_h);
n_infd1=gammaf(y(9),theta_n,sigma_n);
tau_hd1=0.5+14.0*gammaf(y(9),theta_t_h,sigma_t_h);
tau_nd1=(0.087+11.4*gammaf(y(9),theta_tna,sigma_tna))*(0.087+11.4*gammaf(y(9),theta_tnb,sigma_tnb));
i_nad1=gNa*m_infd1^3*y(10)*(y(9)-V_Na);
i_kd1=gK*(y(11)^2)*(y(9)-V_K);
i_ld1=gL*(y(9)-V_L);

%apical dendrite 2
m_infd2=gammaf(y(12),theta_m,sigma_m);
h_infd2=gammaf(y(12),theta_h,sigma_h);
n_infd2=gammaf(y(12),theta_n,sigma_n);
tau_hd2=0.5+14.0*gammaf(y(12),theta_t_h,sigma_t_h);
tau_nd2=(0.087+11.4*gammaf(y(12),theta_tna,sigma_tna))*(0.087+11.4*gammaf(y(12),theta_tnb,sigma_tnb));
i_nad2=gNa*m_infd2^3*y(13)*(y(12)-V_Na);
i_kd2=gK*(y(14)^2)*(y(12)-V_K);
i_ld2=gL*(y(12)-V_L);

%external stimulation
if t<50
    Iapp=0;
elseif t<150
    Iapp=200;
else Iapp=0;
end

%define stimulation application
Iapp1=Iapp;
Iapp2=0;
Iappd=0;
Iapps=0;

%ODE system
dydt = [
    %soma
    -i_na-i_k-i_a-i_l-g*(y(1)-y(6))+Iapps %(1)
    (h_inf-y(2))/tau_h %(2)
    (n_inf-y(3))/tau_n %(3)
    (a_inf-y(4))/tau_a %(4)
    (b_inf-y(5))/tau_b %(5)
    % basal dendrite
    -i_nad-i_kd-i_ld-g*(y(6)-y(1))-g1*(y(6)-y(9))-g1*(y(6)-y(12))+Iappd %(6)
    (h_infd-y(7))/tau_hd %(7)
    (n_infd-y(8))/tau_nd %(8)
    % apical dendrite 1
    -i_nad1-i_kd1-i_ld1-g1*(y(9)-y(6))+Iapp1 %(9)
    (h_infd1-y(10))/tau_hd1 %(10)
    (n_infd1-y(11))/tau_nd1 %(11)
    % apical dendrite 2
    -i_nad2-i_kd2-i_ld2-g1*(y(12)-y(6))+Iapp2 %(12)
    (h_infd2-y(13))/tau_hd2 %(13)
    (n_infd2-y(14))/tau_nd2 %(14)
    ];

end
VgFilePath='C:\pawel\nauka\Neurostim-3\Vgp.csv';
N=64;

a=importdata('Vgp.csv');
InputCurrent=a.data(:,1);
Vg=a.data(:,2);

figure(1)
plot(InputCurrent,Vg,'bd-')
grid on

W=0.4; % in microns
L=20; % in microns

I=1e-6; %na wyjsciu DACa jest 4 x wiecej, na jednym tranzystorze jest 16 x mniej

SigmaV=15.2/sqrt(W*L)*1e-3;

Vg1=interp1(InputCurrent,Vg,I);
Iplus=interp1(Vg,InputCurrent,Vg1+SigmaV);
Iminus=interp1(Vg,InputCurrent,Vg1-SigmaV);

CurrentError=max(abs(Iplus-I),abs(Iminus-I))/I*100
MaxCurrentError=CurrentError*sqrt(N)
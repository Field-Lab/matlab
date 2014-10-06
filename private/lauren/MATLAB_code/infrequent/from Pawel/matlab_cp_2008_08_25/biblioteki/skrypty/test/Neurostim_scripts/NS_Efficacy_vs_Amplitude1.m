%od 21
AmpTri=[0.27 0.3 0.32 0.37 0.4 0.45 0.49 0.53 0.6 0.67 0.75 0.8 0.89];
AmpBi=AmpTri;
N=98;
Fs=14;
Channel=17;

%Channel14:
Neuron1Tri=[0 0 0 0 0 0 50 95 98 98 98 98 98]./N*100;
Neuron1Bi=[0 0 7 44 90 98 98 98 98 98 98 98 98]./N*100;

Neuron2Tri=[0 0 0 0 0 0 0 0 0 0 27 68 98]./N*100;
Neuron2Bi=[0 0 0 0 0 0 0 0 23 81 98 98 98]./N*100;

%Channel15:
Neuron1Tri=[2 6 10 31 37 57 62 92 98 98 98 98 98];
Neuron1Bi=[1 6 6 17 31 54 75 95 98 98 98 98 98];

%Channel17:
Neuron1Tri=[0 0 0 0 0 0 0 1 6 13 52 95 98];
Neuron1Bi=[0 0 4 5 13 27 55 88 98 98 98 98 98];


figure(1)
plot(AmpTri,Neuron1Tri,'b*-',AmpBi,Neuron1Bi,'k*-')
grid on;
h=gca;
set(h,'FontSize',Fs);
h=xlabel('Amplitude [uA]');
set(h,'FontSize',Fs);
h=ylabel('Eficcacy [%]');
set(h,'FontSize',Fs);
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\processed\' num2str(Channel) 'neuron1 Eff_vs_Amp'];
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
print(hj, '-dtiff', name);

figure(2)
plot(AmpTri,Neuron2Tri,'b*-',AmpBi,Neuron2Bi,'k*-')
grid on;
h=gca;
set(h,'FontSize',Fs);
h=xlabel('Amplitude [uA]');
set(h,'FontSize',Fs);
h=ylabel('Eficcacy [%]');
set(h,'FontSize',Fs);
name=['E:\analysis\2008-03-21-0\2008_05_05\report\subtr\100us_tri\processed\' num2str(Channel) 'neuron2 Eff_vs_Amp'];
hj=gcf;
set(hj, 'PaperOrientation', 'portrait');
print(hj, '-dtiff', name);

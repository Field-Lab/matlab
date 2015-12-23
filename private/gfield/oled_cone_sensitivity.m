cd /snle/lab/Experiments/Calibration/Data/2011-06-17/
load CorrectedRGBMatrix

% plot OLED spectra
figure(1); clf;
subplot(2,1,1)
wavelengths = [380:2:778];
plot(wavelengths, CorrectedRGBMatrix(:,1), 'r')
hold on
plot(wavelengths, CorrectedRGBMatrix(:,2), 'g')
plot(wavelengths, CorrectedRGBMatrix(:,3), 'b')
title('OLED spectra')
hold off


% compute the three cone spectra
% Baylor photoreceptor spectral nomogram 
% constants for a polynomial fit
WaveLengthRange = 778:-2:380; % in nm
WL = WaveLengthRange.*0.001;
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;

% L Cone
lambdaMax = 561;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
LConePhotonSensitivity = 10.^LogPhotonSensitivity;

% M Cone
lambdaMax = 531;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
MConePhotonSensitivity = 10.^LogPhotonSensitivity;

% S Cone
lambdaMax = 430;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
SConePhotonSensitivity = 10.^LogPhotonSensitivity;


% plot cone spectra
subplot(2,1,2)
plot(WaveLengthRange, LConePhotonSensitivity, 'r')
hold on
plot(WaveLengthRange, MConePhotonSensitivity, 'g')
plot(WaveLengthRange, SConePhotonSensitivity, 'b')
title('cone spectra')


cone_spectra = [LConePhotonSensitivity(200:-1:1); MConePhotonSensitivity(200:-1:1); SConePhotonSensitivity(200:-1:1)];

sensitivity_matrix = cone_spectra * CorrectedRGBMatrix





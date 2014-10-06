%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START HERE

% Set physical constants
SpeedOfLight = 3e8;
PlanksConstant = 6.626e-34;
WaveLengthRange = 778:-2:380; % in nm
OpsinPhotosensitivity = 9.6e-21; % m^2 
wavelengths = 400:1:800;

%% UDT
%UDT Information
cd /Volumes/lab/Experiments/Calibration/2013-07-19/ 
load cal_udt_spectrum;


%% OLED
Powers = [11.82e-9 19.32e-9 11.11e-9];
spot_area = 0.003 * 0.00225; % meters (aperture size) "optics from above"
cd /Volumes/lab/Experiments/Calibration/2013-07-19/ 
load oled_rgb_matrix
oled_rgb_matrix(:,2) = oled_rgb_matrix(:,2)./norm(oled_rgb_matrix(:,2));
oled_rgb_matrix(:,3) = oled_rgb_matrix(:,3)./norm(oled_rgb_matrix(:,3));
oled_rgb_matrix(:,4) = oled_rgb_matrix(:,4)./norm(oled_rgb_matrix(:,4));

RGBMatrix = oled_rgb_matrix;


%% LED
Powers = 34.94e-9; % Watts: Power made with 
spot_area = 0.011^2 * pi;
cd /Volumes/lab/Experiments/Calibration/2013-07-19/ 
load LED_spectrum
LED_spectrum = LED_spectrum ./ norm(LED_spectrum);

RGBMatrix = LED_spectrum;


%% LED calibration

% Normalize and scale the RGBMAtrix
TruePowerScalers = Powers ./ (cal_udt_spectrum(:,2)' * RGBMatrix');
CalibratedRGBMatrix = TruePowerScalers .* RGBMatrix;

SummedMonitorSpectra = CalibratedRGBMatrix'; % this is to make it compatible with a multi-primary source

%Sanity Check: CheckPower should equal the sum of Powers
Powers
CheckPower = dot(CalibratedRGBMatrix, cal_udt_spectrum(:,2))

%% OLED calibration

TruePowerScalers = Powers ./ (cal_udt_spectrum(:,2)' * RGBMatrix(:,2:4));
TruePowerScalers = repmat(TruePowerScalers, 401, 1); 
CalibratedRGBMatrix = TruePowerScalers .* RGBMatrix(:,2:4);

SummedMonitorSpectra = sum(CalibratedRGBMatrix,2);

%Sanity Check: CheckPower should equal the sum of Powers
SummedPowers = sum(Powers)
CheckPower = dot(SummedMonitorSpectra,  cal_udt_spectrum(:,2))


%% Rod calculation begins here:
% Rod spectral sensitivty fit with Baylor nomogram
WL = 800:-1:400;
WL = WL.*0.001;
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;
lambdaMax = 491;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
RodPhotonSensitivity = 10.^LogPhotonSensitivity;

WavelengthsCorrected = wavelengths * 1e-9; % nm -> m
Intensity = (SummedMonitorSpectra .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
PhotonFlux = Intensity ./ spot_area;
AbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux,RodPhotonSensitivity([401:-1:1])');

%RodCollectingArea = 1.2e-12; % MONKEY (Baylor 1984, Table 1) m^2
RodCollectingArea = 0.5e-12; % MOUSE m^2
EffectivePhotonFlux = dot(PhotonFlux,RodPhotonSensitivity([401:-1:1])');
PhotonCatchRate = EffectivePhotonFlux * RodCollectingArea

%% compute based on time
Time = 10 * 1;
UnBleachedPigment = exp(-AbsorptionRate * Time)
 
% compute based on amount of unbleached pigment
%UnBleachedPigment = 0.99;
%Time = -1*log(UnbleachedPigment) ./ AbsorptionRate
 
% calculate bleaching rate
TotalRodPigment = 1.4e8; % pigment content is thought to be ~1e8 to 3e8
BleachingRate = TotalRodPigment * (1-UnBleachedPigment) / Time;





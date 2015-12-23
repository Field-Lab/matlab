%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% START HERE

% Set physical constants
SpeedOfLight = 3e8;
PlanksConstant = 6.626e-34;
WaveLengthRange = 778:-2:380; % in nm
OpsinPhotosensitivity = 9.6e-21; % m^2 

Rig = 'A'; % A or B
LightPath = 'inverted'; % inverted of conventional
Array = '1501'; % array specification only used for inverted optics
Objective = '6.5x'; % objective specification only used for inverted optics
NDF = 0; 
display_type = 'oled'; % supported 'oled' or 'crt'

% Rig B Bleaching Parameters and monitor spectra
if strcmp(Rig, 'B')
    if strcmp(display_type, 'oled')
        Powers = [35.5e-9 59.3e-9 40.3e-9] ./ 2; % 6.5x in Watts
        %Powers = [35.55e-9 59.39e-9 40.36e-9]; % in Watts, for 10x objective 
        %Powers = [27.70e-9 46.26e-9 31.44e-9]; % in Watts, for 4x objective
        Powers = Powers ./ (10^NDF);
        spot_area = 0.0016 * 0.002133; % 6.5x 600px x 800 px, rectangular display
        %spot_area = 0.00111 * 0.00148; % 10x 600px x 800 px, rectangular display
        %spot_area = 0.0028 * 0.00373; % 4x600px x 800 px, rectangular display
        cd /snle/lab/Experiments/Calibration/Data/2011-06-17/
        load CorrectedRGBMatrix
    elseif strcmp(display_type, 'crt')

        if strcmp(LightPath, 'conventional')
            Powers = [3.55e-9 4.28e-9 4.67e-9]; % in Watts conventional optics
            AdjustedNDF = NDF-0.6;
            Powers = Powers / (10^AdjustedNDF);
            SpotRadius = 0.002800/2; % meters (aperture size) "optics from above"
            % spot size calibration in 2008-08-19
            cd /snle/lab/Experiments/Calibration/Data/2006-04-23/ % rig B conventional
            load RGBMatrix
            CorrectedRGBMatrix = RGBMatrix(1:200,:);
        elseif strcmp(LightPath, 'inverted')
            if strcmp(Objective, '10x')
                if strcmp(Array, '1501')
                    Powers = [((3.15e-9)*2.34) ((3.33e-9)*2.39) ((2.40e-9)*2.45)];
                    SpotRadius = 0.00168/2; %10X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/10x/Array1501/ 
                    load CorrectedRGBMatrix
                else
                    Error('only array 1501 has been calibrated with the 10x objective')
                    return
                end
            elseif strcmp(Objective, '6.5x')
                if strcmp(Array, '504') 
                    Powers = [4.73e-9 4.97e-9 2.91e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.00228/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/6.5x/Array504 
                    load CorrectedRGBMatrix
                elseif strcmp(Array, '508')
                    Powers = [2.61e-9 3.11e-9 3.00e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.00228/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/6.5x/Array508
                    load CorrectedRGBMatrix
                elseif strcmp(Array, '1501')
                    Powers = [2.98e-9 3.27e-9 2.78e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.00228/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/6.5x/Array1501
                    load CorrectedRGBMatrix
                end
            end
        end
    else
        error('display type not supported, must be oled, or crt')
    end
% Rig A
elseif strcmp(Rig, 'A')
    if strcmp(display_type, 'oled')
        Powers = [15.64e-9 25.4e-9 16.14e-9] ./ 2; % 6.5x in Watts
        %Powers = [21.12e-9 34.29e-9 21.79e-9]; % 10x in Watts
        Powers = Powers ./ (10^NDF);
        %spot_area = 0.00102 * 0.00136; % 6.5x, 600px x 800 px, rectangular display
        spot_area = 0.00075 * 0.001; % 10x, 600px x 800 px, rectangular display
        %cd /snle/lab/Experiments/Calibration/Data/2011-06-17/
        cd ~/Desktop/2011-06-17/
        load CorrectedRGBMatrix
    elseif strcmp(display_type, 'conventional')
        
        if strcmp(LightPath, 'conventional')
            Powers = [4.51e-9 6.87e-9 6.43e-9]; % in Watts conventional optics
            AdjustedNDF = NDF-0.6;
            Powers = Powers / (10^AdjustedNDF);
            SpotRadius = 0.003300/2; % meters (aperture size) "optics from above"
            cd /snle/lab/Experiments/Calibration/Data/2006-04-22/ % rig A conventional
            load RGBMatrix
            CorrectedRGBMatrix = RGBMatrix(1:200,:);
        elseif strcmp(LightPath, 'inverted')
            if strcmp(Objective, '10x')
                Error('there is no 10x objective calibrated for rig A')
                return
            elseif strcmp(Objective, '6.5x')
                if strcmp(Array, '508') 
                    Powers = [3.25e-9 5.48e-9 4.88e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.0022/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigA/6.5x/Array508 
                    load CorrectedRGBMatrix
                elseif strcmp(Array, '503')
                    Powers = [3.84e-9 5.25e-9 3.49e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.0022/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigA/6.5x/Array503
                    load CorrectedRGBMatrix
                elseif strcmp(Array, '1501')
                    Powers = [3.9e-9 5.46e-9 4.18e-9]; % in Watts inverted optics
                    Powers = Powers ./ (10^NDF);
                    SpotRadius = 0.0022/2; %6.5X objective from below
                    cd /snle/lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigA/6.5x/Array1501
                    load CorrectedRGBMatrix
                end
            end
        end
    elseif strcmp(Rig, 'C')
        if strcmp(display_type, 'conventional')
            if strcmp(LightPath, 'conventional')
                Powers = [12.81e-9 24.2e-9 22.0e-9]; % in Watts conventional optics
                AdjustedNDF = NDF;
                Powers = Powers / (10^AdjustedNDF);
                SpotRadius = 0.003300/2; % meters (aperture size) "optics from above"
                cd /snle/lab/Experiments/Calibration/Data/2011-09-16/ % rig C conventional
                load RGBMatrix
                CorrectedRGBMatrix = RGBMatrix(1:200,:);
            end
        else
            error('display type not supported, must be oled, or crt')
        end
    end
end
      
RGBMatrix = CorrectedRGBMatrix;

%UDT Information
%cd /snle/lab/Experiments/Calibration/Data/2005-03-26/UDT/
cd ~/Desktop/UDT/
UDT = load('UDT.2');
Wavelengths = 380:2:778;
Index = find(Wavelengths == 500);
CalibratedUDT = UDT ./ UDT(Index); 

% Normalize and scale the RGBMAtrix
TruePowerScalers = Powers ./ (CalibratedUDT' * RGBMatrix);
TruePowerScalers = repmat(TruePowerScalers, 200, 1); 
CalibratedRGBMatrix = TruePowerScalers .* RGBMatrix;

%check influence of UDT spectra
%TruePowerScalers = Powers ./ (ones(200,1)' * RGBMatrix);
%CalibratedRGBMatrix = repmat(TruePowerScalers, 200, 1) .* RGBMatrix; 

SummedMonitorSpectra = sum(CalibratedRGBMatrix,2);

%CheckPower should equal the sum of Powers
SummedPowers = sum(Powers)
CheckPower = dot(SummedMonitorSpectra, CalibratedUDT)


% Rod spectral sensitivty fit with Baylor nomogram
WL = 778:-2:380;
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


%%%%%%%%
% actualy calculation begins here:

% polychromatic source
if ~exist('spot_area')
    spot_area = (pi*SpotRadius.^2);
end

WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (SummedMonitorSpectra .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
PhotonFlux = Intensity ./ spot_area;
AbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux,RodPhotonSensitivity([200:-1:1])');

RodCollectingArea = 1.2e-12; %(Baylor 1984, Table 1) m^2
EffectivePhotonFlux = dot(PhotonFlux,RodPhotonSensitivity([200:-1:1])');
PhotonCatchRate = EffectivePhotonFlux * RodCollectingArea

%plot(SummedMonitorSpectra ./ max(SummedMonitorSpectra))
%hold on
%plot(PhotonFlux ./ max(PhotonFlux), 'r')
%hold off

% compute based on time
Time = 10 * 1;
UnBleachedPigment = exp(-AbsorptionRate * Time)
 
% compute based on amount of unbleached pigment
%UnBleachedPigment = 0.99;
%Time = -1*log(UnbleachedPigment) ./ AbsorptionRate
 

% calculate bleaching rate
TotalRodPigment = 1.4e8; % pigment content is thought to be ~1e8 to 3e8
BleachingRate = TotalRodPigment * (1-UnBleachedPigment) / Time




%%%%%
% relative sensitivity of rhodopsin to each monitor spectra
UnitConversionMatrix = repmat(WavelengthsCorrected ./ (PlanksConstant * SpeedOfLight), 3, 1); % converts from power to photons
RelativeScales = RodPhotonSensitivity([200:-1:1]) * (CalibratedRGBMatrix .* UnitConversionMatrix')
RodRelativeScales = RelativeScales ./ max(RelativeScales)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cone Bleaching Rates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the three cone spectra
% Baylor photoreceptor spectral nomogram 
% constants for a polynomial fit
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
%correct for pigment self-screening
LConePhotonSensitivity = 10.^LogPhotonSensitivity;
%LConePhotonSensitivity = 1-(10.^(-0.17 .* LConePhotonSensitivity));
%LConePhotonSensitivity = LConePhotonSensitivity ./ max(LConePhotonSensitivity);
%LConeEnergySensitivity = LConePhotonSensitivity .* WaveLengthRange; % convert to energy units
%LConeEnergySensitivity = LConeEnergySensitivity ./ max(LConeEnergySensitivity);  % normalize

% M Cone
lambdaMax = 531;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
MConePhotonSensitivity = 10.^LogPhotonSensitivity;
%MConePhotonSensitivity = 1-(10.^(-0.17 .* MConePhotonSensitivity));
%MConePhotonSensitivity = MConePhotonSensitivity ./ max(MConePhotonSensitivity);
%MConeEnergySensitivity = MConePhotonSensitivity .* WaveLengthRange; % convert to energy units
%MConeEnergySensitivity = MConeEnergySensitivity ./ max(MConeEnergySensitivity); % normalize

% S Cone
lambdaMax = 430;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
SConePhotonSensitivity = 10.^LogPhotonSensitivity;
%SConePhotonSensitivity = 1-(10.^(-5 .* SConePhotonSensitivity));
%SConePhotonSensitivity = SConePhotonSensitivity ./ max(SConePhotonSensitivity);
%SConeEnergySensitivity = SConePhotonSensitivity .* WaveLengthRange; % convert to energy units
%SConeEnergySensitivity = SConeEnergySensitivity ./ max(SConeEnergySensitivity); % 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ConeCollectingArea = 0.37e-12;
LConeCatch = dot(PhotonFlux, LConePhotonSensitivity([200:-1:1])') * ConeCollectingArea
MConeCatch = dot(PhotonFlux, MConePhotonSensitivity([200:-1:1])') * ConeCollectingArea
SConeCatch = dot(PhotonFlux, SConePhotonSensitivity([200:-1:1])') * ConeCollectingArea

% polychromatic source
LConeAbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux, LConePhotonSensitivity([200:-1:1])');
MConeAbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux, MConePhotonSensitivity([200:-1:1])');
SConeAbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux, SConePhotonSensitivity([200:-1:1])');
 
% compute based on time
Time = 10;
LConeUnBleachedPigment = exp(-LConeAbsorptionRate * Time);
MConeUnBleachedPigment = exp(-MConeAbsorptionRate * Time);
SConeUnBleachedPigment = exp(-SConeAbsorptionRate * Time);
 
% calculate bleaching rate
TotalConePigment = 5.5e7; % pigment content is thought to be ~1e8 to 3e8 from Rods
LConeBleachingRate = TotalConePigment * (1-LConeUnBleachedPigment) / Time
MConeBleachingRate = TotalConePigment * (1-MConeUnBleachedPigment) / Time
SConeBleachingRate = TotalConePigment * (1-SConeUnBleachedPigment) / Time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L cone fluxes broken down by monitor primary

% red source on L cone
RGBColor = 1; % red
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
RedPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
LConeAbsorptionRate = OpsinPhotosensitivity * dot(RedPhotonFlux, LConePhotonSensitivity([200:-1:1])');
LConeUnBleachedPigment = exp(-LConeAbsorptionRate * Time);
LConeBleachingRateRed = TotalConePigment * (1-LConeUnBleachedPigment) / Time

% green source on L cone
RGBColor = 2; % green
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
GreenPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
LConeAbsorptionRate = OpsinPhotosensitivity * dot(GreenPhotonFlux, LConePhotonSensitivity([200:-1:1])');
LConeUnBleachedPigment = exp(-LConeAbsorptionRate * Time);
LConeBleachingRateGreen = TotalConePigment * (1-LConeUnBleachedPigment) / Time


% blue source on L cone
RGBColor = 3; % blue
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
BluePhotonFlux = Intensity ./ (pi*SpotRadius.^2);
LConeAbsorptionRate = OpsinPhotosensitivity * dot(BluePhotonFlux, LConePhotonSensitivity([200:-1:1])');
LConeUnBleachedPigment = exp(-LConeAbsorptionRate * Time);
LConeBleachingRateBlue = TotalConePigment * (1-LConeUnBleachedPigment) / Time


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% M cone fluxes broken down by monitor primary

% red source on M cone
RGBColor = 1; % red
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
RedPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
MConeAbsorptionRate = OpsinPhotosensitivity * dot(RedPhotonFlux, MConePhotonSensitivity([200:-1:1])');
MConeUnBleachedPigment = exp(-MConeAbsorptionRate * Time);
MConeBleachingRateRed = TotalConePigment * (1-MConeUnBleachedPigment) / Time

% green source on M cone
RGBColor = 2; % green
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
GreenPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
MConeAbsorptionRate = OpsinPhotosensitivity * dot(GreenPhotonFlux, MConePhotonSensitivity([200:-1:1])');
MConeUnBleachedPigment = exp(-MConeAbsorptionRate * Time);
MConeBleachingRateGreen = TotalConePigment * (1-MConeUnBleachedPigment) / Time

% blue source on M cone
RGBColor = 3; % blue
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
BluePhotonFlux = Intensity ./ (pi*SpotRadius.^2);
MConeAbsorptionRate = OpsinPhotosensitivity * dot(BluePhotonFlux, MConePhotonSensitivity([200:-1:1])');
MConeUnBleachedPigment = exp(-MConeAbsorptionRate * Time);
MConeBleachingRateBlue = TotalConePigment * (1-MConeUnBleachedPigment) / Time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S cone fluxes broken down by monitor primary
% red source on M cone

% red source on S cone
RGBColor = 1; % red
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
RedPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
SConeAbsorptionRate = OpsinPhotosensitivity * dot(RedPhotonFlux, SConePhotonSensitivity([200:-1:1])');
SConeUnBleachedPigment = exp(-SConeAbsorptionRate * Time);
SConeBleachingRateRed = TotalConePigment * (1-SConeUnBleachedPigment) / Time

% green source on S cone
RGBColor = 2; % green
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
GreenPhotonFlux = Intensity ./ (pi*SpotRadius.^2);
SConeAbsorptionRate = OpsinPhotosensitivity * dot(GreenPhotonFlux, SConePhotonSensitivity([200:-1:1])');
SConeUnBleachedPigment = exp(-SConeAbsorptionRate * Time);
SConeBleachingRateGreen = TotalConePigment * (1-SConeUnBleachedPigment) / Time

% blue source on S cone
RGBColor = 3; % blue
WavelengthsCorrected = Wavelengths * 1e-9; % nm -> m
Intensity = (CalibratedRGBMatrix(:,RGBColor) .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
BluePhotonFlux = Intensity ./ (pi*SpotRadius.^2);
SConeAbsorptionRate = OpsinPhotosensitivity * dot(BluePhotonFlux, SConePhotonSensitivity([200:-1:1])');
SConeUnBleachedPigment = exp(-SConeAbsorptionRate * Time);
SConeBleachingRateBlue = TotalConePigment * (1-SConeUnBleachedPigment) / Time


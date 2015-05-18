function iso = isomerizations(power_reading, area, varargin)

% TO DO
% Fix units + get constants values
% check file loading + wavelengths of files

% UNITS
% area should be in um^2
% power should be in milliwatts

% TAKING THE MEASUREMENTS
% The power reading should be taken at Ch2 on the meter reading, which is
%   calibrated for 500nm. Meter sensitivity spectrum should be 1 at 500nm.
% When you put the power meter up, shift it around to make sure the reading
%   doesn't change. If the monitor is well centered on the meter, the
%   reading shouldn't change.


%% LOAD AND DEFINE STUFF

p = inputParser;
p.addParameter('meter_sensitivity_file', '/Volumes/Lab/Experiments/Calibration/Data/udt.1') % eg udt.1
p.addParameter('monitor_emission_file', '/Volumes/Lab/Experiments/Calibration/Data/rgb.1') % eg RGB.1
p.addParameter('cone_absorption_file', '/Volumes/Lab/Experiments/Calibration/Data/lms.1') % eg LMS.1
p.parse

% Constants
h = 1; % Planck's constant
c = 1; % speed of light

% Load up spectra
meter_sensitivity_spectrum = load(p.Results.meter_sensitivity_file);
monitor_emission_spectrum = load(p.Results.monitor_emission_file);
cone_absorption_spectrum = load(p.Results.cone_absorption_file);

% Check this or load it up or something?
wavelengths = 370:1:730;

%% CALCULATE THE POWER EMITTED FROM THE MONITOR

% We already have the spectrum in the rgb file, and it is is accurate 
% up to a scale factor alpha. 
% power reading = < meter sensitivity | alpha * monitor emission > so
alpha = power_reading/(meter_sensitivity_spectrum'*monitor_emission_spectrum);

% The true emission spectrum then is 
power_spectrum = alpha * monitor_emission_spectrum; 

%% CALCULATE THE FLUX: photons per second absorbed by the cones
% Flux in photons / second is given by
% Flux = sum over wavelengths of Power(wavelength)*Wavelength / (h*c)
flux = (power_spectrum .* wavelengths)' * cone_absorption_spectrum / (h*c);

%% CALCULATE THE INTENSITY: photons per second per area
% Intensity in photons/(area*second) is given by
% Intensity = flux / area
intensity = flux/area;

%% CALCULATE ISOMERIZATIONS: photons absorbed per cone per second
% Isomerizations in photons/(cone * sec)
% iso = intensity * area/molecule * molecules/cone
% Area/molecule * molecules/cone ~ 1 um^2/cone
iso = intensity;

end
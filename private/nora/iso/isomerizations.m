function iso = isomerizations(power_reading, area, rig, meter, varargin)
% This function takes in the power reading, area, rig, and meter,
% and spits out the isomerizations for L, M, S cones, and Rods, in that
% order. It relies on the subroutine spectra_loading and access to the
% Volumes.

% INPUTS

% Required
% power_reading: in milliwatts either [w], [r g b], or [r g b w]
% area in um^2 microns^2
% rig: currently the only monitor emission file coded in is for Rig 4.
% meter: which UDT you did the measurements with, as a string. EG '18707'

% Optional
% meter_sensitivity_file: use a file other than the default meter file, which is ['/Volumes/Lab/Experiments/Calibration/Data/UDT/UDT-' meter '.1']
% meter_wavelengths: if you change the file, you should input the
%   wavelengths for the new file, default is 381:2:779;
% monitor_emission_file: use a file other than the rig default file
% monitor_wavelengths: if you change the file, input the new wavelengths.
%   Default is 370:730;

% UNITS
% area should be in um^2
% power should be in milliwatts

% TAKING THE MEASUREMENTS
% Get a microscope slide stage. They look like CDs and should be in the
%   optical drawer of the toolbox in a plastic bag.
% Get the microscope slide ruler. It is in a wooden box in the measuring
%   drawer of the toolbox. Put the slide on the stage on the microscope. 
% Look at the stimulus from above. Get the stimulus and the slide ruler in
%   focus in the same plane. Measure stim size. The units on this ruler are
%   100 um per division
% Get a power meter. Put the head of the meter on the microscope stage. 
%   Turn it on, Ch2, Power reading, try 10^-2 to start. Write down which
%   UDT you are using.
% Turn on 50% white. 
% Move the power meter around until you are confident that the whole
%   stimulus is being captured. When you put the power meter up, 
%   shift it around to make sure the reading
%   doesn't change. If the monitor is well centered on the meter, the
%   reading shouldn't change. Then tape it down.
% Record the power for 50% white, and separately for 50% red, blue, and
% 	green. mW is standard for the UDT, so you can just read the number off the
%   screen and the setting (eg 10^-2)
% 
% Your input to the code should look something like 
% iso = isomerizations(5 * 10^(-2), 24*25*100*100, 4, '18707');

% Constants
h = 6.62607*10^-31; % Planck's constant in mW * seconds^2
c = 2.998*10^14; % microns per second

%% LOAD UP VARIOUS SPECTRA
p = inputParser;
p.addParameter('meter_sensitivity_file', ['/Volumes/Lab/Experiments/Calibration/Data/UDT/UDT-' meter '.1']) % eg udt.1
p.addParameter('monitor_emission_file', 0) % eg RGB.1
p.addParameter('monitor_wavelengths', 381:2:779)
p.addParameter('meter_wavelengths', 370:730)
p.parse

% Find the monitor emission for the rig
if ~p.Results.monitor_emission_file
    if rig==4
        monitor_emission_file = '/Volumes/Lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/6.5x/Array504/CorrectedRGBMatrix.mat';
    else
        error('Must input a monitor emission file for your rig')
    end
else
    monitor_emission_file = p.Results.monitor_emission_file;
end

% Load up the spectra
spectra = spectra_loading(monitor_emission_file, p);

%% CALCULATE THE POWER EMITTED FROM THE MONITOR

% We already have the spectrum in the rgb file, and it is is accurate 
% up to a scale factor alpha. 
% power reading = < meter sensitivity | alpha * monitor emission > 
% so alpha = power_reading / < meter sensitivity | monitor emission > 
% and power_spectrum = alpha | monitor emission >

% If you took RGB measurements, it calculates alpha for each color gun, and
% finds the true spectrum by adding all three colors times their alphas
if length(power_reading)>1
    alpha = power_reading(1:3)./(spectra.meter*spectra.monitor);
    power_spectrum = alpha(1:3) * spectra.monitor';
end

% If you took white and RGB measurements, it compares them here, but
% ultimately uses the RGB measurements
if length(power_reading) == 4
    alpha_white = power_reading(4)./(spectra.meter*sum(spectra.monitor, 2)); % take the sum because used white light
    power_spectrum_white = alpha_white * sum(spectra.monitor, 2)';
    plot(spectra.wavelengths, power_spectrum)
    hold on
    plot(spectra.wavelengths, power_spectrum_white)
    legend('Alpha calc separately for each channel', 'Alpha calc from white')
    text(400, 0.06, ['RGB:  ' num2str(alpha)])
    text(400, 0.065, ['W:  ' num2str(alpha_white)])
    text(400, 0.07, 'Alpha')
    title('White versus RGB Readings')
end

% If you ONLY took a white measurement, it just calculates one alpha, and
% assumes all three colors add with the same weight
if length(power_reading) == 1
    alpha_white = power_reading./(spectra.meter*sum(spectra.monitor, 2)); % take the sum because used white light
    power_spectrum = alpha_white * sum(spectra.monitor, 2)';
end

%% CALCULATE THE FLUX: photons per second absorbed by the photoreceptors
% Flux in photons / second is given by
% Flux = sum over wavelengths of Power(wavelength)*Wavelength / (h*c) for
% each of the three colors and for rods
wavelengths_um = spectra.wavelengths*10^-3; % nanometers to micrometers
flux = 10^-3*(power_spectrum .* wavelengths_um) * spectra.PR / (h*c); % 10^-3 because thats the wavelength spacing for the integral

%% CALCULATE THE INTENSITY: photons per second per area
% Intensity in photons/(area*second) is given by
% Intensity = flux / area
intensity = flux/area;

%% CALCULATE ISOMERIZATIONS: photons absorbed per photoreceptor per second
% Isomerizations in photons/(cone * sec)
% iso = intensity * area/molecule * molecules/cone
% Area/molecule * molecules/cone ~ 1 um^2/cone
iso = intensity;

end

function spectra = spectra_loading(monitor_emission_file,p)
% This function finds the files you need, loads them in, and matches the
%   wavelengths across the different spectra. 
% The cone file is the same always, the monitor depends on the rig, and the
%   meter depends on the meter
% 
% INPUTS
% Required
% rig: currently the only monitor emission file coded in is for Rig 4.
% meter: which UDT you did the measurements with, as a string. EG '18707' 

%% Define common wavelengths
monitor_wavelengths = p.Results.monitor_wavelengths;
meter_wavelengths = p.Results.meter_wavelengths; 
PR_wavelengths = 370:730; % PR = photoreceptor
% Set the common range
spectra.wavelengths = max([min(meter_wavelengths) min(monitor_wavelengths) min(PR_wavelengths)]):min([max(meter_wavelengths) max(monitor_wavelengths) max(PR_wavelengths)]);

%% Load files and check that they are the right size
load('macaque_photoreceptor_spectra.mat') % cone
PR_spectra = photoreceptor.spectra; % PR = photoreceptor
if length(PR_spectra) ~= length(PR_wavelengths)
    error('Cone wavelengths do not match the cone spectrum length')
end
meter_spectrum = dlmread(p.Results.meter_sensitivity_file, '\n', 1, 0); % meter 
if length(meter_spectrum) ~= length(meter_wavelengths)
    meter_spectrum = meter_spectrum(meter_spectrum ~= 0);
    if length(meter_spectrum) ~= length(meter_wavelengths)
        error('Meter wavelengths do not match the meter spectrum length')
    end
end
if strcmp(monitor_emission_file((end-2):end), 'mat')
    temp = who('-file', monitor_emission_file);
    load(monitor_emission_file); % monitor
    eval(['monitor_spectra =' temp{1} ';'])
else
    try
        monitor_spectra = dlmread(monitor_emission_file, '\n', 1, 0);
        warn('The monitor spectrum was not loaded from a mat file so it might be wrong or include a header')
    catch
        error('Cannot read monitor file')
    end
end
if length(monitor_spectra) ~= length(monitor_wavelengths)
    error('Monitor wavelengths do not match the monitor spectrum length')
end

%% Get files over common wavelengths
spectra.monitor = interp1(monitor_wavelengths, monitor_spectra, spectra.wavelengths);
spectra.PR = interp1(PR_wavelengths, PR_spectra, spectra.wavelengths); % PR = photoreceptor
spectra.meter = interp1(meter_wavelengths, meter_spectrum, spectra.wavelengths);
spectra.meter = spectra.meter/ spectra.meter(spectra.wavelengths == 500); % normalize meter reading

end
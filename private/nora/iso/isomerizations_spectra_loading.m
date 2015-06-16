function spectra = isomerizations_spectra_loading(rig,meter,varargin)
% This function finds the files you need, loads them in, and matches the
%   wavelengths across the different spectra. 
% The cone file is the same always, the monitor depends on the rig, and the
%   meter depends on the meter
% 
% INPUTS
% Required
% rig: currently the only monitor emission file coded in is for Rig 4.
% meter: which UDT you did the measurements with, as a string. EG '18707' 

% TAKING THE MEASUREMENTS
% The power reading should be taken at Ch2 on the meter reading, which is
%   calibrated for 500nm. Meter sensitivity spectrum should be 1 at 500nm.
% When you put the power meter up, shift it around to make sure the reading
%   doesn't change. If the monitor is well centered on the meter, the
%   reading shouldn't change.

% QUESTIONS
% cone absportion units?
% should be 201 but is only 200... why

%% Get filenames etc
p = inputParser;
p.addParameter('meter_sensitivity_file', ['/Volumes/Lab/Experiments/Calibration/Data/UDT/UDT-' meter '.1']) % eg udt.1
p.addParameter('monitor_emission_file', 0) % eg RGB.1
p.addParameter('monitor_wavelengths', 381:2:779)
p.addParameter('meter_wavelengths', 370:730)
p.parse
if ~p.Results.monitor_emission_file
    if rig==4
        monitor_emission_file = '/Volumes/Lab/Experiments/Calibration/Data/2008-08-19/Matlab-files/RigB/6.5x/Array504/CorrectedRGBMatrix.mat';
    else
        error('Must input a monitor emission file for your rig')
    end
else
    monitor_emission_file = p.Results.monitor_emission_file;
end

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
function iso = isomerizations(power_reading, area, spectra)
% This function takes in the spectra spit out by
% isomerizations_spectra_loading.m, as well as the power reading and area
% measurements, and spits out the isomerizations for L, M, S cones, and Rods, in that
% order

% INPUTS
% Required
% power_reading: [r g b w] white is optional in milliwatts
% area in um^2 microns^2

% UNITS
% area should be in um^2
% power should be in milliwatts

% TAKING THE MEASUREMENTS
% The power reading should be taken at Ch2 on the meter reading, which is
%   calibrated for 500nm. Meter sensitivity spectrum should be 1 at 500nm.
% When you put the power meter up, shift it around to make sure the reading
%   doesn't change. If the monitor is well centered on the meter, the
%   reading shouldn't change.

% QUESTIONS
% cone absportion units?
% still not sure about the unit scaling for the first part??
% should be 201 but is only 200... why

% Constants
h = 6.62607*10^-31; % Planck's constant in mW * seconds^2
c = 2.998*10^14; % microns per second

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
wavelengths_um = spectra.wavelengths*10^-3;
flux = (power_spectrum .* wavelengths_um) * spectra.PR / (h*c); % units here??

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
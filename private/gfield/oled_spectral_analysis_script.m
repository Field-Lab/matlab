% set directory to location with calibration files
cd /snle/lab/Experiments/Calibration/Data/2011-06-17/

% load information about each primary and white
white_primary = dlmread('white.txt');
red_primary = dlmread('red.txt');
green_primary = dlmread('green.txt');
blue_primary = dlmread('blue.txt');

% normalize these to have unit magnitude;
red_primary = red_primary ./ norm(red_primary);
green_primary = green_primary ./ norm(green_primary);
blue_primary = blue_primary ./ norm(blue_primary);

% set wavelength range of calibration
wavelengths = 380:2:780;

% plot the spectra of the primaries
figure(1); clf; hold on;
plot(wavelengths, red_primary, 'r')
plot(wavelengths, green_primary, 'g')
plot(wavelengths, blue_primary, 'b')
hold off


CorrectedRGBMatrix(:,1) = red_primary(1:200);
CorrectedRGBMatrix(:,2) = green_primary(1:200);
CorrectedRGBMatrix(:,3) = blue_primary(1:200);


save CorrectedRGBMatrix CorrectedRGBMatrix


% AKHEITMAN 2013-12-18
% Reworking starting 2014-01-05
% Redocumenting for handover 2015-04-04
%%%%%%%%%%%%%%%%%%%%%%%


% Wavelengths in nanometers, unless otherwise specified
% integrating areas tracked explictly in microns squared

%% Load up scalar power measurements, Cone/UDT/CRT spectral profiles

restoredefaultpath
clear; close all; clc

% PARAMETERS
% Effective Cone absorption area: from the Rieke Lab
cone_effectivearea_micronsq =  0.6;  % in microns squared  
planck           = 6.625* (10^-34); % plancks constant
c_um             = 3*(10^8)*(10^6); % speed light in microns
c_nm             = 3*(10^8)*(10^9);
udt_veridical_nm = 500; % wavelength at which udt recordings veridical


% MEASUREMENTS OF SCALAR POWER FROM AKHEITMAN DECEMBER 2013
% CONSULT ScalarPower_nanoWatts_AKHeitman_2013December.pdf
RGB_recordedWatts_rigA = (10^-3)*(10^-6)*[3.8 5.5 4.3];  % 240 by 240 pixels 
RGB_recordedWatts_rigB = (10^-3)*(10^-6)*[4.0 4.7 3.6];  % 240 by240 pixels
stimarea_micronsq =  (240* 4.5) ^2; %
% Light from below, with Cylinder, 6.5x condensor, through the board
% Stimulus was 240 by 240 pixels [maxR, maxG, maxB]

% OLD SPECTRAL MEASUREMENTS
load spectra_CRTemission.mat
load spectra_UDTsensitivity.mat
load spectra_LMSenergyabsorption.mat


%% Unpack relative spectra, onto common wavelength domain
% xtick refers to wavelength in nm for which we recorded udt and CRTemission
udtxtick_nm          = spectra_UDTsensitivity.wavelengths_nm;
crtxtick_nm          = spectra_CRTemission.wavelengths_nm;
lmsxtick_nm          = spectra_LMSenergyabsorption.wavelengths_nm;

% SET COMMON WAVELENGTH DOMAIN
xtick_nm = max([min(udtxtick_nm), min(crtxtick_nm), min(lmsxtick_nm)])...
    :1:min([max(udtxtick_nm),max(crtxtick_nm),max(lmsxtick_nm)]);
% note: nonoverlapping domains of CRTemissions and udt_sensitivty very small


% INTERPOLATE CRT UDT LMS TO COMMON DOMAIN
RGB_relativespectralpower_rigA        = interp1(crtxtick_nm,spectra_CRTemission.RGBMatrix_RigA,xtick_nm);
RGB_relativespectralpower_rigB        = interp1(crtxtick_nm,spectra_CRTemission.RGBMatrix_RigB,xtick_nm);
LMS             = interp1(lmsxtick_nm, spectra_LMSenergyabsorption.vals, xtick_nm);
UDT_head_1      = interp1(udtxtick_nm,spectra_UDTsensitivity.head_1,xtick_nm);
UDT_head_2      = interp1(udtxtick_nm,spectra_UDTsensitivity.head_2,xtick_nm);
UDT_head_3      = interp1(udtxtick_nm,spectra_UDTsensitivity.head_3,xtick_nm);
UDT_avg         = (UDT_head_1)/sum(UDT_head_1) + (UDT_head_2)/sum(UDT_head_2) + (UDT_head_3)/sum(UDT_head_3); 
UDT_calibrated  = (1/UDT_avg(find(xtick_nm == udt_veridical_nm))) * UDT_avg;


% VERIFICATION PLOTS
figure;
subplot(2,2,1); xlim([350 750]); hold on;
plot(xtick_nm, fliplr(RGB_relativespectralpower_rigA)); 
xlabel('normed RigA emission'); hold off;
subplot(2,2,3); xlim([350 750]); hold on;
plot(xtick_nm, fliplr(RGB_relativespectralpower_rigB));
xlabel('normed RigB emission'); hold off;
subplot(2,2,2); xlim([350 750]); hold on;
for i_udt = 1:3
    switch i_udt
        case 1
            rawudt = UDT_head_1;
            plot(xtick_nm, rawudt, 'k'); hold on;
        case 2
            rawudt = UDT_head_2;
            plot(xtick_nm, rawudt, 'g'); 
        case 3
            rawudt = UDT_head_3;
            plot(xtick_nm, rawudt, 'm'); 
    end
end
xlabel('raw udt sensitivities'); hold off;
subplot(2,2,4); xlim([350 750]); hold on;
plot(xtick_nm,UDT_calibrated); 
xlabel('avg udt, veridical at 500'); hold off;

%% Actual Computation

for i_rig = 1:2
    if i_rig == 1
        RGB_recordedpower_Watts     = RGB_recordedWatts_rigA;
        RGB_relativespectralpower   = RGB_relativespectralpower_rigA; 
    elseif i_rig == 2
        RGB_recordedpower_Watts     = RGB_recordedWatts_rigB;
        RGB_relativespectralpower   = RGB_relativespectralpower_rigB;
    end
    RGB_conversion_constant = RGB_recordedpower_Watts ./ (UDT_calibrated* RGB_relativespectralpower);
    RGB_spectralpower_Watts = RGB_relativespectralpower * diag(RGB_conversion_constant);
    FullWhite_spectralpower_Watts = sum(RGB_spectralpower_Watts,2) ;

    [~,LMS_maxind]   = max(LMS);
    LMS_peakwaves_nm = xtick_nm(LMS_maxind);


    LMS_absorbedwatts             = (cone_effectivearea_micronsq/stimarea_micronsq) * RGB_spectralpower_Watts' * LMS; 
    LMS_photonsabsorbed_matrix    =    (1 /(planck*c_nm))  * LMS_absorbedwatts *  diag(LMS_peakwaves_nm);
    
    LMS_photonsabsorbed.red   = LMS_photonsabsorbed_matrix(1,:);
    LMS_photonsabsorbed.green = LMS_photonsabsorbed_matrix(2,:);
    LMS_photonsabsorbed.blue  = LMS_photonsabsorbed_matrix(3,:);
    LMS_photonsabsorbed.white = sum(LMS_photonsabsorbed_matrix);
    
    if i_rig == 1
        LMS_photonsabsorbed_rigA = LMS_photonsabsorbed
    elseif i_rig == 2
        LMS_photonsabsorbed_rigB = LMS_photonsabsorbed
    end
end






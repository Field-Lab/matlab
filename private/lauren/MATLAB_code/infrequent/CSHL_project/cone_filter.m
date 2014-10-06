function [cone_response cone_response_profile] = cone_filter(cone_pos, cone_type, stimulus_fg, stimulus_bg, fg_spec, bg_spec)

% cone_pos: x,y position (microns)
% fg_spec = power at wavelengths 380:5:780 nm (1x81)
% bg_spec = power spectrum of background
%
% stimulus_fg = matrix of values from 0 to 1 specifying scaling factor to apply to
% fg_spec
%
%

if ~(size(fg_spec,1) == 1)
    fg_spec = fg_spec';
end
if ~(size(bg_spec,1) == 1)
    bg_spec = bg_spec';
end


%% for testing
% clear all
% cone_pos = [20.5 50];
% cone_type = 'L';
% spotCenter = [30, 60];
% spotRadius = 20;
% stimSize = [300 250];
% 
% % make a achromatic stimulus: equal amount of r,g,b phosphors
% load B_monitor
% phosphors = B_monitor;
% %spectrum = 380:5:780;
% gray_spec = (phosphors(:,1) + phosphors(:,2) + phosphors(:,3))';
% 
% stimulus = generateSpotStimulus(spotCenter, spotRadius, stimSize);
% 
% stimulus_fg = psf_filter(stimulus, cone_type);
% stimulus_bg = psf_filter(~stimulus, cone_type);
% 
% fg_spec = (5/3)*gray_spec;
% bg_spec = gray_spec;


%% apply spectral cone filter to spectral stimulus (dot product of
%stimulus spectrum and sensitivity spectrum)

load T_cones_sp % first index corresponds to L,Mstimulus_spec_filtered,S cones

if strcmpi(cone_type, 'L')
    cone_response_fg = sum(T_cones_sp(1,:).*fg_spec);
    cone_response_bg = sum(T_cones_sp(1,:).*bg_spec);
elseif strcmpi(cone_type, 'M')
    cone_response_fg = sum(T_cones_sp(2,:).*fg_spec);
    cone_response_bg = sum(T_cones_sp(2,:).*bg_spec);
else
    error('cone type must be L or M')
end

%% apply foreground and background cone responses to spatial stimulus

stimulus_spec_filtered = cone_response_fg*stimulus_fg + cone_response_bg*stimulus_bg;

%%
% estimated cone inner segment diameter at 12 mm eccentricity = 8.5 microns
% estimate cone spatial filter as full-width at half-max = 61.5% of cone
% aperture = 5.2275

sigma = 5.2275/(2*sqrt(2*log(2))); %formula relating standard deviation to FWHM
covariance = [sigma, 0;
              0,     sigma];

remainder(1) = mod(cone_pos(1),1);
remainder(2) = mod(cone_pos(2),1);

%calculate multivariate normal for a subset of the space--need to make this
%more flexible
grid = zeros(21*21, 2);
for ii = 1:21; %x-position
    for jj = 1:21; %y-position
        grid((jj-1)*21+ii,1) = ii-11;
        grid((jj-1)*21+ii,2) = jj-11;
    end
end
filter_center = mvnpdf(grid, remainder, covariance); %centered on remainder
%value order: f(x1,y1), f(x2,y1),....f(x1,y2), f(x2,y2)....

filter_center = reshape(filter_center, 21, []);
%first index corresponds with x position, second index corresponds with y
%position (in microns, starting from 1)

filter = zeros(size(stimulus_spec_filtered));
filter(floor(cone_pos(1))-10:floor(cone_pos(1)+10), floor(cone_pos(2))...
    -10:floor(cone_pos(2)+10)) = filter_center;

%figure %calculated cone filter
%imagesc(filter)
%colormap(gray)

%% apply spatial cone filter to spectrally-filtered stimulus (dot product)

cone_response = sum(sum(stimulus_spec_filtered.*filter));
cone_response_profile = stimulus_spec_filtered.*filter; %contribution of each spatial location to cone response
% figure
% imagesc(filter)
% colormap(gray)



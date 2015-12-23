function cone_rgb = estimate_cone_rgb(datarun, cell_spec, cone_max_stixels, cone_locations, all_sig_stixels, Wc, varargin)
%
% Estimates the RGB signals at each cone location.  The function takes a
% datarun structure, cell_spec and several inputs related to cone finding.
% At each cone location the function calculates from the RF the weighted
% sum of RGB values.  The weighting is determined by the spatial provile of
% the cone -- handed to the function with the input Wc.
%
% usage: cone_rgb = estimate_cone_rgb(datarun, cell_spec, cone_max_stixels, cone_locations, all_sig_stixels, Wc, varargin)
%
% inputs:
%   datarun             standard datarun structure
%   cell_spec           specify which cells to used
%   cone_max_stixels    peak cone locations
%   cone_locations      optimized cone locations
%   all_sig_stixels     MxN matrix of summed RFs -- produced by local max
%                           cone finding
%   Wc                  PxQ Matrix of spatial cone Rfs where P = MxN and Q
%                           is the number of cones
%   
% optional inputs
%   rad         1       radius around the peak cone location for
%                       calculating the RGB sensitivity of the cone
% output
%   cone_rgb            Qx3 matrix of cones x RGB triplet.  These triplets
%                          are intended to be used to classify the cones
%
% GDF 2013-02
%

% parse inputs
p = inputParser;
p.addParamValue('rad', 1, @isnumeric);
p.parse(varargin{:});
rad = p.Results.rad;

% FUNCTION BEGINS HERE

% get cell indices and determine some numbers
rgc_indices = get_cell_indices(datarun, cell_spec);
num_cones = size(cone_locations,1);

% setup progress dots
dotter = floor(num_cones/20); % will put down 20 dots to indicate progress
fprintf('getting rgb values for each cone \n')

if strcmp(datarun.stimulus.independent,'nil')
    cone_rgb = ones(size(cone_locations, 1), 3);
else
    % initialize cone_rgb
    cone_rgb = zeros(size(cone_locations, 1), 3);
    % get rf info for each cone
    for cn = 1:num_cones

        if mod(cn, dotter) == 0
            fprintf('.')
        end

        % cone_location
        cone_peak = cone_max_stixels(cn,:);

        % get ROI around the cone

        % note the x range and y range of pixels within this square
        % take that to be the ROI
        xrng = max(round(cone_peak(1)-rad),1):min(round(cone_peak(1)+rad),datarun.stimulus.field_height);
        yrng = max(round(cone_peak(2)-rad),1):min(round(cone_peak(2)+rad),datarun.stimulus.field_width);    

        % translate these subscripts to a matrix index (coordinate pair ->
        % integer)
        cone_stixel = sub2ind([datarun.stimulus.field_width,datarun.stimulus.field_height], cone_peak(1),cone_peak(2));

        % identify rgcs the have a significant stixel at this cone center
        cells_of_interest = find(all_sig_stixels(cone_stixel,:) > 0);

        % note number of these rgcs
        num_coi_per_cone = length(cells_of_interest);

        % get rf patch around the cone from these rgcs
        cone_images = zeros(length(xrng), length(yrng), 3, num_coi_per_cone);
        for rgc = 1:num_coi_per_cone

            % get the correct index into datarun for the RGC
            rgc_index = rgc_indices(cells_of_interest(rgc));

            % get the spatial RF from the STA
            rf = get_rf(datarun, datarun.cell_ids(rgc_index));

            % normalized by robust_std to get into SNR units
            rf = rf ./ robust_std(reshape(rf, [],1))^2;

            % extract pixels over and around cone from RG
            cone_images(:,:,:,rgc) = rf(xrng,yrng,:);
        end

        % combine info across RGCs that sample from the cone
        cone_image = sum(cone_images,4);

        % reshape cone image for each color to a vector
        cone_image_r = reshape(cone_image(:,:,1), [],1);
        cone_image_g = reshape(cone_image(:,:,2), [],1);
        cone_image_b = reshape(cone_image(:,:,3), [],1);
        
        % Get the same cone image
        cone_fit = reshape(Wc(:,cn), datarun.stimulus.field_width, datarun.stimulus.field_height);
        cone_fit = cone_fit(xrng,yrng);
        cone_fit = reshape(cone_fit, 1,[]);
        
        % scale each stixel by the RF 
 
        % compute a spatial weighted sum the color signals
        red_sig = cone_fit * cone_image_r;
        green_sig = cone_fit* cone_image_g;
        blue_sig = cone_fit* cone_image_b;  
        
        % make the vector unit length.
        cone_rgb(cn,1:3) = [red_sig, green_sig, blue_sig] ./ norm([red_sig, green_sig, blue_sig]);
    end
end
fprintf('\n')
    
    
    
    
    
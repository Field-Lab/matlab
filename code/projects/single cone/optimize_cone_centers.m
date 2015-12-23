function new_centers = optimize_cone_centers(datarun, cell_spec, cone_locations, all_sig_stixels, varargin)
%
% The function finds the centers of mass (COMs) of cones given RFs and
% rough cone locations found from local-max calculation
%
% usage: new_centers = optimize_cone_centers(datarun, cell_spec, cone_locations, all_sig_stixels, varargin)
%
% inputs
%   datarun                 standard datarun structure
%   cell_spec               specify RGCs
%   cone_locations          Nx2 matrix of center points for N cones
%   all_sig_stixels         PxQ image of summed significant stixels
%
% optional inputs
%   rad             1       radius around which peak to calculate COM
%   center_method   com     function only supports COM now
%   verbose         false   print and plot extra stuff
%
% outputs
%   new_centers             Nx2 matrix of optimized center points for N
%                           cones
%
%  GDF: 2013-02
%

% parse inputs
p = inputParser;

p.addParamValue('verbose', false, @islogical)
p.addParamValue('rad', 1, @isnumeric)  % ROI radius around the cone
p.addParamValue('center_method', 'com') % method for optimizing the center point

p.parse(varargin{:});
verbose = p.Results.verbose;
rad = p.Results.rad;
center_method = p.Results.center_method;

% print shit to the screen
fprintf('\n optimizing cone locations \n')



% FUNCTION BEGINS HERE

% note cone number and rgc indices
num_cones = size(cone_locations,1);
rgc_indices = get_cell_indices(datarun, cell_spec);

% check to make sure the number of RGCs specified here is the same as that
% in all_sig_stixels
if length(rgc_indices) ~= size(all_sig_stixels,2)
    error('cell_spec does not specify the same number of RGCs as contained in all_sig_stixels')
end

% initialize variables (see their purposes below)
num_coi_per_cone = zeros(num_cones,1);
new_centers = zeros(size(cone_locations));

% setup status bar
dotter = floor(num_cones/20); % will put down 20 dots to indicate progress

% go through list of cones
for nn = 1:num_cones
    
    % advance "." every so often
    if mod(nn, dotter) == 0
        fprintf('.')
    end
    
    % get cone location
    cone_peak = cone_locations(nn,:);
    
    % If some cone positions have already been optimized, skip them
    if mod(cone_peak(1),1) ~= 0 || mod(cone_peak(2),1) ~=0
        new_centers(nn,:) = cone_peak;
        continue
    end

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
    num_coi_per_cone(nn) = length(cells_of_interest);
    
    % get rf patch around the cone from these rgcs
    cone_images = zeros(length(xrng), length(yrng), num_coi_per_cone(nn));
    for rgc = 1:num_coi_per_cone(nn)
        
        % get the correct index into datarun for the RGC
        rgc_index = rgc_indices(cells_of_interest(rgc));
        
        % get the spatial RF from the STA
        rf = get_rf(datarun, datarun.cell_ids(rgc_index));
        
        % make RGB into BW RF
        if size(rf, 3) == 3
            rf = sum(rf,3);
        end
        
        % normalized by robust_std to get into SNR units
        rf = rf ./ robust_std(reshape(rf, [],1))^2;
        
        % extract pixels over and around cone from RG
        cone_images(:,:,rgc) = rf(xrng,yrng);
    end
        
    % sum the cone_images together
    summed_cone_image = sum(cone_images,3);

    %% use the the info from all RGCs sampling this cone to determine its location
    
    % translate cone center point to coordinates of cone_map
    ctr = [cone_peak(1)-min(xrng)+1  cone_peak(2)-min(yrng)+1];
    
    % chose the method for determining cone center
    switch center_method
        case 'com'  % center of mass method
            
            % get COM
            [new_x,new_y] = ait_centroid(summed_cone_image);

            % plot change in center location
            if verbose
                figure(1)
                imagesc(summed_cone_image); colormap gray; hold on;
                plot(ctr(2),ctr(1), 'ro', new_x,new_y, 'go'); hold off;
            end
    end
            
    % translate back to display space
    new_centers(nn,:) = [new_y + min(xrng)-1, new_x + min(yrng)-1]; 

end    

% tell user its done
fprintf('\n finished optimizing locations \n');
    
    
    

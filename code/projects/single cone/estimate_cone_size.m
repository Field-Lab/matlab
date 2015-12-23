function cone_std = estimate_cone_size(datarun, cell_spec, cone_max_stixels, cone_locations, all_sig_stixels, varargin)
%
% Gets images for all the cones from RFs of cells that sample them strongly,
% up-samples these images, and superimposes them upon one another.  
% The "average" cone is then fit with a 1-D gaussian to estimate its width.  
%
% usage: cone_std = estimate_cone_size(datarun, cell_spec, cone_max_stixels, cone_locations, all_sig_stixels, varargin)
%
% inputs:
%   datarun             standard datarun structure
%   cell_spec           specify RGC
%   cone_max_stixels    Mx2 matrix of M cone peak stixels
%   cone_locations      Mx2 matrix of M cone coms
%   all_sig_stixels     PxQ image of summed significant pixels across RGCs
%
% optional inputs
%   rad             3       radius around cone center to cut cone image
%   upsampling      10      factor by which to upsample cone images
%   num_cones       []      if empty, use all, otherwise you can specify
%                           number
%   verbose         false   plot and print extra stuff
%
% outputs
%   cone_std            STD of the width of a gaussian fit to the average
%                           cone image
%
% GDF: 2013-02
%

% parse inputs
p = inputParser;
p.addParamValue('rad', 3, @isnumeric);
p.addParamValue('upsampling', 10, @isnumeric);
p.addParamValue('verbose', false, @islogical);
p.addParamValue('num_cones', [])
p.parse(varargin{:});


rad = p.Results.rad;
upsampling = p.Results.upsampling;
verbose = p.Results.verbose;


% get number of RGCs handed
rgc_indices = get_cell_indices(datarun, cell_spec);

% note number of cones
if isempty(p.Results.num_cones)
    num_cones = size(cone_max_stixels,1);
else
    num_cones = p.Results.num_cones;
end

% initialize a canvas upon which to put the upsampled cone RFs
mean_cone = zeros((rad*4+1)*upsampling);

% print to command line
fprintf('\n averaging cones \n')
% setup status bar
dotter = floor(num_cones/20); % will put down 20 dots to indicate progress


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
    cone_images = zeros(length(xrng), length(yrng), num_coi_per_cone);
    for rgc = 1:num_coi_per_cone
        
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
    
    % combine infor across RGCs that sample from the cone
    cone_image = sum(cone_images,3);

    % upsample image of cone
    cone_image = matrix_scaled_up(cone_image, upsampling);

    % translate cone center point to coordinates of cone_map
    ctr = [cone_locations(cn,1)-min(xrng)+1  cone_locations(cn,2)-min(yrng)+1];
    ctr = (ctr * upsampling) - (upsampling/2);
    ctr = round(ctr);
    
    % put cone image on canvas
    xdim = length(xrng)*upsampling;
    ydim = length(yrng)*upsampling;
    cvs_mdpnt = ((rad*4+1)*upsampling)/2;
    xbgn = cvs_mdpnt - ctr(1)+1;
    xend = cvs_mdpnt + (xdim-ctr(1));
    ybgn = cvs_mdpnt - ctr(2)+1;
    yend = cvs_mdpnt + (ydim-ctr(2));
    

    mean_cone(xbgn:xend,ybgn:yend) = mean_cone(xbgn:xend,ybgn:yend) + cone_image;
    
    if verbose
        imagesc(mean_cone); colormap gray
        pause(0.1)
    end
    
end
    
fprintf('\n finished optimizing locations \n');

% cut out an image from the canvas
ave_cone_image = mean_cone(cvs_mdpnt - (rad+0.5)*upsampling+1:cvs_mdpnt+(rad+0.5)*upsampling,cvs_mdpnt)...
                + mean_cone(cvs_mdpnt, cvs_mdpnt - (rad+0.5)*upsampling+1:cvs_mdpnt+(rad+0.5)*upsampling)';

norm_cone_image = ave_cone_image./max(ave_cone_image);
coef = [1,35,10,0.2];
myfun = @(a,x)a(1) * exp(-(x-a(2)).^2./(2*a(3).^2)) + a(4);


fitcoef = nlinfit(1:length(ave_cone_image),norm_cone_image', myfun, coef);

if verbose
    plot(norm_cone_image, 'b')
    fit_vals = myfun(fitcoef,1:length(ave_cone_image));
    hold on
    plot(fit_vals, 'k')
end

cone_std = fitcoef(3) ./ upsampling;

fprintf('finished \n')


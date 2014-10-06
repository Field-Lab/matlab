function [triplet,cone_kernel] = rgb_from_cones(cone_map, params)
% RGB_FROM_CONES     Extract RGB values of a cone using multiple stixels
%  
%   basic strategy is to combine information from multiple sources in two steps
%
% step 1 options (get multiple RGB triplets)
%   * regress each cell using specified shape (e.g. gaussian)
%   * collect all stixels in the region
%
% step 2 options (combine to get single triplet)
%   * pca
%   * sum
%
%
%
% usage:  triplet = rgb_from_cones(cone_map, ctr, params)
%
% arguments:  cone_map - 4d matrix of cone rf from multiple cells (y,x,color,cell)
%               params - struct of optional parameters (see below)
%
% outputs:     triplet - RGB summary of cone weights
%          cone_kernel - cone kernel that was used for regression
%
%
% optional fields in params, their default values, and what they specify:
%
% figure            []          which figure to plot in.  if 0, make new.  if empty, don't plot
% cone_expected     []          struct of expected RGB weights for L, M, and S cones.
%                                   cone_expected.L, etc
%                                   only used to plot in figure
%                                   if empty, not plotted
%
% regress           struct('type','dog')
%                               this provides the option to extract an RGB value from each cell before combining
%                               if empty, no regression is performed
%
%                               struct of parameters to generate the cone kernel which are regressed against the cells
%
%                                   type        how to generate cone kernel
%                                                   'dog'     - difference of gaussian
%                                                   'given'   - use the specified kernel
%
%
%                                   only required field is 'type'
%                                   all other fields are optional, with different default values depending on the type of kernel
%
%                                   for type = 'dog', see 'make_gaussian' for optional fields
%                                           x_size          3       size of cone kernel in x
%                                           y_size          3       in y
%                                           center_radius	1       center gaussian radius
%                                           surround_radius	2.5     surround gaussian radius
%                                           surround_scale  0       surround gaussian scale relative to center
%                                           center          [2 2]   location of center
%
%                                   for type = 'given', these are the fields and their default values
%                                           profile         [1]     the kernel to use.
%
% combine       'pca'           method to combine info across RGB values
%                                   'pca' - do pca on the RGB values
%                                   'sum' - add up all the values
%
%
%
% examples:
%
%   [triplet,ctr] = rgb_from_cones(cone_map, ctr, struct('regress',[]));
%
%
%
% gauthier 2008-09
%




% SET UP OPTIONAL ARGUMENTS

% if not specified, make params empty
if ~exist('params','var');params = [];end

% specify default parameters
defaults.figure = [];
defaults.cone_expected = [];
defaults.regress = struct('type','dog');
defaults.combine = 'pca';

% combine user and default parameters
params = default_params( defaults, params);


% get default kernel parameters
clear defaults
if ~isempty(params.regress)
    switch params.regress.type

        case 'dog'
            defaults.type = 'dog';
            defaults.x_size = 3;
            defaults.y_size = 3;
            defaults.center_radius = 1;
            defaults.surround_radius = 2.5;
            defaults.surround_scale = 0;
            defaults.center = [1 1]*2;

        case 'given'
            defaults.type = 'given';
            defaults.profile = 1;
            
        otherwise
            error('Kernel type ''%s'' not valid',params.filter.type)
    end

    % combine user and default parameters
    params.regress = default_params( defaults, params.regress);
end




% set up figure
if ~isempty(params.figure)
    if params.figure > 0
        figure(params.figure)
    else
        params.figure = figure;
    end
    figure(params.figure);clf
end



% GET RGB VALUES FROM ALL CELLS

if ~isempty(params.regress)

    % if a regression kernel is specified, get it and do regression
    switch params.regress.type
        
        case 'dog' % generate kernel from given parameters
            cone_kernel = make_gaussian(rmfield(params.regress,'type'));
            
        case 'given' % use given kernel for regression
            cone_kernel = params.regress.profile;

    end
    

    % make variable to store RGB values
    rgb_vals = zeros(size(cone_map,4),3);

    % regress the cone kernel against the data
    for cc = 1:size(cone_map,4)
        rgb_vals(cc,:) = reshape(cone_map(:,:,:,cc),[],3)'/reshape(cone_kernel,1,[]);
    end
    
    % take absolute value (to account for surround cones)
    rgb_vals = abs(rgb_vals);
    
else
    
    % if no regression kernel is specified, use all stixels
    rgb_vals = reshape(permute(cone_map,[1 2 4 3]),[],3);
    
    % give something to return for cone kernel
    cone_kernel = [];
    
end


% SUMMARIZE INFO ACROSS ALL RGB VALUES

switch params.combine
    
    case 'pca' % do PCA to extract vector direction
        
        % only do pca if there is more than one set of RGB values
        if size(rgb_vals,1) > 1

            temp = princomp(rgb_vals);
            triplet = temp(:,1)';

            % make positive
            if sum(triplet)<0;triplet = -1*triplet;end
        
        else
            % if only one set, use it
            triplet = rgb_vals;
        end
            

    case 'sum' % add up all RGB values
        triplet = sum(rgb_vals,1);
        
end


% plot result
if params.figure
    
    % go to figure for plotting cell
    figure(params.figure)

    % plot RGB values
    plot(rgb_vals(:,1),rgb_vals(:,2),'.k',[0 triplet(1)],[0 triplet(2)],'k')
    xlabel('red');
    ylabel('green');
    
    if ~isempty(params.cone_expected)
        
        % normalize weights
        Ln = params.cone_expected.L/sqrt(sum(params.cone_expected.L.^2));
        Mn = params.cone_expected.M/sqrt(sum(params.cone_expected.M.^2));
        Sn = params.cone_expected.S/sqrt(sum(params.cone_expected.S.^2));

        % add L-M cone lines to plot
        hold on;plot([0 10*Ln(1)],[0 10*Ln(2)],'r',[0 10*Mn(1)],[0 10*Mn(2)],'g');
    end

    drawnow
end



% plot cone profile averaged across cells
if 0 && length(cells_of_interest) > 1

    figure(105);clf;
    subplot(1,3,1);imagesc(cones_labeled(yrng,xrng)==nn);axis image;
    subplot(1,3,2);imagesc(frame_image(mean(cone_map,4)));axis image
    subplot(1,3,3);imagesc(ideal_cone);axis image
    pause

end


% plot cone profile from each cell
if 0 && length(cells_of_interest) > 0
    for cc=1:length(cells_of_interest)
        %figure(101);clf;imagesc(reshape(cone_vector,320,640))
        %figure(102);clf;imagesc(reshape(all_sig_stixels(:,cells_of_interest(cc)),320,640))

        figure(102);clf;
        subplot(1,3,1);imagesc(cones_labeled(yrng,xrng)==nn);axis image;
        subplot(1,3,2);imagesc(frame_image(cone_map(:,:,:,cc)));axis image
        subplot(1,3,3);imagesc(ideal_cone);axis image
        %pause
    end
end


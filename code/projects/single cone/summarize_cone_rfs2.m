function [cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_rfs, cone_ids] =...
    summarize_cone_rfs2(datarun, cluster_centers, cluster_sources, kernel_params, varargin)
% summarize_cone_rfs     Identify precise location and RGB profile of each cone
%
% usage:  [cone_rgb, cone_spatial_profiles, cone_ideal_shapes, cone_ids] =...
%    summarize_cone_rfs(datarun, cell_spec, all_cone_centers, cones_labeled, all_sig_stixels, params)
%
%
% arguments:  datarun - datarun struct with field specifying X
%     cluster_centers - Nx2 matrix: x,y locations of centers of clusters of local maxima
%     cluster_sources - N length vector: elements in 1:C describing which cell ids gave rise to which cone centers
%       cone_clusters - N vector, which cluster each cone belongs to.
%       kernel_params - struct describing cone kernel.  the cone kernel is used to generate the
%                           cone spatial profile, and possibly to regress against cones.
%                           the required fields depend on the value of field 'type'
%                           for 'type' = 'dog', these are required fields
%                               kernel_params.center_radius
%                               kernel_params.surround_radius
%                               kernel_params.surround_scale
%
%              params - struct of optional parameters (see below)
%
%
% outputs:
%
%               cone_rgb - Nx3 matrix: RGB triplets for each cone
%  cone_spatial_profiles - MxN matrix: pixel maps of the ideal shape of each cone
%      cone_ideal_shapes - N length cell array: analytic descriptions of the ideal shape of each cone
%                               e.g. struct with parameters for a 2D gaussian
%               cone_ids - N length vector: which cone ids are described in the preceding variables.  if no roi, will be [1:num_cones]
%
%
% optional fields in params, their default values, and what they specify:
%
% centers               'fit'       how to refine cone center points
%                                       'fit' - fit a gaussian
%                                       'com' - center of mass
%                                       'fixed' - don't refine
% sensitivity           struct      params struct to pass to find_cones_in_rf, used for identifying cone center
% roi                   []          region in which to use cones
% regress               true        whether to regress the cone spatial profile against measured cone RFs to extract RGB strengths
% combine               'sum'       how to combine multiple RGB triplets into a single triplet for each cone
%
% single_cone_figure    []          figure number to plot fits of each cone serially.  if 0, make new.  if empty, don't plot
% cone_weights_figure   []          figure number to plot RGB from all cones.  if 0, make new.  if empty, don't plot
% cone_remap            []          structure with fields:
%                                       fcn       - function mapping from 3D to 2D to plot RGB from all cones
%                                       x_caption - caption for x axis
%                                       y_caption - caption for y axis
% verbose               false       display in text what's going on
%
%
% note: if it exists, datarun.cones.rgb_expected will be used to plot expected cone weights compared to observed weights
%
% gauthier 2008-09
%


% * add color tag to each cone




% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('verbose', false);
p.addParamValue('centers', 'fit');
p.addParamValue('sensitivity', struct);
p.addParamValue('roi', []);
p.addParamValue('regress', true);
p.addParamValue('combine', 'sum');
p.addParamValue('single_cone_figure', []);
p.addParamValue('foa_fit', []);
p.addParamValue('cone_weights_figure', []);
p.addParamValue('cone_remap', []);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;




% BODY OF THE FUNCTION

% set up figure
pa_fit = set_up_fig_or_axes(params.foa_fit);

% get expected cone RGB values
rgb_expected = cone_rgb_expected(datarun);




% all_cone_centers lists cone center points for all cones which were found,
% but only those within the ROI will be anlayzed

% determine which clusters have at least one cone within the ROI.
% this information is stored as a list of indices to the rows all_cone_centers
if isempty(params.roi)
    % if roi is empty, use all cones
    cone_ids = 1:size(cluster_centers,1);
else
    % if there is a ROI, get the indices of cones within it
    cone_ids = find(diag(params.roi(round(cluster_centers(:,2)),round(cluster_centers(:,1)))))';

    % pare down to just these clusters
    %cluster_centers = cluster_centers(cone_ids,1:2);
    %cluster_sources = {cluster_sources{cone_ids}};
end


% note how many cones there will be (one cone per cluster)
num_cones = length(cone_ids);


if params.verbose
    fprintf('\nLooking through %d clusters...',num_cones)
    start_loading = clock; % note when it started
end


% initialize variables (see their purposes below)
initial_rgb = zeros(num_cones,3);
cone_spatial_profiles = sparse(datarun.stimulus.field_height*datarun.stimulus.field_width,num_cones);
cone_rgb = zeros(num_cones,3);
cone_rfs = cell(num_cones,1);
cone_ideal_shapes = cell(num_cones,1);


% get RGB estimate for each cone, to seed subsequent fit

% go through list of cones
for nn = []%1:num_cones

    % get the id number of this cone
    cone_id = cone_ids(nn);

    % get center
    center = cluster_centers(cone_id,:);

    % get RFs
    [rfs,x_offset,y_offset] = rfs_in_roi(datarun,cluster_sources{cone_id},center,1);

    % ensure at least some were found
    if isempty(rfs)
        error('cone id %d was not sampled significantly by any ganglion cell with an RF.',cone_id)
    end

    % get cone fit
    switch kernel_params.type
        case 'dog'
            cone_fit = make_gaussian('center',center - [x_offset y_offset],...
                'y_size',size(rfs,1),'x_size',size(rfs,2),'normalize','sum',...
                rmfield(kernel_params,'type'));
        otherwise
            error('cone kernel type %s not recognized',cone_kernel_params.type)
    end

    % identify cone's RGB
    initial_rgb(nn,:) = sum(reshape(reshape(rfs,[],size(rfs,3)*size(rfs,4))'*reshape(cone_fit,[],1),size(rfs,3),size(rfs,4)),2);

    % identify which cone types those are
    %cone_types = classify_cones(initial_rgb(nn,:), kernel_params.cone_rgb, 'algorithm','nearest line');

end





% go through list of cones
for nn = 1:num_cones

    % make progress tick
    %fprintf('.')
    %fprintf('%d',nn)
    
    % get the id number of this cone
    cone_id = cone_ids(nn);

    % get center point
    center = cluster_centers(cone_id,:);
    
    


    % GATHER CONE RFS (i.e. the stixels around this cone in each RGC RF)

    % ROI radius
    rad = 1;

    % get RFs
    [rfs,x_offset,y_offset,xrng,yrng] = rfs_in_roi(datarun,cluster_sources{cone_id},center,rad);

    % get cone fit parameters
    switch kernel_params.type
        case 'dog'

    end

    % put sum of individual RFs of this cone into a sparse matrix the size of an entire RGC RF

    % initialize variable
    cone_map_sum_r = sparse(datarun.stimulus.field_height,datarun.stimulus.field_width);
    cone_map_sum_g = cone_map_sum_r;
    cone_map_sum_b = cone_map_sum_r;

    % enter the sum of cone_map in the appropriate location
    cone_map_sum_r(yrng,xrng) = sum(rfs(:,:,1,:),4);
    cone_map_sum_g(yrng,xrng) = sum(rfs(:,:,2,:),4);
    cone_map_sum_b(yrng,xrng) = sum(rfs(:,:,3,:),4);

    % save to cone_rfs
    cone_rfs{nn}.cone_map_sum_r = cone_map_sum_r;
    cone_rfs{nn}.cone_map_sum_g = cone_map_sum_g;
    cone_rfs{nn}.cone_map_sum_b = cone_map_sum_b;





    % choose cone profile shape


    % note the ideal cone shape parameters
    switch kernel_params.type
        case 'dog'
            % save DOG parameters
            cone_ideal = rmfield(kernel_params,'type');

        otherwise
            error('cone kernel type ''%s'' not recognized',kernel_params.type)
    end

    
    
    switch params.centers
        case 'fit neighborhood'  % fit all cones in the neighborhood simultaneously

            % parameters
            fit_radius = 3;



            % COLLECT RELEVANT RFs AND CLUSTER CENTERS
            
            % identify clusters in the ROI
            in_clusters =  find(round(cluster_centers(:,1))>=center(1)-fit_radius &...
                round(cluster_centers(:,1))<=center(1)+fit_radius &...
                round(cluster_centers(:,2))>=center(2)-fit_radius &...
                round(cluster_centers(:,2))<=center(2)+fit_radius);
            
            % stick to the overall cone ROI
            in_clusters = intersect(in_clusters,cone_ids);
            
            % note which is this one
            this_clust = find(in_clusters == cone_id);
            
            % get these RFs which gave rise to these clusters
            [rfs,x_offset,y_offset,xrng,yrng] = rfs_in_roi(datarun,unique([cluster_sources{in_clusters}]),center,fit_radius);
            
            
            % convert center points to  revelant coordinates
            these_centers = cluster_centers(in_clusters,:);
            these_centers(:,1) = these_centers(:,1) - x_offset;
            these_centers(:,2) = these_centers(:,2) - y_offset;
            
            % plot them
            if 0
                for rr=1:size(rfs,4)
                    figure(10);clf;
                    % plot rf
                    imagesc(norm_image(rfs(:,:,:,rr)));axis image; hold on;
                    % plot all clusters
                    plot(these_centers(:,1),these_centers(:,2),'.r')
                    % plot the cluster of interest
                    plot(these_centers(this_clust,1),these_centers(this_clust,2),'.m')
                    axis image
                    pause
                end
            end
            
            
            % FIT THEM
            [fit_cone_info,rf_fits] = fit_cones_in_many_RFs4(rfs,these_centers,'figure',[]);
            %pause
            %continue
            
            %[fit_cone_info,rf_fit] = fit_cones_in_many_RFs(rfs,cone_initial,'figure',[],...
            %    'cone_kernel',kernel_params,'optim',{'TolX',0.1,'Display','off','MaxIter',50});




            % save the parameters

            % use the fit for the cone kernel, normalizing the sum to one
            cone_profile_small  = make_gaussian('center',fit_cone_info(this_clust,1:2),...
                'y_size',length(yrng),'x_size',length(xrng),'normalize','sum',...
                rmfield(kernel_params,'type'));

            % save center point
            cone_ctr = fit_cone_info(this_clust,1:2) + [x_offset y_offset];

            
            % get RGB triplet

            cone_rgb(cone_id,:) = fit_cone_info(this_clust,3:5);




            % plot RGB values from all cones so far

            % be sure plotting is desired, and a function to remap to 2D exists
            if ~isempty(params.cone_weights_figure) && ~isempty(params.cone_remap)

                % map 3D cone weights to 2D
                remapped_rgb = params.cone_remap.fcn(cone_rgb(1:nn,1:3));

                % plot in 2D
                figure(params.cone_weights_figure);clf
                plot(remapped_rgb(:,1),remapped_rgb(:,2),'.k')
                xlabel(params.cone_remap.x_caption)
                ylabel(params.cone_remap.y_caption)

                % plot expected gun strengths (if present in datarun)
                if isfield(datarun,'cones') && isfield(datarun.cones,'rgb_expected') && ~isempty(datarun.cones.rgb_expected)
                    % get L and M center in remapped coordinates
                    L_ctr = params.cone_remap.fcn(rgb_expected.L);
                    M_ctr = params.cone_remap.fcn(rgb_expected.M);
                    % add them to the plot
                    hold on
                    plot(L_ctr(1),L_ctr(2),'.r','MarkerSize',25)
                    plot(M_ctr(1),M_ctr(2),'.g','MarkerSize',25)
                end

                % set axes
                set(gca,'XLim',[-.5 1.5],'YLim',[-.5 1.5])

                drawnow
            end






        case 'fit'  % fit gaussian to the cone_sensitivity combined across cells


            % get RF to fit
            rf_to_fit = sum(cone_sensitivity,3);

            if 0
                % collapse cone_map across cells and color

                % add up the version of the cone from each cell
                %rf_to_fit = sum(cone_map,4);

                % add up across colors
                %rf_to_fit = sum(rf_to_fit,3);

                % OLD: get vector length of each stixel
                %rf_to_fit = sqrt(sum(rf_to_fit.^2,3));
            end


            % set up parameters of the shape to fit to the cone

            switch kernel_params.type
                case 'dog'
                    % save DOG parameters
                    fit_params = kernel_params;
                    fit_params.effective_radius = Inf;
            end

            % get the sum of rf_to_fit to use as the initial guess about cone height
            rf_sum = sum(sum(rf_to_fit));

            % get the max of rf_to_fit to use as the initial guess about cone height
            rf_max = max(max(rf_to_fit));



            % fit the cone location

            % compute fit, using the cone center as initial center and rf_sum as the initial height
            [fit_cone_info,rf_fit] = fit_cones(rf_to_fit,[ctr rf_sum],....
                struct('cone_kernel',fit_params,'figure',[],'fig_final',[]));



            % save the parameters

            % use the fit for the cone kernel, normalizing the sum to one
            cone_profile_small  = make_gaussian('center',fit_cone_info(1:2),...
                'y_size',size(rf_to_fit,1),'x_size',size(rf_to_fit,2),'normalize','sum',...
                rmfield(fit_params,'type'));

            % save center point
            cone_ctr = fit_cone_info(1:2);



            % plot the fit and center point before and after
            if 0
                figure(11);clf;imagesc(cone_profile_small);colormap gray
                hold on; plot(fit_cone_info(:,1),fit_cone_info(:,2),'.r',ctr(1),ctr(2),'.g');
                pause
            end


        case 'com'  % take the center of mass of cone_sensitivity


            % get RF to fit
            rf_for_com = sum(cone_sensitivity,3);

            % get COM
            [com_x,com_y] = ait_centroid(rf_for_com);

            switch kernel_params.type
                case 'dog'
                    params_temp = rmfield(kernel_params,'type');
                    params_temp.center = [com_x com_y];
                    params_temp.x_size = size(cone_map,2);
                    params_temp.y_size = size(cone_map,1);
            end

            % generate cone profile from DOG parameters
            cone_profile_small = make_gaussian(params_temp);

            % note center point
            cone_ctr = [com_x com_y];


        case 'fixed'  % don't fit, just use the location from cone_centers

            switch kernel_params.type
                case 'dog'
                    params_temp = rmfield(kernel_params,'type');
                    params_temp.center = ctr;
                    params_temp.x_size = size(cone_map,2);
                    params_temp.y_size = size(cone_map,1);
            end

            % generate cone profile from DOG parameters
            cone_profile_small = make_gaussian(params_temp);

            % note center point
            cone_ctr = ctr;








        otherwise
            error('center point finding algorithm ''%s'' not recognized.',params.centers)

    end



    % save the cone spatial profile, both in pixels and analytically
    if 1

        % put the ideal cone shape (cone_profile_small) in a sparse matrix as big as a whole sta
        ideal_cone_large = sparse(datarun.stimulus.field_height,datarun.stimulus.field_width);
        ideal_cone_large(yrng,xrng) = cone_profile_small;

        % save this ideal cone profile in the final output matrix
        cone_spatial_profiles(:,nn) = reshape(ideal_cone_large,[],1);

        % note the cone profile center point (translated to stimulus coordinates)
        cone_ideal.center = cone_ctr;

        % put the analytic description in the final output cell array
        cone_ideal_shapes{nn} = cone_ideal;


        % compare before and after
        if 0
            fprintf('cone id %d: from [%0.2f %0.2f] to [%0.2f %0.2f]\n',cone_id,...
                cluster_centers(cone_id,1),cluster_centers(cone_id,2),cone_ideal.center(1),cone_ideal.center(2))
        end

    end

    


    if mod(nn,50)==0
        fprintf('finished %d in %0.1f sec (projected %0.1f sec remaining)\n',nn,etime(clock,start_loading),...
            (num_cones-nn)*etime(clock,start_loading)/nn)
    end


end



if params.verbose
    fprintf(' done (%0.1f seconds)\n',etime(clock,start_loading));
end






function [rfs,x_offset,y_offset,xrng,yrng] = rfs_in_roi(datarun,cell_ids,center,rad)
% get RFs in a ROI for each cell, and the specified center point in local coordinates
% also returns the x and y offset of coordinates, so that e.g. center points
% can be remapped into the x-y coordinates of 'rfs'

% store in variable called rfs, which has four dimensions
% first three dimensions = each cell's rf in the ROI (x,y,color)
% fourth dimension indexes the RGCs the rf come from

% get region
xrng = max(round(center(1)-rad),1) : min(round(center(1)+rad),datarun.stimulus.field_width);
yrng = max(round(center(2)-rad),1) : min(round(center(2)+rad),datarun.stimulus.field_height);

% note offsets
x_offset = xrng(1) - 1;
y_offset = yrng(1) - 1;

% initialize variables
rfs = zeros(length(yrng),length(xrng),3,length(cell_ids));
has_rf = true(length(cell_ids),1);

% get stixels of the ROI from the relevant RGCs
for cc=1:length(cell_ids)

    % get the rf frame
    rf = get_rf(datarun,cell_ids(cc));

    % skip if empty
    if isempty(rf);
        has_rf(cc) = false;
        continue
    end

    % get pixels from the ROI, normalized by the noise level
    rfs(:,:,:,cc) = rf(yrng,xrng,:)/ robust_std(reshape(rf,1,[]));
end

% toss cells which didn't have an RF
rfs = rfs(:,:,:,has_rf);

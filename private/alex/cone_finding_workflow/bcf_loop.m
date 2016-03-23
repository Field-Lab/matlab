function bcf = bayesian_cone_finding_loop(datarun,bcf_params)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   BAYESIAN CONE FINDING STEP 3 OF 5    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% bayesian_cone_finding     perform bayesian cone finding on many small patches
%
%
% usage:  bcf = bayesian_cone_finding(datarun,bcf_params)
%
% arguments:     datarun - datarun struct
%             bcf_params - struct of parameters
%
% outputs:     bcf - struct of results
%
%
%
% 2010-02  gauthier
%





% EXPAND PARAMETERS

fn = fieldnames(bcf_params);
for ff=1:length(fn)
    eval(sprintf('%s = bcf_params.%s;',fn{ff},fn{ff}))
end

distance_prior = C_C;

% choose shape of sigmoid function
% sigmoid(0) = 0, sigmoid(1) = 1
switch 2
    case 1 % exponential
        sigmoid = @ (x)(1 ./ (1 + exp(-10*x+5)));
    case 2 % parabolas
        sigmoid = @ (x) 2*(  (x<0.5).*(x.^2) + (x>=0.5).*(-(x-1).^2+0.5)  );
end
% if x is between 0 and 1 then use sigmoid, else x is 1 or realmin (i.e. very close to zero, while still having its log defined)
fn  = @ (x,c1,c2)((x<=c1).*(realmin) + (x>=c2).*(1) + ((x>c1).*(x<c2)).*(sigmoid((x-c1)./(c2-c1))));


% make filter based on distance prior
% this filter (a 2D matrix) will be convolved with the cone locations
% to determine the repulsion radius between cones of each types

% stores the sigmoid function specific to the each cone type interaction
dist_filt = cell(1);
% get distances, in units of kernels
c1 = distance_prior(1,1,1)/kernel_spacing;
c2 = distance_prior(1,2,1)/kernel_spacing;

% set up parameters for the filter
largest_rad = max(reshape(distance_prior/kernel_spacing,[],1));
filt_diam = (1+2*ceil(largest_rad));
filt_center = (1+ceil(largest_rad));

% make filter for these cone types, load into dist_filt
dist_filt{1} = fn(distance_from_point([filt_diam filt_diam],[filt_center filt_center]),c1,c2);



cell_ids = get_cell_indices(datarun, cone_finding_cell_spec);
num_frames = datarun.duration * datarun.stimulus.monitor_refresh / datarun.stimulus.interval;
stim_variance = (.5*.5*.96)^2;

rel_filt = fspecial('gaussian',rel_radius*4,rel_radius);
rel_filt = double(rel_filt > 0.5 * max(max(rel_filt)));


% load each STA in the ROI
tic
fprintf('\nProcessing STAs... \n\n')
wid = size(datarun.stas.stas{cell_ids(1)},1);
hei =size(datarun.stas.stas{cell_ids(1)},2);
processed_stas = zeros(wid,hei,length(cell_ids));
stas_matrix = zeros(wid,hei,length(cell_ids)); % the one which wwill be used
cell_constants = zeros(length(cell_ids),1);
for cc = 1:length(cell_ids)
    cell_index = cell_ids(cc);
    if ~isempty(datarun.stas.time_courses{cell_index})
        
        % get constant
        fit_b = datarun.stas.snls{cell_index}.fit_params.b;
        cell_constant = num_frames * stim_variance * exp(fit_b);
        n_spikes = length(datarun.spikes{cell_index});
        cell_constants(cc) = cell_constant;
        % get sta
        sta = datarun.stas.stas{cell_index};        
        % force RGB sta into BW format
        if size(sta,3) == 3 % if STA is RGB
            sta = mean(sta,3);
        end
        % reshape, so that first dim is space-color, second dim is time
        sta_r = reshape(sta,[],size(sta,4));
        % get timecourse, and set its norm to 1
        tc = sum(datarun.stas.time_courses{cell_index},2);
        tc = tc/sqrt(sum(tc.^2));
        % extract spatial component, and reshape to standard 3d matrix (y,x,color)
        sta_spatial = sta_r*tc;
        sta_spatial = reshape(sta_spatial,datarun.stimulus.field_height,datarun.stimulus.field_width,[]);
        % use this line to reconstruct the STA (as a product of the spatial structure and timecourse
        %sta_recon = kron(sta_spatial,tc');
        % find significant stixels
        %sig_stix = full(significant_stixels(get_sta(datarun,cell_id)));
        sig_stix = full(datarun.stas.marks{cell_index});
        % zero out where stixels are not significant
        sta_spatial = sta_spatial .* imfilter(sig_stix,rel_filt);
        
        processed_stas(:,:,cc) = sta_spatial;
        stas_matrix(:,:,cc) = (n_spikes / cell_constant) *sta_spatial;
    end
end
toc

% CUT UP INTO SMALL REGIONS

% count ROIs needed
num_roi_x = ceil( diff(relevant_region_x) / (roi_x_size - 2*padding_x) );
num_roi_y = ceil( diff(relevant_region_y) / (roi_y_size - 2*padding_y) );

roi_x = [relevant_region_x(1) relevant_region_x(1) + roi_x_size - 1];
roi_y = [relevant_region_y(1) relevant_region_y(1) + roi_y_size - 1];
% get cone centers: regular lattice, in stixel coordinates
num_kern_x = length(roi_x(1):kernel_spacing:roi_x(2));
num_kern_y = length(roi_y(1):kernel_spacing:roi_y(2));
[kernel_x, kernel_y] = meshgrid(roi_x(1):kernel_spacing:roi_x(2),roi_y(1):kernel_spacing:roi_y(2));
kernel_x = reshape(kernel_x,[],1);
kernel_y = reshape(kernel_y,[],1);
% size of the ROI
rf_size = [diff(roi_y)+1 diff(roi_x)+1];
% generate spec for each kernel
% locations
kernel_centers = 1 + [kernel_x-roi_x(1) kernel_y-roi_y(1)];
kernel_spec = struct('center',mat2cell(kernel_centers,ones(size(kernel_centers,1),1)));
[kernel_spec.radius] = deal(kernel_radii(1));

tic
[W,kernel_norms] = make_cone_weights_matrix_bw(rf_size,kernel_spec);
toc

% initialize variable
rois_x = zeros(num_roi_x*num_roi_y,2);
rois_y = rois_x;
num_rfs_roi=cell(num_roi_x*num_roi_y,1);
loop_cones = cell(num_roi_x*num_roi_y,1);
dlls = cell(num_roi_x*num_roi_y,1);


% if ~isempty(cones_fig)
%     figure(cones_fig);clf;
%     image(0.5*ones(datarun.stimulus.field_height,datarun.stimulus.field_width,3)); axis image; hold on;
%     cone_axes = gca;
%     drawnow
% end

% loop through each non-zero
rr = 1;
for xx = 1:num_roi_x
    for yy = 1:num_roi_y
        start_x = relevant_region_x(1) + (xx-1)*(roi_x_size - 2*padding_x);
        start_y = relevant_region_y(1) + (yy-1)*(roi_y_size - 2*padding_y);
        rois_x(rr,:) = [start_x start_x + roi_x_size - 1];
        rois_y(rr,:) = [start_y start_y + roi_y_size - 1];
        
        stas = stas_matrix(rois_y(rr,1):rois_y(rr,2),rois_x(rr,1):rois_x(rr,2),:);
        [~,~,d] = ind2sub(size(stas), find(stas));
        stas = stas(:,:,unique(d));
        stas = reshape(stas,[],length(unique(d)));
        
        if ~isempty(stas)
            
            [added_cones,first_dll] = bayesian_cone_finding_bw(bcf_params,...
                rois_x(rr,:),rois_y(rr,:),W,kernel_norms, cell_constants(unique(d)), stas, dist_filt);
            
            % save results
            loop_cones{rr} = added_cones;
            dlls{rr} = first_dll;
            
            % add cones to plot
%             if ~isempty(added_cones)
%                 cns = find(added_cones(:,3));
%                 if ~isempty(cones_fig)
%                     plot(added_cones(cns,4),added_cones(cns,5),'.','Color',kernel_plot_colors,'Parent',cone_axes)
%                 end
%                 drawnow
%             end
        end
        num_rfs_roi{rr} = unique(d);
        rr = rr + 1;
    end
end

my_rois=find(~cellfun(@isempty, num_rfs_roi));


% CONCATENATE DELTA-LOG-LIKELIHOOD IMAGES and ACCUMULATE CONES
if 1
    
    % initialize
    dll = zeros((datarun.stimulus.field_height-1)/kernel_spacing+1,(datarun.stimulus.field_width-1)/kernel_spacing+1,size(kernel_plot_colors,1));
    all_added_cones = [];
    
    % go through each subset
    for rr = my_rois' %1:size(rois_x,1)
        
        % likelihood
        
        % identify where this region fits
        x_range = [rois_x(rr,1) + padding_x - 1   rois_x(rr,2) - padding_x - kernel_spacing]/kernel_spacing;
        y_range = [rois_y(rr,1) + padding_y - 1   rois_y(rr,2) - padding_y - kernel_spacing]/kernel_spacing;
        x_range = round(x_range);
        y_range = round(y_range);
        
        % reshape dll
        temp = reshape(dlls{rr},num_kern_y,num_kern_x);
        
        % put it in
        dll(y_range(1):y_range(2),x_range(1):x_range(2)) = ...
            temp(1+(padding_y-1)/kernel_spacing:end-padding_y/kernel_spacing-1,...
            1+(padding_x-1)/kernel_spacing:end-padding_x/kernel_spacing-1);
        
        
        % cones
        
        % set up ROI
        roi_x = rois_x(rr,:) + [padding_x-1 -padding_x];
        roi_y = rois_y(rr,:) + [padding_y-1 -padding_y];
        
        % identify cones in the ROI
        if ~isempty(loop_cones{rr})
            added_cones = loop_cones{rr}( ...
                (loop_cones{rr}(:,4)>=roi_x(1)) & (loop_cones{rr}(:,4)<roi_x(2)) & ...
                (loop_cones{rr}(:,5)>=roi_y(1))  & (loop_cones{rr}(:,5)<roi_y(2)) ...
                ,:);
            
            % keep this in growing list
            all_added_cones = [all_added_cones ; added_cones]; %#ok<AGROW>
        end
        
        % add boundaries of the roi
        %plot([roi_x roi_x(2) roi_x(1) roi_x(1)],[roi_y(1) roi_y(1) roi_y(2) roi_y(2) roi_y(1)],'w')
        
    end
    
end



% package up the results
bcf.all_added_cones = all_added_cones;
bcf.dll = dll;
bcf.rois_x = rois_x;
bcf.rois_y = rois_y;


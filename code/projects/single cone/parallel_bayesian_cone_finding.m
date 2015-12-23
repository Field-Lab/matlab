function found_cones = parallel_bayesian_cone_finding(varargin)
%
%
% bayesian_cone_finding     perform bayesian cone finding on a small patch
%
% usage:  bcf = bayesian_cone_finding(datarun,bcf_params)
%
% arguments:     datarun - datarun struct
%             bcf_params - struct of parameters
%
% outputs:     result - struct of results
%
%
%
% 2010-02  gauthier
%

p = inputParser;
p.addParamValue('bcf_params', [], @isstruct);
p.addParamValue('roi_x', [], @isnumeric);
p.addParamValue('roi_y', [], @isnumeric);
p.addParamValue('W', [], @isnumeric);
p.addParamValue('kernel_norms', [], @isnumeric);
p.addParamValue('datarun_path', [], @ischar);

p.parse(varargin{:});

bcf_params = p.Results.bcf_params;
roi_y = p.Results.roi_y;
roi_x = p.Results.roi_x;
W = p.Results.W;
kernel_norms = p.Results.kernel_norms;

% load datarun
load(p.Results.datarun_path);

make_plots = false;
loop_mode = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EXPAND PARAMETERS


% get passed parameters
fn = fieldnames(bcf_params);
for ff=1:length(fn)
    eval(sprintf('%s = bcf_params.%s;',fn{ff},fn{ff}))
end
    


% set up distance prior, a 3D matrix specifying d_min & d_max for each pair of cone types
% dim1,dim3: which cone types
% dim2: [c1 c2] = [d_min d_max]
% e.g. distance_prior(2,1,3) = distance_prior(3,1,2) = d_min for M-S cone interaction

clear distance_prior
if length(fieldnames(bcf_params.kernel_colors)) == 1  % BW stimulus
    distance_prior = C_C;
else
    distance_prior(:,:,1) = [LM_MM; LM_MM; LM_S];
    distance_prior(:,:,2) = [LM_MM; LM_MM; LM_S];
    distance_prior(:,:,3) = [LM_S; LM_S; S_S];
end

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

% number of cone kernel colors
kern_c = size(distance_prior,3);
% stores the sigmoid function specific to the each cone type interaction
dist_filt = cell(kern_c,kern_c);
% loop through cone types
for source_color = 1:kern_c
    for repulsed_color = 1:kern_c
        % get distances, in units of kernels
        c1 = distance_prior(source_color,1,repulsed_color)/kernel_spacing;
        c2 = distance_prior(source_color,2,repulsed_color)/kernel_spacing;

        % set up parameters for the filter
        largest_rad = max(reshape(distance_prior/kernel_spacing,[],1));
        filt_diam = (1+2*ceil(largest_rad));
        filt_center = (1+ceil(largest_rad));

        % make filter for these cone types, load into dist_filt
        dist_filt{source_color,repulsed_color} = ...
            fn(distance_from_point([filt_diam filt_diam],[filt_center filt_center]),c1,c2);
    end
end

% make ROI in which to use pixels for analysis
% initialize to all 0s
roi = false(datarun.stimulus.field_height,datarun.stimulus.field_width);
% set a block of values to 1
roi(roi_y(1):roi_y(2),roi_x(1):roi_x(2)) = true;
% count how many values are in there and
% multiply by 3 (STAs must be RGB)
num_stixels = sum(sum(roi)) * 3;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GENERATE MATRICES

% note locations of kernels in original cone space
if new_W || ~exist('W','var') || isempty(W)  || (exist('loop_mode','var') && loop_mode == 1)
    
    % note the available cone colors
    kernel_color_names = fieldnames(kernel_colors);
    %kern_c = length(kernel_color_names);

    % get cone centers
    switch 1
        case 1 % regular lattice, in stixel coordinates
            num_kern_x = length(roi_x(1):kernel_spacing:roi_x(2));
            num_kern_y = length(roi_y(1):kernel_spacing:roi_y(2));
            [kernel_x, kernel_y] = meshgrid(roi_x(1):kernel_spacing:roi_x(2),roi_y(1):kernel_spacing:roi_y(2));
            kernel_x = reshape(kernel_x,[],1);
            kernel_y = reshape(kernel_y,[],1);
    end
    num_kernels = length(kernel_x)*kern_c;
end



% generate W
num_cone_types = length(fieldnames(kernel_colors));

if new_W || ~exist('W','var') || isempty(W)

    
    % size of the ROI
    rf_size = [diff(roi_y)+1 diff(roi_x)+1];

    % generate spec for each kernel
    
    % locations
    kernel_centers = 1 + [kernel_x-roi_x(1) kernel_y-roi_y(1)];
    % replicate, once for each cone type
    kernel_centers = reshape(reshape(repmat(kernel_centers,1,kern_c)',[],1),2,[])';
    % load into struct
    %kernel_spec = struct('center',mat2cell(kernel_centers,repmat(1,size(kernel_centers,1),1)));
    kernel_spec = struct('center',mat2cell(kernel_centers,ones(size(kernel_centers,1),1)));
    
    % types
    for cc = 1:kern_c; [kernel_spec(cc:num_cone_types:end).type] = deal(kernel_color_names{cc}); end
    
    % radii
    for cc = 1:kern_c; [kernel_spec(cc:num_cone_types:end).radius] = deal(kernel_radii(cc)); end
    
    % compute W and norm of each kernel
    [W,kernel_norms] = make_cone_weights_matrix(rf_size,kernel_spec,kernel_colors);

end



% generate STAs matrix

% get cell IDs in ROI
cell_ids = datarun.cell_ids(get_cell_indices_roi(datarun,cone_finding_cell_spec,roi,'cell_locations','marks'));
num_rfs = length(cell_ids);

% get SNLs (compute from scratch if needed)
datarun = get_snls(datarun, cell_ids,'frames',-2:0,'start_time',start_time,'stimuli',10000,'new',false);




% exit loop if no RFs are in this ROI
if num_rfs == 0
    disp('no RFs found in region')
    
    added_cones = [];
    first_dll = zeros(1,num_kern_y*num_kern_x*kern_c);

    output_info.added_cones = added_cones;
    output_info.first_dll = first_dll;
    found_cones{1} = output_info;
    
    return
end


% fill in STAS matrix
if new_STAs || ~exist('stas_matrix','var')
    
    % put relevant parts of RFs into matrix form

    % compute constant
    num_frames = datarun.duration * datarun.stimulus.monitor_refresh / datarun.stimulus.interval;
    stim_variance = (.5*.5*.96)^2;

    % initialize matrices
    raw_stas_matrix = sparse(num_stixels,num_rfs);
    stas_matrix = sparse(num_stixels,num_rfs);
    cell_constants = zeros(num_rfs,1);
    fit_b = zeros(num_rfs,1);
    n_spikes = zeros(num_rfs,1);

    % make circle of radius rel_radius
    if rel_radius < Inf
        rel_filt = fspecial('gaussian',rel_radius*4,rel_radius);
        rel_filt = double(rel_filt > 0.5 * max(max(rel_filt)));
    else
        rel_filt = ones(roi_y(2)-roi_y(1)+1,roi_x(2)-roi_x(1)+1);
    end

    %T=text_waitbar('Loading STAs...');

    % load each STA in the ROI
    for cc = 1:num_rfs
        cell_id = cell_ids(cc);
        cell_index = get_cell_indices(datarun,cell_id);
        % get constant
        fit_b(cc) = datarun.stas.snls{cell_index}.fit_params.b;
        cell_constants(cc) = num_frames * stim_variance * exp(fit_b(cc));
        % get sta
        sta = get_sta(datarun,cell_id);
        % force a BW sta into RGB format 
        if size(sta,3) == 1 % if STA is BW
            sta = repmat(sta, [1,1,3,1]);
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
        sig_stix = full(datarun.stas.marks{get_cell_indices(datarun,cell_id)});
        % zero out where stixels are not significant
        sta_spatial = sta_spatial .* imfilter(repmat(sig_stix,[1 1 3]),rel_filt);
        % pare to ROI
        sta_spatial = sta_spatial(roi_y(1):roi_y(2),roi_x(1):roi_x(2),:);
        % get number of spikes
        n_spikes(cc) = length(datarun.spikes{get_cell_indices(datarun,cell_id)});
        % reshape in ROI, and put into STAs matrix
        raw_stas_matrix(:,cc) = reshape(sta_spatial,1,[]);
        stas_matrix(:,cc) = (n_spikes(cc) / cell_constants(cc)) * reshape(sta_spatial,1,[]);

        % update waitbar
        %T=text_waitbar(T,cc/length(cell_ids));
    end

    % view each STA
    %for cc=1:num_rfs; figure(2);clf;imagesc(norm_image(reshape(full(stas_matrix(:,cc)),diff(roi_y)+1,diff(roi_x)+1,3)));title(num2str(cell_ids(cc)));pause; end
end




% ITERATE

%tic


% initialize RGC weights matrix
B = sparse(num_kernels,num_rfs);

% initialize storage variables
best_dll = zeros(num_iter,1);
best_kernel = zeros(num_iter,1);
added_cones = zeros(num_iter,9);
added_cones_mat = zeros(num_kern_y,num_kern_x,kern_c);
likelihood = zeros(num_iter,1);

% effective cone density
eff_q = q * (kernel_spacing^2);


% begin looping through iterations
for ii =1:num_iter

    % compute intermediate value (used in a couple computations below)
    temp_value = (stas_matrix - W*B)'*W;

    % compute delta_log_likelihood
    delta_log_likelihood =  cell_constants' * (full(temp_value).^2) ./ (2 * kernel_norms);
    
    % this calculation was used before each cell had a different constant:
    %delta_log_likelihood =    sum((full(temp_value).^2),1) ./ (2 * kernel_norms);
    
    % store the dll of the first iteration for plotting below
    if ii==1;first_dll = delta_log_likelihood;end
    
    
    % compute the term in delta prior for cone repulsion
    % initialize
    repulsion_force = zeros(num_kern_y,num_kern_x,kern_c);
    % for each pair of cone types...
    for new_cone_c = 1:kern_c
        for cc=1:kern_c
            % compute the repulsion from other cones
            repulsion_force(:,:,new_cone_c) = repulsion_force(:,:,new_cone_c) + ...
                imfilter(added_cones_mat(:,:,cc),log(dist_filt{new_cone_c,cc}));
        end
    end
    
    % sum up terms in delta_prior
    delta_prior = log(eff_q) - log(1-eff_q) + reshape(permute(repulsion_force,[3 1 2]),1,[]);
    
    
    % arbitrary scale factor!
    delta_prior_ = delta_prior; % save a copy with no multiplier for plotting below
    delta_prior = delta_prior * magic_number;


    % find best kernel to weights matrix
    [best_dll(ii), best_kernel(ii)] = max(delta_log_likelihood + delta_prior);

    % if it is less than 0, finish
    if best_dll(ii) < 0
        break
    end
    
    % add weights to it in each RGC
    B(best_kernel(ii),:) = temp_value(:,best_kernel(ii))' ./ kernel_norms(best_kernel(ii));
    
    % note the location of the new cone in a list of coordinates
    bk_c = mod(best_kernel(ii)-1,kern_c) + 1;
    bk_y = mod( ceil(best_kernel(ii)/kern_c) - 1, num_kern_y) + 1;
    bk_x = mod( ceil(best_kernel(ii)/(kern_c*num_kern_y)) - 1, num_kern_x) + 1;
    bk_real_x = kernel_x(ceil(best_kernel(ii)/kern_c));
    bk_real_y  = kernel_y(ceil(best_kernel(ii)/kern_c));
    added_cones(ii,:) = [bk_x bk_y bk_c bk_real_x bk_real_y...
        best_dll(ii) delta_log_likelihood(best_kernel(ii)) delta_prior(best_kernel(ii)) best_kernel(ii)];
    
    % and in a matrix
    added_cones_mat(bk_y,bk_x,bk_c) = 1;
    
    
    % compute likelihood
    % used to evaluate whether delta-likelihood is correct
    if 0
        % compute fits
        fit_stas = W*B;
        % initialize to 0
        likelihood_ = 0;
        % add in contribution of each cell
        for cc=1:length(cell_ids)
            likelihood_ = likelihood_ + n_spikes(cc)*(fit_b(cc) + fit_stas(:,cc)'*raw_stas_matrix(:,cc)) + ...
                num_frames * exp(fit_b(cc) + 0.5 * stim_variance * (fit_stas(:,cc)'*fit_stas(:,cc)));
        end
        % store
        likelihood(ii) = likelihood_;
        
        % show 
        if ii > 1
            disp(log(likelihood(ii)) - log(likelihood(ii-1)))
        end
    end

    
    
    % plot update
    if 1 && (ii == 1 || ii > 35 || 1) && make_plots
        figure(plot_fig);clf;set(plot_fig,'color','w')
        la = []; % to store axes to be linked
        
        % dll
        subplot(2,3,6);
        plot(best_dll)
        title(sprintf('iter %d, best = %0.1f, magic number = %0.1f',ii,best_dll(ii),magic_number))
        
        % remaining delta_log_likelihood
        subplot(2,3,1); la = [la gca];
        imagesc((norm_image(   permute(reshape(delta_log_likelihood,kern_c,num_kern_y,num_kern_x),[2 3 1])   )-0.5).^0.4); axis image;
        title(datarun.names.nickname)
        
        % prior
        subplot(2,3,2); la = [la gca];
        imagesc((norm_image(   permute(reshape(exp(delta_prior_),kern_c,num_kern_y,num_kern_x),[2 3 1])   )-0.5).^0.4); axis image;
        
        % prior + dll
        subplot(2,3,3); la = [la gca];
        %imagesc((norm_image(   permute(reshape(delta_prior + delta_log_likelihood,kern_c,num_kern_y,num_kern_x),[2 3 1])   )-0.5).^0.4)
        imagesc((norm_image(...
            permute(reshape((delta_prior + delta_log_likelihood).*((delta_prior + delta_log_likelihood)>0),kern_c,num_kern_y,num_kern_x),[2 3 1])...
            )-0.5).^0.4)
        axis image;
        
        % all cones found so far
        subplot(2,3,4); la = [la gca];
        imagesc((norm_image(   permute(reshape(sum(abs(full(B)),2),kern_c,num_kern_y,num_kern_x),[2 3 1])   )-0.5).^0.4); axis image;
        
        % vector representation of cones
        subplot(2,3,5); la = [la gca];
        image(zeros(num_kern_y,num_kern_x)); hold on;axis image
        %figure(5);clf; image(zeros(datarun.stimulus.field_height,datarun.stimulus.field_width)); hold on;axis image 
        imagesc((norm_image(   permute(reshape(first_dll,kern_c,num_kern_y,num_kern_x),[2 3 1])   )-0.5).^0.4); axis image;
        for cc=1:kern_c
            cns = added_cones(:,3) == cc;
            figure(plot_fig);subplot(2,3,5);
            plot(added_cones(cns,1),added_cones(cns,2),'.','Color',kernel_plot_colors(cc,:))
            %figure(1);
            %plot(added_cones(cns,4),added_cones(cns,5),'o','Color',kernel_plot_colors(cc,:),'MarkerSize',15)
        end
        
        linkaxes(la)
        
        %if ii > 15; pause;end
    end
    drawnow

end

%toc


% VIEW FITS
if make_plots && 1
    
    fits = W*B;
    
    start_index=1; index_min=1; index_max=num_rfs;

    figure(plot_fig + 1);clf
    slider = make_loop_slider_list(start_index,index_min,index_max);
    while 1
        cc = round(get(slider,'Value'));
        cell_id = cell_ids(cc);
        
        % STA
        sta = reshape(full(stas_matrix(:,cc)),diff(roi_y)+1,diff(roi_x)+1,3);
        subplot(131);imagesc(norm_image(sta));axis image;title(sprintf('sta (%d)',cell_id))
        
        % fit
        fit = reshape(full(fits(:,cc)),diff(roi_y)+1,diff(roi_x)+1,3);
        subplot(132);imagesc(norm_image(fit));axis image;title('fit')
        
        % residual
        subplot(133);imagesc(norm_image(sta- fit));axis image;title('sta - fit')
        
        linkaxes([subplot(131) subplot(132) subplot(133)])
        
        uiwait;
    end

end

output_info.added_cones = added_cones;
output_info.first_dll = first_dll;

found_cones{1} = output_info;

%,first_dll,W,kernel_norms,num_kern_y,num_kern_x







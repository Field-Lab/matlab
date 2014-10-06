function bcf = bayesian_cone_finding_multicore(datarun,bcf_params)
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


% CUT UP INTO SMALL REGIONS

% count ROIs needed
num_roi_x = ceil( diff(relevant_region_x) / (roi_x_size - 2*padding_x) );
num_roi_y = ceil( diff(relevant_region_y) / (roi_y_size - 2*padding_y) );

% initialize variable
rois_x = zeros(num_roi_x*num_roi_y,2);
rois_y = rois_x;

% specify each
rr = 1;
for xx = 1:num_roi_x
    for yy = 1:num_roi_y
        start_x = relevant_region_x(1) + (xx-1)*(roi_x_size - 2*padding_x);
        start_y = relevant_region_y(1) + (yy-1)*(roi_y_size - 2*padding_y);
        rois_x(rr,:) = [start_x start_x + roi_x_size - 1];
        rois_y(rr,:) = [start_y start_y + roi_y_size - 1];
        rr = rr + 1;
    end
end

%rois_x = rois_x(63,:);rois_y = rois_y(63,:);
% show maximum range
fprintf('\nY range: %d to %d, X range: %d to %d\n\n',min(rois_x(:,1)),max(rois_x(:,2)),min(rois_y(:,1)),max(rois_y(:,2)))


% ANALYZE EACH SEPARATELY
if ~exist('loop_cones','var') || length(loop_cones) < size(rois_x,1)
    
    make_plots = 0;

    % remove clutter
    loop_mode = true;
    clear W

    % initialize storage
    loop_cones = cell(0);
    dlls = cell(0);

    figure(cones_fig);clf;
    image(0.5*ones(datarun.stimulus.field_height,datarun.stimulus.field_width,3));axis image; hold on;
    cone_axes = gca;
    title(datarun.names.nickname)
    drawnow

    % note start
    start_time_all = clock;
    
    % initialize W as empty
    W = [];
    kernel_norms = [];

    
    %% Multicore enabled

    
    if isempty(W)

        clear distance_prior
        if length(fieldnames(bcf_params.kernel_colors)) == 1  % BW stimulus
            distance_prior = C_C;
        else
            distance_prior(:,:,1) = [LM_MM; LM_MM; LM_S];
            distance_prior(:,:,2) = [LM_MM; LM_MM; LM_S];
            distance_prior(:,:,3) = [LM_S; LM_S; S_S];
        end

        kern_c = size(distance_prior,3);
        
        roi_x = rois_x(1,:);
        roi_y = rois_y(1,:);
        
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

    if isempty(W)


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
    
    
    cd ~/test/
    datarun.stas = rmfield(datarun.stas, 'java_sta');
    save datarun datarun
    
    bcf_cell = cell(size(rois_x,1),1);
    num_chunks = size(rois_x,1)
    for rr = 1:size(rois_x,1)
        bcf_stuff.bcf_params = bcf_params;
        bcf_stuff.roi_x = rois_x(rr,:);
        bcf_stuff.roi_y = rois_y(rr,:);
        bcf_stuff.W = W;
        bcf_stuff.kernel_norms = kernel_norms;
        bcf_stuff.datarun_path = '~/test/datarun.mat';
        bcf_cell{rr} = bcf_stuff;
    end    
    
    
    settings.multicoreDir = '~/multicore-matlab/';
    settings.useWaitbar = true;
    
    found_cones = startmulticoremaster(@parallel_bayesian_cone_finding, bcf_cell, settings);
    
    for rr = 1:size(rois_x,1)
        tmp = found_cones{rr};
        loop_cones{rr} = tmp.output_info.added_cones;
        dlls{rr} = tmp.output_info.first_dll;
    end
        

    % add cones to plot
    if ~isempty(bcf_info.added_cones)
        for cc=1:size(kernel_plot_colors,1)
            cns = added_cones(:,3) == cc;
            plot(added_cones(cns,4),added_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'Parent',cone_axes) %#ok<COLND>
        end
        drawnow
    end
    %%
    
%     % go through each subset
%     for rr = 1:size(rois_x,1)
%         start_time_one = clock;
% 
%         % set up ROI
%         roi_x = rois_x(rr,:);
%         roi_y = rois_y(rr,:);
% 
%         % remove clutter
%         clear stas_matrix
% 
%         % analysis
%         [added_cones,first_dll,W,kernel_norms,num_kern_y,num_kern_x] = bayesian_cone_finding(datarun,bcf_params,roi_x,roi_y,W,kernel_norms);
% 
%         % save results
%         loop_cones{rr} = added_cones;
%         dlls{rr} = first_dll;
%         
%         % add cones to plot
%         if ~isempty(added_cones)
%             for cc=1:size(kernel_plot_colors,1)
%                 cns = added_cones(:,3) == cc;
%                 plot(added_cones(cns,4),added_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'Parent',cone_axes) %#ok<COLND>
%             end
%             drawnow
%         end
% 
%         % note duration
%         fprintf('\nfinished iteration %d of %d in %0.1f sec\n',rr,size(rois_x,1),etime(clock,start_time_one))
%     end
    fprintf('\nfinished everything in %0.1f sec\n\n\n',etime(clock,start_time_all))

end



% CONCATENATE DELTA-LOG-LIKELIHOOD IMAGES and ACCUMULATE CONES
if 1

    % initialize
    dll = zeros((datarun.stimulus.field_height-1)/kernel_spacing+1,(datarun.stimulus.field_width-1)/kernel_spacing+1,size(kernel_plot_colors,1));
    all_added_cones = [];

    % go through each subset
    for rr = 1:size(rois_x,1)
        
        % likelihood
        
        % identify where this region fits
        x_range = [rois_x(rr,1) + padding_x - 1   rois_x(rr,2) - padding_x - kernel_spacing]/kernel_spacing;
        y_range = [rois_y(rr,1) + padding_y - 1   rois_y(rr,2) - padding_y - kernel_spacing]/kernel_spacing;
        x_range = round(x_range);
        y_range = round(y_range);

        % reshape dll
        temp = permute(reshape(dlls{rr},size(kernel_plot_colors,1),num_kern_y,num_kern_x),[2 3 1]);

        % put it in
        dll(y_range(1):y_range(2),x_range(1):x_range(2),:) = ...
            temp(1+(padding_y-1)/kernel_spacing:end-padding_y/kernel_spacing-1,...
            1+(padding_x-1)/kernel_spacing:end-padding_x/kernel_spacing-1,:);
        
        
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


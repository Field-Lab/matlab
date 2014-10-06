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


% CUT UP INTO SMALL REGIONS

% count ROIs needed
num_roi_x = ceil( diff(relevant_region_x) / (roi_x_size - 2*padding_x) );
num_roi_y = ceil( diff(relevant_region_y) / (roi_y_size - 2*padding_y) );

% initialize variable
rois_x = zeros(num_roi_x*num_roi_y,2);
rois_y = rois_x;
num_rfs_roi=zeros(num_roi_x*num_roi_y,1);

% specify each
rr = 1;
for xx = 1:num_roi_x
    for yy = 1:num_roi_y
        start_x = relevant_region_x(1) + (xx-1)*(roi_x_size - 2*padding_x);
        start_y = relevant_region_y(1) + (yy-1)*(roi_y_size - 2*padding_y);
        rois_x(rr,:) = [start_x start_x + roi_x_size - 1];
        rois_y(rr,:) = [start_y start_y + roi_y_size - 1];
        
        
        roi = false(datarun.stimulus.field_height,datarun.stimulus.field_width);
        % set a block of values to 1
        roi(rois_y(rr,1):rois_y(rr,2),rois_x(rr,1):rois_x(rr,2)) = true;        
        cell_ids = datarun.cell_ids(get_cell_indices_roi(datarun,cone_finding_cell_spec,roi,'cell_locations','marks'));
        num_rfs_roi(rr) = length(cell_ids);
        
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

    if ~isempty(cones_fig)
        figure(cones_fig);clf;
        image(0.5*ones(datarun.stimulus.field_height,datarun.stimulus.field_width,3)); axis image; hold on;
        cone_axes = gca;
        title(datarun.names.nickname)
        drawnow
    end

    % note start
    start_time_all = clock;
    
    % initialize W as empty
    W = [];
    kernel_norms = [];

    % go through each subset
    
    my_rois=find(num_rfs_roi);
    
    for rr = my_rois' %1:size(rois_x,1)
        start_time_one = clock;

        % set up ROI
        roi_x = rois_x(rr,:);
        roi_y = rois_y(rr,:);
 
        % remove clutter
        clear stas_matrix

        % analysis
        [added_cones,first_dll,W,kernel_norms,num_kern_y,num_kern_x] = bayesian_cone_finding(datarun,bcf_params,roi_x,roi_y,W,kernel_norms);
      
        % save results
        loop_cones{rr} = added_cones;
        dlls{rr} = first_dll;
        
        % add cones to plot
        if ~isempty(added_cones)
            for cc=1:size(kernel_plot_colors,1)
                cns = added_cones(:,3) == cc;
                
                if ~isempty(cones_fig)
                    plot(added_cones(cns,4),added_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'Parent',cone_axes) %#ok<COLND>
                end
            end
            drawnow
        end

        % note duration
        fprintf('\nfinished iteration %d of %d in %0.1f sec\n',rr,size(rois_x,1),etime(clock,start_time_one))
    end
    fprintf('\nfinished everything in %0.1f sec\n\n\n',etime(clock,start_time_all))

end



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


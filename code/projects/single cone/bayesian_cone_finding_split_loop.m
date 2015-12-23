function bcf = bayesian_cone_finding_split_loop(datarun, bcf_params, varargin)
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
% 2010-10  phli, copied from BCF_loop, with added options for parallelization
%                should be completely back compatible with BCF_loop, so
%                ideally will copy this into that namespace once absolutely
%                sure no problems.
%

opts = inputParser;
opts.addParamValue('plot', true);
opts.addParamValue('threads', 1);
opts.addParamValue('thread',  1);
opts.addParamValue('thread_roi_set_ratios', []);
opts.parse(varargin{:});
opts = opts.Results;


% EXPAND PARAMETERS
fn = fieldnames(bcf_params);
for ff=1:length(fn)
    eval(sprintf('%s = bcf_params.%s;',fn{ff},fn{ff}))
end


% CUT UP INTO SMALL REGIONS

% count ROIs needed
num_roi_x = ceil( diff(relevant_region_x) / (roi_x_size - 2*padding_x) );
num_roi_y = ceil( diff(relevant_region_y) / (roi_y_size - 2*padding_y) );
num_rois = num_roi_x * num_roi_y;

% initialize variables
rois_x = zeros(num_rois, 2);
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


% pick subset of rois for this parallel thread
rois_per_set = num_rois ./ opts.threads;
if isempty(opts.thread_roi_set_ratios)
    opts.thread_roi_set_ratios = 0:1/opts.threads:1;
end
set_edges = round(opts.thread_roi_set_ratios .* num_rois);
set_starts = set_edges(1:opts.threads) + 1;
set_ends   = [set_edges(2:opts.threads) num_rois];
start_roi = set_starts(opts.thread);
end_roi   = set_ends(opts.thread);
num_rois_this_thread = end_roi - start_roi + 1;
disp(['This thread will run rois ' num2str(start_roi) ' through ' num2str(end_roi) ' out of ' num2str(num_rois) ' total']);


% ANALYZE EACH SEPARATELY
if ~exist('loop_cones','var') || length(loop_cones) < size(rois_x,1)
    
    make_plots = 0;

    % remove clutter
    loop_mode = true;
    clear W

    % initialize storage
    loop_cones = cell(0);
    dlls = cell(0);

    if opts.plot
        figure(cones_fig);clf;
        image(0.5*ones(datarun.stimulus.field_height,datarun.stimulus.field_width,3));axis image; hold on;
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
    rr_this_thread = 1;
    for rr = start_roi:end_roi
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
        if opts.plot && ~isempty(added_cones)
            for cc=1:size(kernel_plot_colors,1)
                cns = added_cones(:,3) == cc;
                plot(added_cones(cns,4),added_cones(cns,5),'.','Color',kernel_plot_colors(cc,:),'Parent',cone_axes);
            end
            drawnow
        end

        % note duration
        fprintf('\nfinished iteration %d of %d in %0.1f sec\n', rr_this_thread, num_rois_this_thread, etime(clock,start_time_one))
        rr_this_thread = rr_this_thread + 1;
    end
    fprintf('\nfinished %d iterations in %0.1f sec\n\n\n', num_rois_this_thread, etime(clock,start_time_all))

end



% CONCATENATE DELTA-LOG-LIKELIHOOD IMAGES and ACCUMULATE CONES

% initialize
dll = zeros((datarun.stimulus.field_height-1)/kernel_spacing+1,(datarun.stimulus.field_width-1)/kernel_spacing+1,size(kernel_plot_colors,1));
all_added_cones = [];

% go through each subset
for rr = start_roi:end_roi

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
            (loop_cones{rr}(:,4) >= roi_x(1)) & (loop_cones{rr}(:,4) < roi_x(2)) & ...
            (loop_cones{rr}(:,5) >= roi_y(1)) & (loop_cones{rr}(:,5) < roi_y(2)) ...
            , :);

        % keep this in growing list
        all_added_cones = [all_added_cones ; added_cones]; %#ok<AGROW>
    end

    % add boundaries of the roi
    %plot([roi_x roi_x(2) roi_x(1) roi_x(1)],[roi_y(1) roi_y(1) roi_y(2) roi_y(2) roi_y(1)],'w')
end



% package up the results
bcf.all_added_cones = all_added_cones;
bcf.dll = dll;
bcf.rois_x = rois_x;
bcf.rois_y = rois_y;
% temp function -- delete


nbr_radius = 3;
%roix=155:185;%167:176;
%roiy=154:185;%157:165;
%roix = round(cone_centers(nn,1))+[-nbr_radius:nbr_radius];
%roiy = round(cone_centers(nn,2))+[-nbr_radius:nbr_radius];

roix = max(round(cone_centers(nn,1)-nbr_radius),1):min(round(cone_centers(nn,1)+nbr_radius),datarun.stimulus.field_width);
roiy = max(round(cone_centers(nn,2)-nbr_radius),1):min(round(cone_centers(nn,2)+nbr_radius),datarun.stimulus.field_height);

gauss_params = {'center_radius',0.75,'effective_radius',3};


% identify RGCs in the ROI

if ~exist('roicells','var')

    % ROI in graphic form
    roi=zeros(datarun.stimulus.field_height,datarun.stimulus.field_width);
    roi(roiy,roix)=1;

    % reshape to compare to all_sig_stixels
    roi = reshape(roi,[],1);

    % find cells in there
    roicells = find(any(all_sig_stixels&repmat(roi,1,size(all_sig_stixels,2)),1));

    % plot to verity
    if 0
        for cc=1:length(roicells)
            figure(2);clf;
            imagesc(reshape(all_sig_stixels(:,roicells(cc)),datarun.stimulus.field_height,datarun.stimulus.field_width))
            pause
        end
    end
end


% get a sensitivity map of each such RGC in this region

if ~(exist('rfs','var') && exist('roicenters','var'))

    rfs = zeros(length(roiy),length(roix),length(roicells));
    cell_indices = get_cell_indices(datarun,cell_ids);

    % go through each found cell
    for cc=1:length(roicells)

        % get the rf
        rf = get_rf(datarun,datarun.cell_ids(cell_indices(roicells(cc))));

        % get the roi
        rf = rf(roiy,roix,:)/robust_std(reshape(rf,[],1));

        % combine across color channels
        rfs(:,:,cc) = sum(rf,3);
    end



    % get cone center points, translate to this region

    % get center points
    %cp=datarun.cones.centers;
    cp=cone_centers;

    % find ones in this region
    center_indices = find(cp(:,1)>min(roix)-0.5 & cp(:,1)<max(roix)+0.5 & cp(:,2)>min(roiy)-0.5 & cp(:,2)<max(roiy)+0.5);

    % keep them
    roicenters=cp(center_indices,1:2);

    % subtract to fit in ROI
    roicenters = roicenters - repmat([min(roix) min(roiy)],size(roicenters,1),1) + 1;

    % plot RGCs with center points overlaid
    if 1
        for cc=1:length(roicells)
            figure(2);clf;colormap gray
            imagesc(rfs(:,:,cc))
            hold on; plot(roicenters(:,1),roicenters(:,2),'.r','MarkerSize',25)
            pause
        end
    end

end


% fit gaussians
if ~exist('fit_cone_info','var')

    % prepare initial cone parameters
    cone_initial = roicenters;
    
    % estimate heights for each cone in each RF
    for cc=1:size(roicenters,1)
        % make little gaussian of the cone
        conerf = make_gaussian('x_size',size(rfs,2),'y_size',size(rfs,1),...
            'center',roicenters(cc,:),gauss_params{:},'normalize','sum');
        
        % estimate height for each RF
        for rr=1:length(roicells)
            cone_initial(cc,rr+2) = reshape(rfs(:,:,rr),[],1)'/reshape(conerf,[],1)';
        end
    end

    % fit them
    [fit_cone_info,rf_fit] = fit_cones_in_many_RFs(rfs,cone_initial,'figure',3,...
        'cone_kernel',struct('type','dog',gauss_params{:}));

end

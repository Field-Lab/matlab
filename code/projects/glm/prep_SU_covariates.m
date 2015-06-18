function SU_cov_vec = prep_SU_covariates(pooling_filter, fitmovie, ROIcoords)

stimsize_x = size(fitmovie, 1);
stimsize_y = size(fitmovie, 2);
bins = size(fitmovie, 3);

SU_cov_vec = zeros(9,bins);

% loop through 9 SU pixels
SU_idx = 0;
for j_SU = 1:3
    for i_SU = 1:3
        SU_idx = SU_idx + 1;
        cov = zeros(1,bins);
        dist_from_SU_center = [i_SU j_SU]-[2 2]; % find the offset from the pooling filter weight on that subunit to the location of the SU pixel
        
        % loop through all possible SU locations
        for i_PF = 1:size(pooling_filter, 1)
            for j_PF = 1:size(pooling_filter, 2)
                
                % find the location of that subunit pixel in the stimulus,
                % for the pixel in the pooling filter
                stim_idx = [ROIcoords.xvals(i_PF)+dist_from_SU_center(1), ROIcoords.yvals(j_PF)+dist_from_SU_center(2)];
                
                % As long as it is a valid location, within the stimulus
                if stim_idx(1)>=1 && stim_idx(1)<=stimsize_x && stim_idx(2)>=1 && stim_idx(2)<=stimsize_y
                    
                    % Find the value of the pooling filter at that location
                    pooling_filter_weight = pooling_filter(i_PF, j_PF);
                    
                    % find the value of the stimulus on the subunit pixel
                    pixel_values = squeeze(fitmovie(stim_idx(1), stim_idx(2), :));
                    
                    % add to cov vec
                    cov = cov + pooling_filter_weight*double(pixel_values');
                end
            end
        end
        
        SU_cov_vec(SU_idx, :) = cov;
    end
   
end
end
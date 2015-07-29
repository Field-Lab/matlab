function [SU_cov_vec, pooling_weights] = prep_SU_covariates(pooling_filter, fitmovie, ROIcoords, inputstats, timefilter)

% might need to move this to before SU? NB
fitmoviestats.mean     =  inputstats.mu_avgIperpix;
fitmoviestats.span     =  inputstats.range;
fitmoviestats.normmean =  inputstats.mu_avgIperpix / inputstats.range;
stim   = double(fitmovie) / double(fitmoviestats.span)-double(fitmoviestats.normmean);

stimsize_x = size(fitmovie, 1);
stimsize_y = size(fitmovie, 2);
bins = size(fitmovie, 3);
n_locations = size(pooling_filter, 1)*size(pooling_filter, 2);
SU_cov_vec = zeros(9,n_locations,bins);
pooling_weights = zeros(n_locations,1);

% loop through 9 SU pixels
SU_idx = 0;
for j_SU = 1:3
    for i_SU = 1:3
        SU_idx = SU_idx + 1;
        dist_from_SU_center = [i_SU j_SU]-[2 2]; % find the offset from the pooling filter weight on that subunit to the location of the SU pixel
        
        % loop through all possible SU locations
        loc_idx = 0;
        for i_PF = 1:size(pooling_filter, 1)
            for j_PF = 1:size(pooling_filter, 2)
                loc_idx = loc_idx+1;
                % find the location of that subunit pixel in the stimulus,
                % for the pixel in the pooling filter
                stim_idx = [ROIcoords.xvals(i_PF)+dist_from_SU_center(1), ROIcoords.yvals(j_PF)+dist_from_SU_center(2)];
                
                % As long as it is a valid location, within the stimulus
                if stim_idx(1)>=1 && stim_idx(1)<=stimsize_x && stim_idx(2)>=1 && stim_idx(2)<=stimsize_y
                    
                    % Find the value of the pooling filter at that location
                    pooling_weights(loc_idx) = pooling_filter(i_PF, j_PF);
                    
                    % find the value of the stimulus on the subunit pixel
                    pixel_values = squeeze(stim(stim_idx(1), stim_idx(2), :));
                    if ~(timefilter==0)
                        stim_temp = conv(pixel_values, timefilter, 'full');
                        pixel_values = stim_temp(1:bins);
                    end
                    
                    % add to cov vec
                    SU_cov_vec(SU_idx, loc_idx, :) = pixel_values;
                end
            end
        end
    end
    
end
end
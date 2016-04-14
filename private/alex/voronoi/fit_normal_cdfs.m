function loglikratio = fit_normal_cdfs(inputs, spikes, center_cones)

if length(spikes)==size(inputs,2)
    spike_rate = spikes;
else
    
    spikes(spikes>size(inputs,2)) = [];
    spikes_tmp = spikes;
    spike_rate=zeros(size(inputs,2),1);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        spikes_tmp(ia)=[];
    end
    clear spikes_tmp
end

%% NORMAL CDF

if 0
    loglikratio =[];
    
    for cone1 = 1:length(center_cones)
        cone1
        other_cones = 1:length(center_cones);
        other_cones(cone1)=[];
        tic
        for cone2 = other_cones
            
            % get the input ranges for conditioning cone
            mean2 = mean(inputs(cone2, :));
            std2 = std(inputs(cone2, :))/2;
            
            inds_min = find(inputs(cone2, :) < (mean2-std2));
            %         inds_mid = find(inputs(cone2, :) > (mean2-std2) & inputs(cone2, :) < (mean2+std2));
            inds_max = find(inputs(cone2, :) > (mean2+std2));
            
            if length(inds_max) > length(inds_min)
                inds_max = inds_max(1:length(inds_min));
            else
                inds_min = inds_min(1:length(inds_max));
            end
            
            % re-sampling
            for sample = 1:9
                resampled_inds=round(rand(length(inds_max),1)*length(inds_max));
                resampled_inds(resampled_inds==0)=1;
                resampled_inds(resampled_inds>length(inds_max)) = length(inds_max);
                resampled_inds_min = inds_min(resampled_inds);
                
                resampled_inds=round(rand(length(inds_max),1)*length(inds_max));
                resampled_inds(resampled_inds==0)=1;
                resampled_inds(resampled_inds>length(inds_max)) = length(inds_max);
                resampled_inds_max = inds_max(resampled_inds);
                
                cone_input_min = inputs(cone1, resampled_inds_min);
                firing_rate_min = spike_rate(resampled_inds_min);
                cone_input_max = inputs(cone1, resampled_inds_max);
                firing_rate_max = spike_rate(resampled_inds_max);
                
                % start points - least square
                [start_points_x, start_points_y] = ...
                    fit_norm_cdf_lsq(cone_input_min',firing_rate_min, cone_input_max', firing_rate_max);
                
                % maximize loglik
                [mllparams_x, fval_x, mllparams_y, fval_y] = ...
                    fit_norm_cdf_mll(cone_input_min',firing_rate_min, cone_input_max', firing_rate_max, ...
                    start_points_x, start_points_y);
                
                loglikratio(cone1, cone2, sample) = fval_x-fval_y;
            end
        end
        toc
    end
end


%% NORMAL CDF 2-D
loglikratio =[];

n_inps = size(inputs,2);

if 1
    for cone1 = 1:length(center_cones)-1
        cone1
        tic
        for cone2 = (cone1+1):length(center_cones)
            cone2
            % re-sampling
            tic
            for sample = 1:9                
                resampled_inds=round(rand(n_inps,1)*n_inps);
                resampled_inds(resampled_inds==0)=1;
                resampled_inds(resampled_inds>n_inps)=n_inps;
                [start_points_x, resn1, start_points_y, resn2] = ...
                    fit_norm_cdf_2d_lsq(inputs(cone1, resampled_inds),inputs(cone2, resampled_inds), ...
                    spike_rate(resampled_inds));
                
                % maximize loglik
                [mllparams_x, fval_x, mllparams_y, fval_y] = ...
                    fit_norm_cdf_2d_mll(inputs(cone1, resampled_inds),inputs(cone2, resampled_inds), ...
                    spike_rate(resampled_inds), start_points_x, start_points_y);
                
%                 save_2d_plots(inputs(cone1, resampled_inds)', inputs(cone2, resampled_inds)',...
%                     spike_rate(resampled_inds), mllparams_x, mllparams_y, ...
%                     center_cones(cone1), center_cones(cone2), sample)
                
                loglikratio(cone1, cone2, sample) = fval_x-fval_y;
            end
            toc
        end
        toc
    end
end


% for cone1 = 1:length(center_cones)-1
%     cone1
%     for cone2 = (cone1+1):length(center_cones)
%         save_2d_plots(inputs(cone1, :)', inputs(cone2, :)',...
%             spike_rate, [], [], center_cones(cone1), center_cones(cone2), [])
%     end
% end
%
%
%
% tmp = loglikratio(4, 5, :);
% tmp=reshape(tmp,3,3);
% mean(tmp(:))/std(tmp(:))
%
%
%
%
% all_array = zeros(length(center_cones)*3);
% for cone1 = 1:length(center_cones)
%
%     for cone2=1:length(center_cones)
%         tmp = loglikratio(cone1, cone2, :);
%         tmp=reshape(tmp,3,3);
%         all_array(cone1*3-2:cone1*3,cone2*3-2:cone2*3) = tmp;
%     end
% end
% figure
% tmp = all_array;
% % tmp(tmp>160) = 160;
% imagesc(tmp);
% set(gca, 'xtick', 2:3:length(center_cones)*3,'xticklabel',int2str(center_cones))
% set(gca, 'ytick', 2:3:length(center_cones)*3,'yticklabel',int2str(center_cones))
% xlabel('conditioning cones')
% ylabel('referent cones')
% set(gca,'dataaspectratio', [1 1 1])
% hold on
% line([0, length(center_cones)*3+1], [0, length(center_cones)*3+1], 'color', 'k')
%
%

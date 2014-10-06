function proj_struct = spike_compute_projection(projection_type,proj_struct,dataset);

% compute projection of a projection structure (proj_struct)
% here only two fields are used:
%
% proj_struct.spike_waveforms
%            .spike_projections

NUM_SPIKES_FOR_COVARIANCE = 10000;
num_spikes_to_use = min(size(proj_struct.spike_waveforms,1),NUM_SPIKES_FOR_COVARIANCE);
%waveforms_for_proj = proj_struct.spike_waveforms(1:num_spikes_to_use,proj_struct.spike_waveform_points_to_use-1);
%all_waveforms = proj_struct.spike_waveforms(:,proj_struct.spike_waveform_points_to_use-1);

switch projection_type
    case 'pca'
        %perform PCA using matlab's built-in function
        clear dataset;
        [coef, pc_projection] = princomp(single(proj_struct.spike_waveforms(1:num_spikes_to_use,proj_struct.spike_waveform_points_to_use-1)));
        proj_struct.spike_projections = single(proj_struct.spike_waveforms(:,proj_struct.spike_waveform_points_to_use-1) * coef);
        %caption = 'PC';
        
    case 'wavelet'
        try
        %code to do wavelet decomposition using the wavelet toolbox
        %this implementation largely from wave_clus framework
        clear dataset;
        
        
        %fewer input?
        %scales parameter wrong?
        scales = 4;
        
        points = proj_struct.spike_waveform_points_to_use-1;
        
        cc=zeros(num_spikes_to_use,length(points));
        
        for i=1:num_spikes_to_use              % Wavelet decomposition
            c = wavedec(single(proj_struct.spike_waveforms(i,proj_struct.spike_waveform_points_to_use-1)),scales,'haar');
            cc(i,points)=c(points);
        end
        
        
        sd = zeros(length(points));
        
        for i=points                         % KS test for coefficient selection   
            thr_dist = std(cc(:,i)) * 3;
            thr_dist_min = mean(cc(:,i)) - thr_dist;
            thr_dist_max = mean(cc(:,i)) + thr_dist;
            aux = cc(find(cc(:,i)>thr_dist_min & cc(:,i)<thr_dist_max),i);
       
            if length(aux) > 10;
                [ksstat]=test_ks(aux);
                sd(i)=ksstat;
            else
                sd(i)=0;
            end
        end
        
        
        [max ind]=sort(sd);
        coef(1:length(points))=ind(length(points):-1:1);
        
        
        %[max coef]=sort(sd, 'descend');
        
        
        
        %figure(9); plot(proj_struct.spike_waveforms(:,find(coef == max(1))),proj_struct.spike_waveforms(:,find(coef == max(2))),'.'); 
        %proj_spike.spike_projections =
        %zeros(num_spikes_to_use,proj_struct.spike_waveform_points_to_use);
        %for i=1:num_spikes_to_use
        %proj_struct.spike_projections(i,:) = single(proj_struct.spike_waveforms(i,proj_struct.spike_waveform_points_to_use-1) .* coef);
        %end
        
        
        proj_spike.spike_projections = zeros(num_spikes_to_use,length(points));
        
        for i=1:num_spikes_to_use
            %for j=points
                proj_struct.spike_projections(i,:)=cc(i,coef(:));
            %end
        end
        
        
        
        disp('done')
        
        
        catch
            err = lasterror;
            err.stack.file
            err.stack.name
            err.stack.line
        end
        
    case {'nwpca-center','nwpca-all','lda', 'nnwpca'}
        %perform LDA using spike_lda_on_dataset
        %99 is the figure number
        LDvectors = single(spike_compute_lda_projection(dataset, single(proj_struct.spike_waveforms(1:num_spikes_to_use,proj_struct.spike_waveform_points_to_use-1)), projection_type,99));
        clear dataset;
        proj_struct.spike_projections = proj_struct.spike_waveforms(:,proj_struct.spike_waveform_points_to_use-1) * LDvectors;
        
end

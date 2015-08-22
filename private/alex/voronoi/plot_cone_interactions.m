function plot_cone_interactions(inputs,sta, n_bins,spikes)

n_cones = size(inputs,1);


spikes(spikes>size(inputs,2)) = []; 
% get spike rate
spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp ia

spike_rate(1:size(sta,2)-1) = [];

% filter all inputs with STA
filt_inputs = zeros(n_cones, size(inputs,2)-size(sta,2)+1);
for i = 1:n_cones
    filt_inputs(i,:) = conv(inputs(i,:),sta(i,:),'valid');
end

bin_size=floor(size(filt_inputs,2)/10);

for cone1 = 1:n_cones
    
    cone1_sorted = sort(filt_inputs(cone1,:));

    figure
    set(gcf, 'position',[1253         679         560         420])
    for cone2 = 1:n_cones
        if cone1 ~= cone2
            tmp = sort(filt_inputs(cone2,:));
            cone2_pos_inputs = filt_inputs(cone2,:)>tmp(end-bin_size);
            cone2_neg_inputs = filt_inputs(cone2,:)<tmp(bin_size);
            
            cone1_cone2 = zeros(n_bins-1,2);
            for i = 1:n_bins-1
                cone1_inputs = filt_inputs(cone1,:)>=cone1_sorted(bin_size*(i-1)+1) & filt_inputs(cone1,:)<cone1_sorted(bin_size*i);
                cone1_cone2(i,1) = mean(spike_rate(cone1_inputs & cone2_pos_inputs)); % cone2 positive
                cone1_cone2(i,2) = mean(spike_rate(cone1_inputs & cone2_neg_inputs)); % cone2 negative
            end
            subplot(1,2,2)
            plot(cone1_cone2(:,1))
            hold on
            title('cone 2 positive')
            subplot(1,2,1)
            plot(cone1_cone2(:,2))
            hold on        
            title('cone 2 negative')
        end
    end    
    
end

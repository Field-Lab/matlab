% construct a model with several cones, all summed up linearly and passed
% through nonlinearity (only 1). Emulate cone-interaction plots.

% compare with another model, where each cone gets through own
% nonlinearity, then summed up linearly with others, and passes through
% final nonlinearity.
far_cones = 1:length(sel_cones.far);
real_center_cones = (1:length(sel_cones.real)) + length(sel_cones.far);
real_sur_cones = (1:length(sel_cones.sur)) + length(sel_cones.far) + length(sel_cones.real);


mstas = mean_sta; % sta for each cone

minputs = inputs(cones,1:frames2take); % old raw inputs
finputs = inputs(cones,frames2take+1:end);% new raw inputs

msr = spike_rate_all(1:frames2take); % old spike rate
fsr = spike_rate_all(frames2take+1:end); % new spike rate

% filter old inputs
mfiltered_inputs=zeros(size(minputs,1),size(minputs,2)-sta_params.length+1);  
for i=1:size(minputs,1)
    mfiltered_inputs(i,:)=conv(minputs(i,:),mstas(i,:),'valid');
end
msr = msr(sta_params.length:end);

% filter new inputs
ffiltered_inputs=zeros(size(finputs,1),size(finputs,2)-sta_params.length+1);  
for i=1:size(minputs,1)
    ffiltered_inputs(i,:)=conv(finputs(i,:),mstas(i,:),'valid');
end
fsr = fsr(sta_params.length:end);

figure
for i=1:16
    subplot(4,4,i)
    plot(mstas(i+20,:))
    axis([1 15 -0.16 0.05])
%     hist(ffiltered_inputs(i+20,:),50);
%     axis([-0.2 0.2 0 1500])
end

% find individual nonlinearities - OLD data
nbins = 8;
my_nl=zeros(nbins,length(real_center_cones));
partial_gs = zeros(size(mfiltered_inputs,2),length(real_center_cones));

% bins - common for all
all_filt_inputs = mfiltered_inputs(real_center_cones,:);
all_filt_inputs=all_filt_inputs(:)';
n_gs=floor(length(all_filt_inputs)/nbins);
tmp=sort(all_filt_inputs);
my_bins=tmp(1:n_gs:n_gs*(nbins+1));
my_bins(end)=my_bins(end)+1;

all_filt_inputs = mfiltered_inputs(real_center_cones,:);
for i=1:size(all_filt_inputs,1)
    tt = all_filt_inputs(i,:);
    other_cones = 1:size(all_filt_inputs,1);
    other_cones(i) = [];
    good_inputs = 1:size(all_filt_inputs,2);
    for j=1:length(other_cones)
        good_inputs = intersect(good_inputs,find(abs(all_filt_inputs(other_cones(j),:)) < -my_bins(1)/3));
    end
        
    for k=1:length(my_bins)-1
        tmp=find(tt>=my_bins(k) & tt<my_bins(k+1)); % find instances when GS had certain value
        tmp = intersect(tmp, good_inputs);
        my_nl(k,i)=mean(msr(tmp)); % find mean firing rate after this certain value        
        partial_gs(tmp, i) = my_nl(k,i);
    end
end

figure
plot(my_bins(1:end-1), nanmean(my_nl'))

figure
plot(my_bins(1:end-1), my_nl)

my_nl = nanmean(my_nl');


% MODEL WITH SINGLETONS and no RGC nonlinearity

% find individual nonlinearities - OLD data
nbins = 8;
my_nl=zeros(nbins,length(real_center_cones));
partial_gs = zeros(size(mfiltered_inputs,2),length(real_center_cones));

% bins - common for all
all_filt_inputs = mfiltered_inputs(real_center_cones,:);
all_filt_inputs=all_filt_inputs(:)';
n_gs=floor(length(all_filt_inputs)/nbins);
tmp=sort(all_filt_inputs);
my_bins=tmp(1:n_gs:n_gs*(nbins+1));
my_bins(end)=my_bins(end)+1;

cnt = 1;
for i=real_center_cones
    tt = mfiltered_inputs(i,:);
    for k=1:length(my_bins)-1
        tmp=find(tt>=my_bins(k) & tt<my_bins(k+1)); % find instances when GS had certain value
        my_nl(k,cnt)=mean(msr(tmp)); % find mean firing rate after this certain value        
        partial_gs(tmp, cnt) = my_nl(k,cnt);
    end
    cnt=cnt+1;
end

figure
plot(my_bins(1:end-1), nanmean(my_nl'))

% mean nonlinearity for each subunits
my_nl = nanmean(my_nl');


% predict on NEW data

my_f_input = ffiltered_inputs(real_center_cones,:);
partial_gs = zeros(size(my_f_input,2),length(real_center_cones));
% pass through nonlinearity individually
for i=1:length(real_center_cones)
    tt = my_f_input(i,:);
   
    for k=1:length(my_bins)-1
        tmp=find(tt>=my_bins(k) & tt<my_bins(k+1)); % find instances when GS had certain value
        partial_gs(tmp, cnt) = my_nl(k);
    end
end

predicted_rate = sum(partial_gs')';

nbins = 100;
n_gs=floor(length(predicted_rate)/nbins);
tmp=sort(predicted_rate);
my_bins=tmp(1:n_gs:n_gs*(nbins+1));
my_bins(end)=my_bins(end)+1;
my_nl_rgc = zeros(nbins,1);
for k=1:length(my_bins)-1
    tmp=find(predicted_rate>=my_bins(k) & predicted_rate<my_bins(k+1)); % find instances when GS had certain value
    my_nl_rgc(k)=mean(fsr(tmp)); % find mean firing rate after this certain value
end
figure
plot(my_bins(1:end-1), my_nl_rgc)

a = fsr - predicted_rate;
mean(a)
predicted_rate = predicted_rate + mean(a);

sse = sum((predicted_rate-fsr).^2);
sst = sum((fsr-mean(fsr)).^2);
r2 = 1 - sse/sst;

figure
plot(predicted_rate)
hold on
plot(fsr)

% Estimate zero response bin in the model
for cone1 = cones_to_plot
    
    figure(cone1)
    set(gcf, 'Name',['Cone1 - ', int2str(cone1)])
    
    cone2_effect_pos = zeros(10, size(inputs,1));% bins, sta_trial, n of cones
    cone2_effect_neg = cone2_effect_pos;
    cone2_effect_aver = cone2_effect_pos;
    
    % bin cone1 input
    cone1_input = ffiltered_inputs(cone1,:);  % this is independent of sta
    n_gs=floor(length(cone1_input)/10);
    tmp=sort(cone1_input);
    my_bins=[min(cone1_input) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(cone1_input)+1;
    
    for cone2 = 1:size(inputs,1)
        if cone1 ~= cone2            
            
            % find 10% of most negative, positive, and near 0 input of cone2
            cone2_input = ffiltered_inputs(cone2,:);
            tmp = sort(cone2_input);
            cone2_neg_gs = find(cone2_input<tmp(ceil(length(tmp)/10)));
            cone2_pos_gs = find(cone2_input>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone2_aver_gs = find( cone2_input>tmp(aver+2*bord) & cone2_input<tmp(aver+4*bord));
            
            
            
            clear max_neg max_pos max_aver
            for i=1:length(my_bins)-1
                
                % find instances of cone1 input in certain range
                cone1_local_gs = find( cone1_input>=my_bins(i) & cone1_input<my_bins(i+1));
                
                % find common instances when cone2 was negative
                tmp = intersect(cone1_local_gs, cone2_neg_gs);
                max_neg(i) = mean(all_sr{sta_trial}(tmp)); % mean firing rate at these instances
                
                % find common instances when cone2 was positive
                tmp = intersect(cone1_local_gs, cone2_pos_gs);
                max_pos(i) = mean(all_sr{sta_trial}(tmp)); % mean firing rate at these instances
                
                % find common instances when cone2 was near 0 by response
                tmp = intersect(cone1_local_gs, cone2_aver_gs);
                max_aver(i) = mean(all_sr{sta_trial}(tmp)); % mean firing rate at these instances
                
            end
            cone2_effect_pos(:,sta_trial, cone2) = max_pos';
            cone2_effect_neg(:,sta_trial, cone2) = max_neg';
            cone2_effect_aver(:,sta_trial, cone2) = max_aver';
        end
    end
    
    my_cones = 1:size(inputs,1);
    my_cones(cone1) = [];
    put_c = intersect(my_cones, far_cones);
    real_cones = intersect(my_cones, real_center_cones);
    sur_cone = intersect(my_cones, real_sur_cones);
    
    subplot(1,3,1)
    hold on
    plot(squeeze(mean(cone2_effect_neg(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_neg(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_neg(:,:,sur_cone),2)),'rx-','linewidth',2)
    plot(mean(squeeze(mean(cone2_effect_neg(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 max negative')
    axis tight
    xlims = get(gca,'XLim');
    ylims = get(gca, 'YLim');
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(1,3,2)
    hold on
    plot(squeeze(mean(cone2_effect_aver(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_aver(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_aver(:,:,sur_cone),2)),'rx-','linewidth',2)
    plot(mean(squeeze(mean(cone2_effect_aver(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 near 0 resp')
    axis tight
    ylims = [ylims get(gca, 'YLim')];
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(1,3,3)
    hold on
    plot(squeeze(mean(cone2_effect_pos(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_pos(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_pos(:,:,sur_cone),2)),'rx-','linewidth',2)
    plot(mean(squeeze(mean(cone2_effect_pos(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 max positive')
    axis tight
    ylims = [ylims get(gca, 'YLim')];
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    for i = 1:3
        subplot(1,3,i)
        axis([xlims(1) xlims(2) min(ylims(:)), max(ylims(:))])
    end
end
    


% MODEL WITH ONE SUBUNIT

% common nonlinearity
nbins = 20;
my_nl=zeros(nbins,1);

my_input = sum(mfiltered_inputs(real_center_cones,:));
events_per_bin=floor(length(my_input)/nbins);
tmp=sort(my_input);
my_bins=[min(my_input) tmp(events_per_bin:events_per_bin:events_per_bin*nbins)];
my_bins(end)=max(my_input)+0.01;

for k=1:length(my_bins)-1
    tmp=find(my_input>=my_bins(k) & my_input<my_bins(k+1)); % find instances when GS had certain value
    my_nl(k)=mean(msr(tmp)); % find mean firing rate after this certain value
end

figure
plot(my_bins(1:end-1), my_nl)



% pass through nonlinearity
my_f_input = sum(ffiltered_inputs(real_center_cones,:));
predicted_rate = zeros(size(ffiltered_inputs,2),1);
for k=1:length(my_bins)-1
    tmp=find(my_f_input>=my_bins(k) & my_f_input<my_bins(k+1)); % find instances when GS had certain value
    predicted_rate(tmp) = my_nl(k);
end

sse = sum((predicted_rate-fsr).^2);
sst = sum((fsr-mean(fsr)).^2);
r2 = 1 - sse/sst;



% Estimate zero response bin in the model
for cone1 = 26:36
    
    figure
    set(gcf, 'Name',['Cone1 - ', int2str(cone1)])
    
    cone2_effect_pos = zeros(10, size(finputs,1));% bins, n of cones
    cone2_effect_neg = cone2_effect_pos;
    cone2_effect_aver = cone2_effect_pos;
    
    % bin cone1 input
    cone1_input = ffiltered_inputs(cone1,:); 
    n_gs=floor(length(cone1_input)/10);
    tmp=sort(cone1_input);
    my_bins=[min(cone1_input) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(cone1_input)+1;
        
    for cone2 = 1:size(finputs,1)
        if cone1 ~= cone2
            
            
            % find 10% of most negative, positive, and near 0 input of cone2
            cone2_input = ffiltered_inputs(cone2,:);
            tmp = sort(cone2_input);
            cone2_neg_gs = find(cone2_input<tmp(ceil(length(tmp)/10)));
            cone2_pos_gs = find(cone2_input>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone2_aver_gs = find( cone2_input>tmp(aver+2*bord) & cone2_input<tmp(aver+4*bord));
            
            
            
            clear max_neg max_pos max_aver
            for i=1:length(my_bins)-1
                
                % find instances of cone1 input in certain range
                cone1_local_gs = find( cone1_input>=my_bins(i) & cone1_input<my_bins(i+1));
                
                % find common instances when cone2 was negative
                tmp = intersect(cone1_local_gs, cone2_neg_gs);
                max_neg(i) = mean(fsr(tmp)); % mean firing rate at these instances
                
                % find common instances when cone2 was positive
                tmp = intersect(cone1_local_gs, cone2_pos_gs);
                max_pos(i) = mean(fsr(tmp)); % mean firing rate at these instances
                
                % find common instances when cone2 was near 0 by response
                tmp = intersect(cone1_local_gs, cone2_aver_gs);
                max_aver(i) = mean(fsr(tmp)); % mean firing rate at these instances
                
            end
            cone2_effect_pos(:,cone2) = max_pos';
            cone2_effect_neg(:,cone2) = max_neg';
            cone2_effect_aver(:,cone2) = max_aver';
        end
    end
    
    my_cones = 1:size(inputs,1);
    my_cones(cone1) = [];
    put_c = intersect(my_cones, far_cones);
    real_cones = intersect(my_cones, real_center_cones);
    sur_cone = intersect(my_cones, real_sur_cones);
    
    subplot(1,3,1)
    hold on
    plot(cone2_effect_neg(:,put_c))
    plot(cone2_effect_neg(:,real_cones),'linewidth',2)
    plot(cone2_effect_neg(:,sur_cone),'rx-','linewidth',2)
    plot(mean(cone2_effect_neg(:,put_c),2),'kx-','linewidth',4)
    title('cone2 max negative')
    axis tight
    xlims = get(gca,'XLim');
    ylims = get(gca, 'YLim');
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(1,3,2)
    hold on
    plot(cone2_effect_aver(:,put_c))
    plot(cone2_effect_aver(:,real_cones),'linewidth',2)
    plot(cone2_effect_aver(:,sur_cone),'rx-','linewidth',2)
    plot(mean(cone2_effect_aver(:,put_c),2),'kx-','linewidth',4)
    title('cone2 near 0 resp')
    axis tight
    ylims = [ylims get(gca, 'YLim')];
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(1,3,3)
    hold on
    plot(cone2_effect_pos(:,put_c))
    plot(cone2_effect_pos(:,real_cones),'linewidth',2)
    plot(cone2_effect_pos(:,sur_cone),'rx-','linewidth',2)
    plot(mean(cone2_effect_pos(:,put_c),2),'kx-','linewidth',4)
    title('cone2 max positive')
    axis tight
    ylims = [ylims get(gca, 'YLim')];
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    for i = 1:3
        subplot(1,3,i)
        axis([xlims(1) xlims(2) min(ylims(:)), max(ylims(:))])
    end


end








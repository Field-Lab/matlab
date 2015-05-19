function [stas, mean_bins, mean_nl, spike_rate]=partial_sta(inputs, spikes, fraction, sta_params, sel_cones, cones_to_plot)

spikes(spikes>size(inputs,2)) = [];

spikes_tmp = spikes;
spike_rate=zeros(size(inputs,2),1);
while ~isempty(spikes_tmp)
    [~, ia, ~] = unique(spikes_tmp);
    spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
    spikes_tmp(ia)=[];
end
clear spikes_tmp

bin_ng = 8;
n_sta = round(1/(1 - fraction));
bin = floor(size(inputs,2)/n_sta);
stas = zeros(size(inputs,1),sta_params.length, n_sta);
all_nl = zeros(size(inputs,1),bin_ng, n_sta);
all_bins = all_nl;
my_filtered_inputs=zeros(size(inputs,1),bin, n_sta);  
all_sr = cell(size(inputs,1),1);
figure
for i=1:n_sta
    spikes_tmp = spikes;
    spikes_tmp(spikes_tmp >= (bin*(i-1)+1) & spikes_tmp < bin*i) = [];
    
    % partial sta
    my_sta=zeros(size(inputs,1),sta_params.length);
    nspikes = numel(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:,sta_params.length-j+1)...
                + sum(inputs(:,spikes_tmp(ia) - sta_params.length + j + sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    stas(:,:,i)=my_sta/nspikes;
    
    %generator signal on the remainder part
    bin_points = (bin*(i-1)+1):bin*i;
    spike_rate_tmp = spike_rate((bin*(i-1)+sta_params.length):bin*i);
    all_sr{i} = spike_rate_tmp;
     
    for j=1:size(inputs,1)
        my_sta = stas(j,:,i);
        my_filtered_inputs(j,1:length(bin_points)-sta_params.length+1,i)=conv(inputs(j,bin_points),my_sta,'valid');
        tt = my_filtered_inputs(j,1:length(bin_points)-sta_params.length+1,i);
        n_gs=floor(length(tt)/bin_ng);
        tmp=sort(tt);
        my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*bin_ng)];
        my_bins(end)=max(tt)+1;
        my_nl=zeros(size(my_bins,2)-1,1);
        for k=1:length(my_bins)-1
            tmp=find(tt>=my_bins(k) & tt<my_bins(k+1)); % find instances when GS had certain value
            my_nl(k)=mean(spike_rate_tmp(tmp)); % find mean firing rate after this certain value
        end
        all_nl(j,:,i) = my_nl;
        all_bins(j,:,i) = my_bins(1:end-1);
        subplot(6,7,j)
        hold on
        plot(my_bins(1:end-1),my_nl)
    end
    
end

mean_bins = zeros(size(inputs,1),bin_ng);
mean_nl = mean_bins;
for i=1:size(inputs,1)
    subplot(6,7,i)
    mean_bins(i,:) = mean(all_bins(i,:,:),3);
    mean_nl(i,:) = mean(all_nl(i,:,:),3);
    plot(mean_bins(i,:), mean_nl(i,:), 'k','linewidth', 2)
    axis tight
end


far_cones = 1:length(sel_cones.far);
real_center_cones = (1:length(sel_cones.real)) + length(sel_cones.far);
real_sur_cones = (1:length(sel_cones.sur)) + length(sel_cones.far) + length(sel_cones.real);
all_cones = [sel_cones.far; sel_cones.real; sel_cones.sur];
bin_ng = 8;
frac = 20; % in %

% Estimate zero response bin
prc = 100/frac;
fit_acc = cell(length(cones_to_plot),10, 6);
fit_params = fit_acc;
for cone1 = 11:17%cones_to_plot
    
    figure(cone1)
    set(gcf, 'Name',['Cone1 - ', int2str(all_cones(cone1))])
    
    max_input = length(bin_points)-sta_params.length-5;
    cone2_effect_pos = zeros(bin_ng, n_sta, size(inputs,1));% bins, sta_trial, n of cones
    cone2_effect_neg = cone2_effect_pos;
    cone2_effect_aver = cone2_effect_pos;
    
    for sta_trial = 1:n_sta
        % bin cone1 input
        cone1_input = my_filtered_inputs(cone1,1:max_input,sta_trial);  % this is independent of sta
        n_gs=floor(length(cone1_input)/bin_ng);
        tmp=sort(cone1_input);
        my_bins=[min(cone1_input) tmp(n_gs:n_gs:n_gs*bin_ng)];
        my_bins(end)=max(cone1_input)+1;
        
        for cone2 = 1:size(inputs,1)
            if cone1 ~= cone2
                
                
                % find 10% of most negative, positive, and near 0 input of cone2
                cone2_input = my_filtered_inputs(cone2,1:max_input,sta_trial);
                tmp = sort(cone2_input);
                cone2_neg_gs = find(cone2_input<tmp(ceil(length(tmp)/prc)));
                cone2_pos_gs = find(cone2_input>tmp(ceil(length(tmp)/prc*(prc-1))));
                aver = ceil(length(tmp)/2);
                bord = ceil(length(tmp)/(prc*2));
                cone2_aver_gs = find( cone2_input>tmp(aver+0.8*bord) & cone2_input<tmp(aver+2.8*bord));
                
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
    end
    
    
    %     a = squeeze(mean(cone2_effect_neg(:,:,put_c),2));
    %     y = a(:,1);
    %     x = 1:8;
    %     save('/Users/alexth/Desktop/tmp.mat', 'a','x')
    
    
    my_cones = 1:size(inputs,1);
    my_cones(cone1) = [];
    put_c = intersect(my_cones, far_cones);
    real_cones = intersect(my_cones, real_center_cones);
    sur_cone = intersect(my_cones, real_sur_cones);
    
    
    
    subplot(2,3,1)
    hold on
    plot(squeeze(mean(cone2_effect_neg(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_neg(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_neg(:,:,sur_cone),2)),'rx-','linewidth',2)
    plot(mean(squeeze(mean(cone2_effect_neg(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 max negative')
    axis tight
    xlims = get(gca,'XLim');
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(2,3,4)
    hold on
    tmp = squeeze(mean(cone2_effect_neg(:,:,put_c),2));
    fit_acc{cone1,1} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(put_c)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[0.41,8.2,1.89]);
        fit_params{cone1,mc,1} =fit_res;
        fit_acc{cone1,1}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,1})
    tmp = squeeze(mean(cone2_effect_neg(:,:,real_cones),2));
    fit_acc{cone1,2} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(real_cones)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[0.41,8.2,1.89]);
        fit_params{cone1,mc,2} =fit_res;
        fit_acc{cone1,2}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,2},'linewidth',2)
    title('cone2 max negative, fit')
    axis tight
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    
    subplot(2,3,2)
    hold on
    plot(squeeze(mean(cone2_effect_aver(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_aver(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_aver(:,:,sur_cone),2)),'rx-','linewidth',2)
    legend(int2str(all_cones([put_c real_cones sur_cone])))
    plot(mean(squeeze(mean(cone2_effect_aver(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 near 0 resp')
    axis tight
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(2,3,5)
    hold on
    tmp = squeeze(mean(cone2_effect_aver(:,:,put_c),2));
    fit_acc{cone1,3} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(put_c)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[0.41,8.2,1.89]);
        fit_params{cone1,mc,3} =fit_res;
        fit_acc{cone1,3}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,3})
    tmp = squeeze(mean(cone2_effect_aver(:,:,real_cones),2));
    fit_acc{cone1,4} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(real_cones)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[0.41,8.2,1.89]);
        fit_params{cone1,mc,4} =fit_res;
        fit_acc{cone1,4}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,4},'linewidth',2)
    title('cone2 max negative, fit')
    axis tight
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    
    subplot(2,3,3)
    hold on
    plot(squeeze(mean(cone2_effect_pos(:,:,put_c),2)))
    plot(squeeze(mean(cone2_effect_pos(:,:,real_cones),2)),'linewidth',2)
    plot(squeeze(mean(cone2_effect_pos(:,:,sur_cone),2)),'rx-','linewidth',2)
    plot(mean(squeeze(mean(cone2_effect_pos(:,:,put_c),2)),2),'kx-','linewidth',4)
    title('cone2 max positive')
    axis tight
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    subplot(2,3,6)
    hold on
    tmp = squeeze(mean(cone2_effect_pos(:,:,put_c),2));
    fit_acc{cone1,5} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(put_c)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[0.41,8.2,1.89]);
        fit_params{cone1,mc,5} =fit_res;
        fit_acc{cone1,5}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,5})
    tmp = squeeze(mean(cone2_effect_pos(:,:,real_cones),2));
    fit_acc{cone1,6} = zeros(size(tmp));
    x = 1:bin_ng;
    for mc=1:length(real_cones)
        y = tmp(:,mc);
        fit_res=fit(x',y,'normcdf(x,mu,sigma)+b','Startpoint',[y(1),8.2,1.89], 'Lower', [0.2, 3,1], 'Upper', [2,12,3]);
        fit_params{cone1,mc,6} =fit_res;
        fit_acc{cone1,6}(:,mc) = normcdf(x,fit_res.mu,fit_res.sigma)+fit_res.b;
    end
    plot(fit_acc{cone1,6},'linewidth',2)
    title('cone2 max negative, fit')
    axis tight
    line([xlims(1) xlims(2)], [mean(spike_rate) mean(spike_rate)], 'color', 'k')
    
    mmin = 100;mmax = -100;
    for i = 1:6
        subplot(2,3,i)
        tt = get(gca, 'YLim');
        mmin = min(mmin, tt(1));
        mmax = max(mmax,tt(2));
    end
    for i = 1:6
        subplot(2,3,i)
        axis([xlims(1) xlims(2) mmin mmax])
    end
    
    drawnow
    
end

for cone1=11:17
    figure
    set(gcf,'Name',int2str(cone1))
    subplot(3,1,1)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,6})
            a = [a fit_params{cone1,i,6}.b];
        end
    end
    plot(a)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,5})
            a = [a fit_params{cone1,i,5}.b];
        end
    end
    hold on
    plot(a)
    axis([1 6 0 1])
    title('b')
    
    
    subplot(3,1,2)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,6})
            a = [a fit_params{cone1,i,6}.mu];
        end
    end
    plot(a)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,5})
            a = [a fit_params{cone1,i,5}.mu];
        end
    end
    hold on
    plot(a)
    axis([1 6 5 9])
    title('mu')
    
    subplot(3,1,3)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,6})
            a = [a fit_params{cone1,i,6}.sigma];
        end
    end
    plot(a)
    a = [];
    for i=1:6
        if ~isempty(fit_params{cone1,i,5})
            a = [a fit_params{cone1,i,5}.sigma];
        end
    end
    hold on
    plot(a)
    axis([1 6 1 5])
    title('sigma')
end



function contrast_response_ej_talk(inputs, spikes, nbins_cone1, nbins_cone2, center_cones, map, cones, comb, datarunID)

% contrast_bin_size in %: conditioning on the second cone
% inputs are filtered!

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
clear spikes_tmp

size_cone1_bin = floor(size(inputs,2)/nbins_cone1);
size_cone2_bin = floor(size(inputs,2)/nbins_cone2);

figure
set(gcf, 'position', [-1008         314         521         762])

colors = [1 0 0; 0 0.8 0; 0.5 0.5 0.5; 0 0 1; 0.8 0.8 0; 0 0.8 0.8; 0.8 0 0.8; 0.8 0.4 0; 0.3 0.8 0.6; 1 0 0.6; 0.6 0.65 0];

cone1 = 5;
subunit_cone = 4;
non_cone = 1; % use 4 for other cones

ax=[-0.28 0.23 0.1 0.95]; % cone 4 cell 168
% ax=[-0.25 0.2 0.1 0.8]; % cone 4 cell 168
% ax=[-0.2 0.15 0.08 0.37]; % cone 3 cell 271
% ax=[-0.2 0.15 0.1 0.3]; % cone 4 cell 271

tmp = sort(inputs(cone1,:));
cone1_contr = tmp(1:size_cone1_bin:end);
for j = 1:nbins_cone1
    cone1_valid = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
    cone1_uncond_rate(j) = mean(spike_rate(cone1_valid));
end

subplot(3,2,1)
hold on
for cone2 = 1:size(inputs,1)
    if cone1 ~= cone2
        cond_rate = zeros(nbins_cone2-1, nbins_cone1-1);
        linewid = 2;
        cols = colors(cone2,:);
        
        % second cone: bin the inputs
        tmp = sort(inputs(cone2,:));
        cone2_contr = tmp(1:size_cone2_bin:end);
        
        % go through contrast bins on cone2 and calculate the response
        % of cone 1
        
        for i=1
            % find instances when cone2 had certain contrast
            cone2_valid = find(inputs(cone2,:)>cone2_contr(i) & inputs(cone2,:)<=cone2_contr(i+1));
            
            % go through contrast range on cone1
            for j = 1:nbins_cone1
                cone1_valid = find(inputs(cone1,:)>cone1_contr(j) & inputs(cone1,:)<=cone1_contr(j+1));
                cone12_valid = intersect(cone1_valid, cone2_valid);
                cond_rate(i,j) = mean(spike_rate(cone12_valid));
            end
            %                  plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cond_rate(i,:)-mean(cond_rate(i,:)), 'linewidth', linewid, 'color', cols);
            plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cond_rate(i,:), 'linewidth', linewid, 'color', cols);
            
        end
        
        [params{cone2} params1{cone2} resp{cone2} resp1{cone2} x1{cone2}, resp2{cone2}, x2{cone2}] = fit_cone_interactions(cone1_contr(1:end-1)+diff(cone1_contr)/2, cone1_uncond_rate, cond_rate(1,:));
        
        %              hold on
        %              plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate)
        
    end
end
plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate,  '-*', 'color',[0.05 0.05 0.05], 'linewidth',2)
axis(ax)
set(gca, 'box','on')
set(gca, 'xtick', -0.1:0.1:0.1, 'ytick', 0:0.2:0.4)
title(['raw, ref cone ', int2str(cone1)])
% set(gca,'dataaspectratio', [1 1 1])


subplot(3,2,2)

plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate,  '-*', 'color',[0.05 0.05 0.05], 'linewidth',2)
hold on
plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,resp{non_cone},'color',colors(non_cone,:), 'linewidth', 2)
axis(ax)
set(gca, 'xtick', -0.1:0.1:0.1, 'ytick', 0:0.2:0.4)
title(['cond cone ',int2str(non_cone), ' y shifted'])
% set(gca,'dataaspectratio', [1 1 1])


subplot(3,2,3)

plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate,  '-*', 'color',[0.05 0.05 0.05], 'linewidth',2)
hold on
for cone2 = 1:size(inputs,1)
    if cone1 ~= cone2
        cols = colors(cone2,:);
        plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,resp{cone2},'color',cols, 'linewidth', 2)
    end
end
axis(ax)
set(gca, 'box','on')
set(gca, 'xtick', -0.1:0.1:0.1, 'ytick', 0:0.2:0.4)
title('all y shifted')
% set(gca,'dataaspectratio', [1 1 1])


subplot(3,2,4)

plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate, '-*', 'color',[0.05 0.05 0.05],'linewidth',2)
hold on
for i=1:length(center_cones)
    if cone1 ~= i && subunit_cone~=i
        cols = colors(i,:);
        plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,resp{i},'color',cols, 'linewidth', 2)
    end
end
plot(x1{subunit_cone},resp1{subunit_cone},'color',colors(subunit_cone,:), 'linewidth', 2)
plot(x2{subunit_cone},resp2{subunit_cone},'-x','color',colors(subunit_cone,:), 'linewidth', 2)
axis(ax)
set(gca, 'box','on')
set(gca, 'xtick', -0.1:0.1:0.1, 'ytick', 0:0.2:0.4)
title(['cond cone ',int2str(subunit_cone) ,' x shifted, others y shifted'])
% set(gca,'dataaspectratio', [1 1 1])


% subplot(3,3,5)
% plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate, '-*', 'color',[0.05 0.05 0.05],'linewidth',2)
% hold on
% plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,resp{subunit_cone},'color',colors(subunit_cone,:), 'linewidth', 2)
% plot(x1{subunit_cone},resp1{subunit_cone},'color',colors(subunit_cone,:), 'linewidth', 2)
%
% axis(ax)

% subplot(3,3,6)
% plot(cone1_contr(1:end-1)+diff(cone1_contr)/2,cone1_uncond_rate, '-*', 'color',[0.05 0.05 0.05],'linewidth',2)
% hold on
% plot(x2{subunit_cone},resp2{subunit_cone},'-o','color',colors(subunit_cone,:), 'linewidth', 2)
% plot(x1{subunit_cone},resp1{subunit_cone},'color',colors(subunit_cone,:), 'linewidth', 2)
% axis(ax)


colored_cones = zeros([size(comb)]);
for cone2 = 1:length(center_cones)
    [a, b] = find(map==center_cones(cone2));
    if cone1 ~= cone2
        for j = 1:length(a)
            colored_cones(a(j),b(j), :) = colors(cone2,:);
        end
    else % ref cone
        for j = 1:length(a)
            colored_cones(a(j),b(j),:) = 1;
        end
    end
end
subplot(3,2,5)
imagesc(colored_cones);
axis([290 360 460 525]) % cell 168
% axis([280 360 460 535]) % cell 271
set(gca,'dataaspectratio', [1 1 1])
% axis off

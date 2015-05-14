%% load stuff
datarun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data011-from-d08_11/data011-from-d08_11');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');


%% find stable cells

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

cell_type = 4;
bords = 30;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(600,600,3);
tmp(:,:,2) = tmp_map;

for i=49%find(ww./pixn<5)
    visionID = datarun.cell_types{cell_type}.cell_ids(i);
    datarunID = find(datarun.cell_ids == visionID);
    sta = imresize(double(datarun.stas.stas{datarunID}(:,:,1,4)),600/datarun.stimulus.field_height,'nearest');
    sta1 = imresize(double(datarun1.stas.stas{datarunID}(:,:,1,4)),600/datarun1.stimulus.field_height,'nearest');
    
    if datarun.stas.polarities{datarunID} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([600,600],ind);
    
    sta = sta/tt;
    sta1 = sta1/tt1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1; 
    figure(1)
    subplot(3,3,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])    
    title(['cell ID ', int2str(datarunID), ', vision ID ', int2str(visionID), ', i=', int2str(i)])
    
%     tt = unique(vormap(find(sta>0.4)));
%     tt(tt==0)=[];
%     clear sta_inds sta1_inds
%     for j = 1:length(tt)
%         a = find(vormap == tt(j));
%         tmp1 = sta(a);
%         tmp1 = tmp1(tmp1>0.4);
%         sta_inds(j)=sum(tmp1);
%         tmp1 = sta1(a);
%         tmp1 = tmp1(tmp1>0.4);
%         sta1_inds(j)=sum(tmp1);
%         [~, ic] = sort(sta_inds);
%     end
%     figure(2)
%     subplot(3,4,cnt)
%     plot(sta_inds(ic), sta1_inds(ic), '*')    
%     title(['cell ID ', int2str(cellInd), ', vision ID ', int2str(visionID)])    
    
    cnt = cnt+1;
end

% good_cells = [27 31 37 43]; % old classification
good_cells = [49];

%% preliminary, biased sta

datarunID = 168;

visionID = datarun.cell_ids(datarunID);

sta_params.length = 15;
sta_params.offset = 0;

my_sta=zeros(size(inputs,1),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
nspikes = numel(spikes);

spike_rate_all=zeros(size(inputs,2),1);
    
while ~isempty(spikes)
    [~, ia, ~] = unique(spikes);
    spike_rate_all(spikes(ia))=spike_rate_all(spikes(ia))+1;
    for j=1:sta_params.length
        my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
            sum(inputs(:,spikes(ia)-sta_params.length+j+sta_params.offset),2);
    end
    spikes(ia)=[];
end
my_sta=my_sta/nspikes;

% get voronoi map sta
tmp_map = vormap;
tt=0:max(tmp_map(:));
vorsta=zeros(600,600,sta_params.length);
coneX = zeros(max(tmp_map(:)),1);
coneY = coneX;
for i=1:937
    [a, b] = find(tmp_map==tt(i+1));
    % find center of this cone
    coneX(i) = mean(a);
    coneY(i) = mean(b);
    if ~isempty(a)
        for j = 1:length(a)
            vorsta(a(j),b(j),:) = my_sta(i,:);
        end
    end
end

figure
colormap gray
imagesc(vorsta(:,:,3))
hold on
for i=1:937
    text(coneY(i),coneX(i), int2str(i), 'color', 'r')
end

% plot voronoi sta frame by frame
figure
for i=1:9
    subplot(3,3,10-i)
    colormap gray
    imagesc(vorsta(:,:,i))
    title(['frame ',int2str(i)])
end

%plot sta to find the threshold
figure
plot(my_sta')
% cones with real center input
sel_cones.real = find(my_sta(:,3)<-0.065);
% cones with surround input
sel_cones.sur = [];%find(my_sta(:,3)>0.017);
% cones with no iput (furthest away)
realX = mean(coneX(sel_cones.real));
realY = mean(coneY(sel_cones.real));
tmp = pdist2([realX, realY], [coneX coneY]);
[~, ic] = sort(tmp);
ic(isnan(tmp(ic))) = [];
sel_cones.far = ic(end-9:end)';

cones = [sel_cones.far; sel_cones.real; sel_cones.sur];
vorsta_select_cones=zeros(600,600);
for i=cones'
    [a, b] = find(tmp_map==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            vorsta_select_cones(a(j),b(j)) = my_sta(i,3);
        end
    end
end

figure
subplot(1,2,1)
plot(my_sta(cones,:)')
title([int2str(length(cones)), ' cones'])

subplot(1,2,2)
figure(30)
colormap gray
imagesc(vorsta_select_cones)
hold on
for i=cones'
    text(coneY(i),coneX(i), int2str(i), 'color', 'r')
end

clear realX realY duration a b tmp tmp_map tt

%% unbiased sta/GS

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
fraction = 0.9;
cones_to_plot = 1:length(cones);
frames2take = floor(size(inputs,2)/2);

[stas, mean_bins, mean_nl, spike_rate]=partial_sta(inputs(cones,1:frames2take), spikes, fraction, sta_params, sel_cones, cones_to_plot);

mean_sta = mean(stas,3);

figure
plot(mean_sta')

load('/Users/alexth/Desktop/tmp.mat')
y = a(:,1);
% a = squeeze(mean(cone2_effect_neg(:,:,put_c),2));
% y = a(:,1);
% x = 1:8;
% save('/Users/alexth/Desktop/tmp.mat', 'a','x')

%% 

% filtered inputs
my_filtered_inputs=zeros(length(cones),size(inputs,2)-sta_params.length+1);  
for i=1:length(cones)
    my_filtered_inputs(i,:)=conv(inputs(cones(i),:),mean_sta(i,:),'valid');
end
spike_rate_tmp = spike_rate(sta_params.length:end);

% get response
resp_length = 20;
pre_spike = 10;
nbins = 15;

for cone1=1:length(cones)
    
    figure
    set(gcf,'Name',['Cone ', int2str(cone1)]);
    
    % bin input
    tt = my_filtered_inputs(cone1,:);
    n_gs=floor(length(tt)/nbins);
    tmp=sort(tt);
    my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*nbins)];
    my_bins(end)=max(tt)+1;
    
    mmin = 100; mmax = -100;
    for i=1:length(my_bins)-1
        tmp = find(tt>=my_bins(i) & tt<my_bins(i+1));
        subplot(3,5,i)
        hold on
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if ~isempty(tmp)
            resp=zeros(resp_length,1);
            for j=1:resp_length
                resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            end
            plot(resp, 'linewidth',2)
            title([num2str(my_bins(i)), ' to ', num2str(my_bins(i+1)), ', n = ', int2str(length(tmp))])
            mmin = min([mmin; resp]);
            mmax = max([mmax; resp]);
            
        end
    end
    
    for i = 1:length(my_bins)-1
        subplot(3,5,i)
        line([pre_spike+1, pre_spike+1], [ mmin mmax], 'color','k')
        line([1 resp_length], [mean(spike_rate_tmp) mean(spike_rate_tmp)])
        axis([1 resp_length mmin mmax])
    end


end



% TRIAL
resp_length = 30;
pre_spike = 15;
same_n = false;
null_second_cone = true;
for cone1 = 1:8
    for cone2 = 1:8
        if cone1 ~= cone2
            % bin cone1 input
            tt = my_filtered_inputs(cone1,:);            
            n_gs=floor(length(tt)/9);
            tmp=sort(tt);
            my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*9)];
            my_bins(end)=max(tt)+1;
            
            % find 10% of most negative, positive, and near 0 input of cone2
            tt1 = my_filtered_inputs(cone2,:);
            tmp = sort(tt1);
            cone2_neg_gs = find(tt1<tmp(ceil(length(tmp)/10)));
            cone2_pos_gs = find(tt1>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone2_aver_gs = find( tt1>tmp(aver-bord) & tt1<tmp(aver+bord));
           
            
            figure
            set(gcf, 'Name',['cone ', int2str(cone1), ' vs cone ', int2str(cone2)])
            
            % cone1 varying input, cone2 always negative (5%)
            for i=1:length(my_bins)-1
                
                subplot(3,3,i)
                hold on
                
                % find cone1 binned input instances
                cone1_local_gs = find( tt>=my_bins(i) & tt<my_bins(i+1));
                
                % find when cone2 was negative                                
                tmp = intersect(cone1_local_gs, cone2_neg_gs);                                
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                plot(resp, 'r+-','linewidth',2)                

                % find when cone2 was positive                                
                tmp = intersect(cone1_local_gs, cone2_pos_gs);
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                plot(resp, 'bx-','linewidth',2)
               
                % find when cone2 was near 0                                
                tmp = intersect(cone1_local_gs, cone2_aver_gs);
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                plot(resp, 'ko-','linewidth',2)
                
                axis tight
                
                title([num2str(my_bins(i)), ' to ', num2str(my_bins(i+1))])
                
                if i==5
                    legend('neg', 'pos', '0')
                end
            end
            subplot(3,3,1)
            ylmin = get(gca,'YLim');
            subplot(3,3,9)
            ylmax = get(gca,'YLim');
            for i=1:length(my_bins)-1               
                subplot(3,3,i)
                axis([1 resp_length ylmin(1) ylmax(2)])
                line([1 resp_length], [0,0], 'color', [1 1 1]*0.3)
                line([pre_spike+1 pre_spike+1], [ylmin(1) ylmax(2)], 'color', [1 1 1]*0.3)
            end
        end
    end
end



% TRIAL SUMMARY
resp_length = 30;
pre_spike = 15;
for cone1 = 1:length(cones)

    figure
    set(gcf, 'Name',['cone ', int2str(cone1)])
      
    for cone2 = 1:length(cones)
        if cone1 ~= cone2
            % bin cone1 input
            tt = my_filtered_inputs(cone1,:);            
            n_gs=floor(length(tt)/9);
            tmp=sort(tt);
            my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*9)];
            my_bins(end)=max(tt)+1;
            
            % find 10% of most negative, positive, and near 0 input of cone2
            tt1 = my_filtered_inputs(cone2,:);
            tmp = sort(tt1);
            cone2_neg_gs = find(tt1<tmp(ceil(length(tmp)/10)));
            cone2_pos_gs = find(tt1>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone2_aver_gs = find( tt1>tmp(aver-bord) & tt1<tmp(aver+bord));
           
            
     
            clear max_neg max_pos max_aver
            for i=1:length(my_bins)-1
      
                % find cone1 binned input instances
                cone1_local_gs = find( tt>=my_bins(i) & tt<my_bins(i+1));
                
                % find when cone2 was negative                                
                tmp = intersect(cone1_local_gs, cone2_neg_gs);                                
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                max_neg(i) = resp(16);             

                % find when cone2 was positive                                
                tmp = intersect(cone1_local_gs, cone2_pos_gs);
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                max_pos(i) = resp(16);
               
                % find when cone2 was near 0                                
                tmp = intersect(cone1_local_gs, cone2_aver_gs);
                tmp(tmp<(pre_spike+1)) = [];
                tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
                n_inps = length(tmp);
                resp=zeros(resp_length,1);
                for j=1:resp_length
                    resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                end
                max_aver(i) = resp(16);
                
            end
            
            subplot(3,ceil(length(cones)/3),cone2)
            hold on
            plot(max_neg,'xr-')
%             plot(max_pos,'+b-')
            plot(max_aver,'ko-')
            line([1 length(my_bins)], [mean(spike_rate_tmp) mean(spike_rate_tmp)], 'color', 'k')
            title(['cone ', int2str(cone2)])
        end
    end
end




% Estimate zero response bin
resp_length = 10;
pre_spike = 5;
for cone1 = 27:30%29:length(cones)

    figure
    set(gcf, 'Name',['cone ', int2str(cone1)])
      
    all_acc=0; all_acc_neg = 0; all_acc_aver = 0;
    cnt = 0;
    for cone2 = 1:length(cones)
        if cone1 ~= cone2
            % bin cone1 input
            cone1_input = my_filtered_inputs(cone1,:);            
            n_gs=floor(length(cone1_input)/10);
            tmp=sort(cone1_input);
            my_bins=[min(cone1_input) tmp(n_gs:n_gs:n_gs*10)];
            my_bins(end)=max(cone1_input)+1;
            
            % find 10% of most negative, positive, and near 0 input of cone2
            cone2_input = my_filtered_inputs(cone2,:);
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
                max_neg(i) = mean(spike_rate_tmp(tmp)); % mean firing rate at these instances

                % find common instances when cone2 was positive                                
                tmp = intersect(cone1_local_gs, cone2_pos_gs);
                max_pos(i) = mean(spike_rate_tmp(tmp)); % mean firing rate at these instances
               
                % find common instances when cone2 was near 0 by response                              
                tmp = intersect(cone1_local_gs, cone2_aver_gs);
                max_aver(i) = mean(spike_rate_tmp(tmp)); % mean firing rate at these instances
               
            end
            
            subplot(1,3,1)
            hold on
            if cone2<length(put_cones) % far away cones
                all_acc_neg = all_acc_neg+max_neg;
                cnt = cnt+1;
                plot(max_neg);
            elseif cone2==length(put_cones) % mean of far away cones
                all_acc_neg = all_acc_neg+max_neg;
                cnt = cnt+1;
                plot(max_neg);
                plot(all_acc_neg/cnt,'o-r', 'linewidth',4)
            elseif cone2<length(cones) % real input cones
                plot(max_neg,'x-', 'linewidth',2)
            else % surround cone
                plot(max_neg,'x-k', 'linewidth',4)
            end
            title('cone2 max negative')
            line([1 length(my_bins)], [mean(spike_rate_tmp) mean(spike_rate_tmp)], 'color', 'k')

            subplot(1,3,2)
            hold on
            if cone2<length(put_cones) % far away cones
                all_acc_aver = all_acc_aver+max_aver;
                plot(max_aver);
            elseif cone2==length(put_cones) % mean of far away cones
                all_acc_aver = all_acc_aver+max_aver;
                plot(max_aver);
                plot(all_acc_aver/cnt,'o-r', 'linewidth',4)
            elseif cone2<length(cones) % real input cones
                plot(max_aver,'x-', 'linewidth',2)
            else % surround cone
                plot(max_aver,'x-k', 'linewidth',4)
            end
            title('cone2 near 0 by response')
            line([1 length(my_bins)], [mean(spike_rate_tmp) mean(spike_rate_tmp)], 'color', 'k')
            
            subplot(1,3,3)
            hold on
            if cone2<length(put_cones) % far away cones
                all_acc = all_acc+max_pos; 
                plot(max_pos);
            elseif cone2==length(put_cones) % mean of far away cones
                all_acc = all_acc+max_pos; 
                plot(max_pos);
                plot(all_acc/cnt,'o-r', 'linewidth',4)
            elseif cone2<length(cones) % real input cones
                plot(max_pos,'x-', 'linewidth',2)
            else % surround cone
                plot(max_pos,'x-k', 'linewidth',4)
            end
            title('cone2 max positive')
            line([1 length(my_bins)], [mean(spike_rate_tmp) mean(spike_rate_tmp)], 'color', 'k')
        end
    end
    
    ylims = [];
    for i = 1:3
        subplot(1,3,i)
        axis tight
        ylims = [ylims; get(gca,'YLim')];
    end
    for i = 1:3
        subplot(1,3,i)
        axis([1 10 min(ylims(:)), max(ylims(:))])
    end
    
end




% joint or linear sum
norm_sr = spike_rate_tmp - mean(spike_rate_tmp);
figure
range_cones = 25:40;
cnt = 1;
for cone1 = range_cones

    subplot(4,4,cnt)
      
    oppos_simult=zeros(1,length(cones)); pos_zero_simult = oppos_simult; zero_neg_simult = oppos_simult;
    for cone2 = range_cones
        if cone1 ~= cone2
            % find 10% of most negative, positive, and near 0 input of cone1
            cone1_input = my_filtered_inputs(cone1,:);            
            tmp = sort(cone1_input);
            cone1_neg_gs = find(cone1_input<tmp(ceil(length(tmp)/10)));
            cone1_pos_gs = find(cone1_input>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone1_aver_gs = find( cone1_input>tmp(aver+2*bord) & cone1_input<tmp(aver+4*bord));
            
            % find 10% of most negative, positive, and near 0 input of cone2
            cone2_input = my_filtered_inputs(cone2,:);
            tmp = sort(cone2_input);
            cone2_neg_gs = find(cone2_input<tmp(ceil(length(tmp)/10)));
            cone2_pos_gs = find(cone2_input>tmp(ceil(length(tmp)/10*9)));
            aver = ceil(length(tmp)/2);
            bord = ceil(length(tmp)/20);
            cone2_aver_gs = find( cone2_input>tmp(aver+2*bord) & cone2_input<tmp(aver+4*bord));
            
            % cone1 max pos, cone2 max neg, co-occuring 
            tmp = intersect(cone1_pos_gs, cone2_neg_gs); 
            oppos_simult(cone2) = mean(norm_sr(tmp)); % mean firing rate
            
            % cone1 max pos, cone2 0, co-occuring
            tmp = intersect(cone1_pos_gs, cone2_aver_gs); 
            pos_zero_simult(cone2) = mean(norm_sr(tmp)); % mean firing rate

            % cone1 0, cone2 max neg, co-occuring
            tmp = intersect(cone1_aver_gs, cone2_neg_gs); 
            zero_neg_simult(cone2) = mean(norm_sr(tmp)); % mean firing rate            
        end
    end
    
    plot(range_cones, oppos_simult(range_cones), 'ro-', 'linewidth',3);
    hold on
%     plot(range_cones, pos_zero_simult(range_cones), 'cx-');
%     plot(range_cones, zero_neg_simult(range_cones), 'mx-');
    plot(range_cones, pos_zero_simult(range_cones) + zero_neg_simult(range_cones), 'bx-', 'linewidth',3);
    if cnt ==1
%         legend('+/-', '+/0', '0/-', '+/0 sum 0/+')
        legend('+/-', '+/0 sum 0/+')
    end
    ylims = get(gca, 'Ylim');
    line([range_cones(1) range_cones(end)], [0 0], 'color', 'k');
    title(['cone ',int2str(cone1)])
    axis([range_cones(1) range_cones(end) min([ylims(1)-0.05 -0.1]), max([ylims(2)+0.05 +0.1])])
    ylims = get(gca, 'Ylim');
    line([cone1 cone1], [ylims(1) ylims(2)], 'color', 'k');
    line([length(put_cones)+0.5 length(put_cones)+0.5], [ylims(1) ylims(2)], 'color', 'k');
    cnt = cnt+1;
end


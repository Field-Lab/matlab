%% load stuff

local_path = '/Users/alexandra/Desktop/datasets/';

datarun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11']);
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);

datarun1 = load_data([local_path, '2011-12-13-2/d08-11-norefit/data011-from-d08_11/data011-from-d08_11']);
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

vormap = load([local_path, '2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt']);
figure
colormap gray
imagesc(vormap)


vorrun = load_data([local_path, '2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11']);
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_sta(vorrun);
vorrun = load_neurons(vorrun);
tic
[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');
toc
% rgb_flag  = 'bw'; bin_flag = 'bin'; rgb = 0.48; seed = 11111; width = 1; height = 937; probability = 1; frames = 54351;
% tic
% my_movie = compute_raw_wn_movie(rgb_flag, bin_flag, rgb, seed, width, height, probability, frames);
% toc
% inputs = squeeze(my_movie(1,:,1,:))/255-0.5;
% clear rgb_flag bin_flag rgb seed width height probability frames my_movie


%% preliminary, biased sta
find(datarun.cell_ids==6577)

% off midgets
datarunID = 168; % 2 pairs? thresh -0.07
datarunID = 162; % 4 pairs (but beware of voronoi): 1/2 (714/718), 2/5 (718/721), 3/8 (719/759), 6/7 (727/729). Threshold -0.08
datarunID = 52; % ID 1321; bad movement
datarunID = 69; % ID 1741; bad voronoi
datarunID = 103; % ID 2476; bad voronoi
datarunID = 108; % ID 2611; bad voronoi
datarunID = 112; % ID 2657; thresh -0.1
datarunID = 125; % ID 2836; possibly 1 pair 2/4 (406/408)
datarunID = 126; % ID 2837; crappy STAs
datarunID = 128; % ID 2971; bad voronoi
datarunID = 148; % ID 3331; no pairs
datarunID = 207; % ID 4621; no pairs
datarunID = 286; % ID 6901; 1 pair 2/8 (481/543). thresh -0.05

datarunID = 280; % ID 6712

datarunID = 273; % ID 6577



% On midgets
datarunID = 36; % ID 857; 2 pairs: cones 2/3 (374/375) and 4/5 (424/426)
datarunID = 64; % big movement
datarunID = 68; % bad voronoi
datarunID = 113; % bad voronoi
datarunID = 120;% ID 2779 bad voronoi
datarunID = 123;% ID 2821; 3 pairs: 3/7 (538/615), 3/4 (538/539) maybe, 1/5 (479/540) - possibly 2 BC connecting to 3 cones. Off - 207, 280
datarunID = 244;% ID 5718; 2 good cones, 3rd moved. 2 good ones are in pair: 1/2 (693/743)
datarunID = 263;% ID 6392; possibly pairs, but too many cones, messed up somewhat
datarunID = 265;% ID 6436; 2 pairs: 2/5 (435/496), 2/4 (435/438) - possibly 2 BC connecting to 3 cones. Off:273
datarunID = 271;% ID 6541; 5 pairs: 6/9 (549/626), 6/7 (549/550), 5/8 (546/623), 5/9 (546/626), 3/4 (493/494). Off cell:168,273
datarunID = 298; % ID 7471; probably no pairs


visionID = datarun.cell_ids(datarunID);
raw_sta = squeeze(vorrun.stas.stas{datarunID});
[full_sta, cones] = expand_voronoi_sta(raw_sta, vormap);

% plot voronoi sta frame by frame
figure
for i=1:30
    subplot(5,6,i)
    colormap gray
    imagesc(full_sta(:,:,i))
    title(['frame ',int2str(i)])
end

figure
colormap gray
imagesc(full_sta(:,:,27))
hold on
for i=1:937
    text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
end


figure
plot(raw_sta')

center_cones = find(raw_sta(:,27)>0.05);
% center_cones = find(raw_sta(:,27)<-0.07);

x = mean(cones(center_cones,1));
y = mean(cones(center_cones,2));
tmp = pdist2([x, y], [cones(:,1) cones(:,2)]);
[~, ic] = sort(tmp);
ic(isnan(tmp(ic))) = [];
far_cones = ic(end-15:end)';

select_cones = [center_cones; far_cones];

figure
colormap gray
imagesc(full_sta(:,:,27))
hold on
for i=select_cones
    text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
end


sta = squeeze(datarun.stas.stas{datarunID});
sta1 = squeeze(datarun1.stas.stas{datarunID});
sta_snippet = -imresize(double(sta(:,:,4)), 2, 'nearest');
sta1_snippet = -imresize(double(sta1(:,:,4)), 2, 'nearest');
voronoi_regions = full_sta(:,:,27);
% voronoi_regions = voronoi_regions/(max(voronoi_regions(:))*1.5);
voronoi_regions = voronoi_regions/(min(voronoi_regions(:))*1.5);
comb = zeros(600,600,3);
comb(:,:,1) = sta_snippet/max(sta_snippet(:));
comb(:,:,3) = sta1_snippet/max(sta1_snippet(:));
comb(:,:,2) = voronoi_regions;
figure
imagesc(comb)
hold on
for i=select_cones
    text(cones(i,2),cones(i,1), int2str(i), 'color', 'r')
end



%% unbiased sta/GS
sta_params.length = 15;
sta_params.offset = 0;
fraction = 0.9;

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];

[unbiased_sta, gensig_bins, nonlinearity]=unbiased_STA(inputs(select_cones,1:50000), spikes, fraction, sta_params);

figure
ax = zeros(length(select_cones),1);
[r,c] = opt_subplots(length(select_cones));
for i=1:length(select_cones)
    tmp = gensig_bins(i,:);
    ax(i) = subplot(r,c,i);
    plot(tmp(1:end-1),nonlinearity(i,:))
    axis tight
end
set(ax, 'YLim',[min(nonlinearity(:)), max(nonlinearity(:))])
figure
plot(unbiased_sta')

filt_inputs = zeros(length(select_cones), size(inputs,2)-sta_params.length+1);
cnt = 1;
for current_cone=select_cones'
    filt_inputs(cnt,:)=conv(inputs(current_cone,:), unbiased_sta(cnt,:),'valid');
    cnt=cnt+1;
end
spikes_tmp = spikes;
spikes_tmp(spikes<sta_params.length) = [];
nbins_cone1 = 5;
nbins_cone2 = 6;
contrast_response_cone(filt_inputs, spikes_tmp-sta_params.length+1, nbins_cone1, nbins_cone2, center_cones, vormap, cones, comb);

%% combination of voronoi and cones

comb2 = comb;
comb1 = comb(:,:,1);
comb1(comb1<0.3) = 0;
comb2(:,:,1) = comb1*1.5;
comb1 = comb(:,:,3);
comb1(comb1<0.3) = 0;
comb2(:,:,3) = comb1*1.5;
comb1 = comb(:,:,2);
comb1(comb1<0.3) = 0;
comb2(:,:,2) = comb1;
figure
imagesc(comb2)
set(gca,'dataaspectratio',[1 1 1])
axis([290 355 470 515])



%% 
spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_params.length-sta_params.offset)=[];
fraction = 0.9;
% cones_to_plot = 1:length(cones);
% frames2take = floor(size(inputs,2)/2);
% [stas, mean_bins, mean_nl, spike_rate]=partial_sta(inputs(cones,1:frames2take), spikes, fraction, sta_params, sel_cones, cones_to_plot);

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


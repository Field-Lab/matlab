datarun = load_data('/Volumes/Analysis/2011-10-25-9/d06-10-norefit/data006-from-d06_10/data006-from-d06_10');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);


datarun1 = load_data('/Volumes/Analysis/2011-10-25-9/d06-10-norefit/data010-from-d06_10/data010-from-d06_10');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);

vorrun = load_data('/Volumes/Analysis/2011-10-25-9/d06-10-norefit/data008-from-d06_10/data008-from-d06_10');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

vormap = load('/Volumes/Data/2011-10-25-9/Visual/2011-10-25-9_f06_vorcones/map-0000.txt');
figure
imagesc(vormap)
max_ind = max(vormap(:));

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

offm = datarun.cell_types{4}.cell_ids;
offp = datarun.cell_types{2}.cell_ids;
sbc = datarun.cell_types{5}.cell_ids;
bords = 30;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(600,600,3);
tmp(:,:,2) = tmp_map;
for i=1:12
    cellInd = find(datarun.cell_ids == offm(i));
    sta = imresize(double(datarun.stas.stas{cellInd}(:,:,1,4)),3,'nearest');
    sta1 = imresize(double(datarun1.stas.stas{cellInd}(:,:,1,4)),3,'nearest');
    
    if datarun.stas.polarities{cellInd} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([600,600],ind);
    
    sta = sta/tt;
%     sta(sta>0.25) = 1;

    sta1 = sta1/tt1;
%     sta1(sta1>0.25) = 1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1;    
    subplot(3,4,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])
    title(['cell ', int2str(cellInd), ' vision ID ', int2str(offm(i))])
    cnt = cnt+1;
end

cells = [10 14 31 33 39 40 41 42 44 45 47 48 49 51 54 55 56 60 65 70 72 78];

cells = [14];

clear bords sta sta1 ind tt tt1 cellInd a b

[inputs, refresh, duration] = get_wn_movie_ath(vorrun, 'BW-1-3-0.50-11111-1264x1-60.35.xml');

offset = 0;
sta_length = 15;
my_sta=zeros(max_ind,sta_length, length(cells));
spike_rate=zeros(duration,length(cells)); 
cnt = 1;
for i = cells
    cellInd = find(datarun.cell_ids == offm(i));
    spikes=ceil((vorrun.spikes{cellInd}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
    spikes(spikes<sta_length-offset)=[];
    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        spike_rate(spikes(ia), cnt)=spike_rate(spikes(ia), cnt)+1;
        for j=1:sta_length
            my_sta(:,sta_length-j+1, cnt)=my_sta(:,sta_length-j+1, cnt)+sum(inputs(:,spikes(ia)-sta_length+j+offset),2);
        end
        spikes(ia)=[];
    end
    my_sta(:,:,cnt)=my_sta(:,:,cnt)/sum(spike_rate(:,cnt));
    cnt = cnt+1;
end

tmp_map = vormap;

my_cell = 1;
tt=0:max_ind;
voronoi_sta=zeros(600,600,sta_length);
for i=1:max_ind
    [a, b] = find(vormap==tt(i+1));
    if ~isempty(a)
        for j = 1:length(a)
            voronoi_sta(a(j),b(j),:) = my_sta(i,:, my_cell);
        end
    end
end

figure
for i=1:9
    subplot(3,3,10-i)
    colormap gray
    imagesc(voronoi_sta(:,:,i))
    title(['frame back ', int2str(i)])
%     axis([300 365 460 520])
end

figure
colormap gray
imagesc(voronoi_sta(:,:,2))

my_cell_sta = my_sta(:,:,my_cell);

cones=find(my_cell_sta(:,2)<-0.03);
figure
subplot(3,1,1)
plot(my_cell_sta')
subplot(3,1,2)
plot(my_cell_sta(cones,:)')
subplot(3,1,3)
plot(my_cell_sta(cones,2))
set(gca, 'xtick', 1:length(cones), 'xticklabel', {int2str(cones)})

cellInd = find(datarun.cell_ids == offm(cells(my_cell)));

spikes=ceil((vorrun.spikes{cellInd}-vorrun.triggers(1))*1000/(refresh)); % spikes in frames
spikes(spikes<sta_length-offset)=[];

% get convolved inputs

my_filtered_inputs=zeros(size(cones,1),size(inputs,2)-size(my_cell_sta,2)+1);
figure
for i=1:length(cones)
    my_filtered_inputs(i,:)=conv(inputs(cones(i),:),my_cell_sta(cones(i),:),'valid');
    subplot(5,6,i)
    hist(my_filtered_inputs(i,:),20)
    title(num2str(std(my_filtered_inputs(i,:))))
end
spike_rate_tmp =spike_rate(sta_length:end, my_cell);
spike_rate_tmp = spike_rate_tmp - mean(spike_rate_tmp);

cone1 = 4;

tt = my_filtered_inputs(cone1,:);
n_gs=floor(length(tt)/10);
tmp=sort(tt);
my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
my_bins(end)=max(tt)+1;

my_nl=zeros(size(my_bins,2)-1,1);
my_std=my_nl;
for i=1:length(my_bins)-1
    tmp=find(tt>=my_bins(i) & tt<my_bins(i+1)); % find instances when GS had certain value
    my_std(i)=std(spike_rate_tmp(tmp));
    my_nl(i)=mean(spike_rate_tmp(tmp)); % find mean firing rate after this certain value
end
figure
plot(my_bins(1:end-1),my_nl)
tmp = mean(spike_rate_tmp);
line([my_bins(1) my_bins(end-1)], [tmp, tmp], 'color','k')
axis tight

resp_length = 30;
pre_spike = 15;
figure
for i=1:length(my_bins)-1
    tmp = find(tt>=my_bins(i) & tt<my_bins(i+1));
    subplot(3,4,i)
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
        end
        plot(resp, 'linewidth',2)
        title([num2str(my_bins(i)), ' to ', num2str(my_bins(i+1)), ', n = ', int2str(length(tmp))])
%         line([pre_spike, pre_spike], [0.2 0.81], 'color','k')
%         axis([1 50 0.2 0.81])
    end
end

% same inputs
cone1 = 3;
cone2 = 7;
tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
n_gs=floor(length(tt)/10);
tmp=sort(tt);
my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
my_bins(end)=max(tt)+1;

resp_length = 30;
pre_spike = 15;
figure
for i=1:length(my_bins)-1
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)>=my_bins(i) & my_filtered_inputs(cone2,:)<my_bins(i+1));
    tmp = intersect(tmp, tmp1);
    subplot(3,4,i)
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
        end
        plot(resp, 'b','linewidth',2)
        title([num2str(my_bins(i)), ' to ', num2str(my_bins(i+1)), ', n = ', int2str(length(tmp))])
        line([1 resp_length], [mean(spike_rate_tmp) mean(spike_rate_tmp)], 'color','k')
        line([pre_spike+1, pre_spike+1], [0 1.3], 'color','k')
        axis([1 resp_length 0 1.3])
    end
end



% opposite inputs against linear sum
cone1 = 5;
cone2 = 3;
tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
n_gs=floor(length(tt)/10);
tmp=sort(tt);
my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
my_bins(end)=max(tt)+1;
% spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
resp_length = 30;
pre_spike = 15;
same_n = false;
figure
for i=1:length(my_bins)-1
    % 2 inputs together
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp = intersect(tmp, tmp1);
    subplot(3,4,i)
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
        end
        plot(resp, 'r+-','linewidth',2)
        title({[num2str(my_bins(i)), ' to ', num2str(my_bins(i+1)), ', n = ', int2str(length(tmp))],...
            [num2str(my_bins(end-i)), ' to ', num2str(my_bins(end-i+1))]})
        %         line([pre_spike, pre_spike], [0.2 0.55], 'color','k')
        %         axis([1 50 0.1 0.65])
    end
    hold on
    n_inps = length(tmp);
    
    % input 1 separately
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if same_n
        tmp = sort(tmp(randperm(length(tmp),n_inps)));
    end
    resp1=zeros(resp_length,1);
    for j=1:resp_length
        resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
    end
    
    % input 2 separately
    tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if same_n
        tmp = sort(tmp(randperm(length(tmp),n_inps)));
    end
    resp2=zeros(resp_length,1);
    for j=1:resp_length
        resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
    end
    
    % linear sum
    plot(resp1 + resp2, 'kx-','linewidth',2);
%      plot(resp1, 'kx-','linewidth',2);
    axis([1 resp_length -0.3 0.4])
    if i==5
        legend('joint','linear')
    end
end


ttt(kk)=resp1(16)+resp2(16);


for cone1 = 1:15

    figure
    set(gcf, 'name', ['ref cone ', int2str(cone1)])
    for cone2 = 1:15
        tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
        n_gs=floor(length(tt)/10);
        tmp=sort(tt);
        my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
        my_bins(end)=max(tt)+1;
        % spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
        resp_length = 30;
        pre_spike = 15;
        same_n = true;
        i=length(my_bins)-1;
        % 2 inputs together
        tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
        tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
        tmp = intersect(tmp, tmp1);
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if ~isempty(tmp)
            resp=zeros(resp_length,1);
            for j=1:resp_length
                resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                resp_std(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
            end
        end
        n_inps = length(tmp);
        
        for kk=1:20
            % input 1 separately
            tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
            tmp(tmp<(pre_spike+1)) = [];
            tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
            if same_n
                tmp = sort(tmp(randperm(length(tmp),n_inps)));
            end
            resp1=zeros(resp_length,1);
            for j=1:resp_length
                resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                resp_std1(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
            end
            
            % input 2 separately
            tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
            tmp(tmp<(pre_spike+1)) = [];
            tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
            if same_n
                tmp = sort(tmp(randperm(length(tmp),n_inps)));
            end
            resp2=zeros(resp_length,1);
            for j=1:resp_length
                resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
                resp_std2(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
            end
            
            ttt(kk)=resp1(16)+resp2(16);
            ttt_ste(kk)=mean(resp_std1(16)+resp_std2(16));
        end
        
        subplot(3,5,cone2)
        plot(ttt, 'linewidth',2)
        hold on
        plot(ttt + ttt_ste, '--b')
        plot(ttt - ttt_ste, '--b')
        plot(ones(20,1)*resp(16), 'linewidth',2)
        plot(ones(20,1)*(resp(16) + resp_std(16)),'--r')
        plot(ones(20,1)*(resp(16) - resp_std(16)),'--r')
        title(['cone ', int2str(cone2), ', n = ', int2str(n_inps)])
        drawnow
    end

end


% for pairs
cone_pairs = [5, 2; 4, 7; 3, 2; 1, 3; 1, 2];
cnt=1;
test_n = 10;
clear ttt ttt_ste
null_second_cone = false
same_n = false
figure
for ii = 1:length(cone_pairs)
    cone1 = cone_pairs(ii,1);
    cone2 = cone_pairs(ii,2);
    tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
    n_gs=floor(length(tt)/10);
    tmp=sort(tt);
    my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(tt)+1;
    % spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
    resp_length = 30;
    pre_spike = 15;
    i=length(my_bins)-1;
    % 2 inputs together
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp = intersect(tmp, tmp1);
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
    end
    n_inps = length(tmp);
        
        
    for kk=1:test_n
        % input 1 separately
        tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone2,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp1=zeros(resp_length,1);
        for j=1:resp_length
            resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std1(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        % input 2 separately
        tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone1,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone1,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp2=zeros(resp_length,1);
        for j=1:resp_length
            resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std2(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        ttt(kk)=resp1(16)+resp2(16);
        ttt_ste(kk)=mean(resp_std1(16)+resp_std2(16));
    end
    
    
    subplot(3,4,cnt)
    plot(ttt, 'linewidth',2)
    hold on
    plot(ttt + ttt_ste, '--b')
    plot(ttt - ttt_ste, '--b')
    plot(ones(test_n,1)*resp(16), 'linewidth',2)
    plot(ones(test_n,1)*(resp(16) + resp_std(16)),'--r')
    plot(ones(test_n,1)*(resp(16) - resp_std(16)),'--r')
    title(['cone ', int2str(cone1),' vs cone ', int2str(cone2), ', n = ', int2str(n_inps)])
    drawnow
    cnt=cnt+1;
    
    cone2 = cone_pairs(ii,1);
    cone1 = cone_pairs(ii,2);
    tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
    n_gs=floor(length(tt)/10);
    tmp=sort(tt);
    my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(tt)+1;
    % spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
    resp_length = 30;
    pre_spike = 15;
    i=length(my_bins)-1;
    % 2 inputs together
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp = intersect(tmp, tmp1);
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
    end
    n_inps = length(tmp);
        
    for kk=1:test_n
        % input 1 separately
        tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone2,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp1=zeros(resp_length,1);
        for j=1:resp_length
            resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std1(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        % input 2 separately
        tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone1,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone1,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp2=zeros(resp_length,1);
        for j=1:resp_length
            resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std2(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        ttt(kk)=resp1(16)+resp2(16);
        ttt_ste(kk)=mean(resp_std1(16)+resp_std2(16));
    end
    
    subplot(3,4,cnt)
    plot(ttt, 'linewidth',2)
    hold on
    plot(ttt + ttt_ste, '--b')
    plot(ttt - ttt_ste, '--b')
    plot(ones(test_n,1)*resp(16), 'linewidth',2)
    plot(ones(test_n,1)*(resp(16) + resp_std(16)),'--r')
    plot(ones(test_n,1)*(resp(16) - resp_std(16)),'--r')
    title(['cone ', int2str(cone1),' vs cone ', int2str(cone2), ', n = ', int2str(n_inps)])
    drawnow
    cnt=cnt+1;
    
end



% for fake pairs
cone_pairs = [5, 1; 4, 10; 3, 8; 1, 10; 1, 9];
cnt=1;
test_n = 2;
clear ttt ttt_ste
null_second_cone = false
same_n = false
figure
for ii = 1:length(cone_pairs)
    cone1 = cone_pairs(ii,1);
    cone2 = cone_pairs(ii,2);
    tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
    n_gs=floor(length(tt)/10);
    tmp=sort(tt);
    my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(tt)+1;
    % spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
    resp_length = 30;
    pre_spike = 15;
    i=length(my_bins)-1;
    % 2 inputs together
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp = intersect(tmp, tmp1);
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
    end
    n_inps = length(tmp);
        
        
    for kk=1:test_n
        % input 1 separately
        tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone2,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp1=zeros(resp_length,1);
        for j=1:resp_length
            resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std1(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        % input 2 separately
        tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone1,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone1,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp2=zeros(resp_length,1);
        for j=1:resp_length
            resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std2(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        ttt(kk)=resp1(16)+resp2(16);
        ttt_ste(kk)=mean(resp_std1(16)+resp_std2(16));
    end
    
    
    subplot(3,4,cnt)
    plot(ttt, 'linewidth',2)
    hold on
    plot(ttt + ttt_ste, '--b')
    plot(ttt - ttt_ste, '--b')
    plot(ones(test_n,1)*resp(16), 'linewidth',2)
    plot(ones(test_n,1)*(resp(16) + resp_std(16)),'--r')
    plot(ones(test_n,1)*(resp(16) - resp_std(16)),'--r')
    title(['cone ', int2str(cone1),' vs cone ', int2str(cone2), ', n = ', int2str(n_inps)])
    drawnow
    cnt=cnt+1;
    
    cone2 = cone_pairs(ii,1);
    cone1 = cone_pairs(ii,2);
    tt = [my_filtered_inputs(cone1,:) my_filtered_inputs(cone2,:)];
    n_gs=floor(length(tt)/10);
    tmp=sort(tt);
    my_bins=[min(tt) tmp(n_gs:n_gs:n_gs*10)];
    my_bins(end)=max(tt)+1;
    % spike_rate_tmp = spike_rate_tmp-mean(spike_rate_tmp);
    resp_length = 30;
    pre_spike = 15;
    i=length(my_bins)-1;
    % 2 inputs together
    tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
    tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
    tmp = intersect(tmp, tmp1);
    tmp(tmp<(pre_spike+1)) = [];
    tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
    if ~isempty(tmp)
        resp=zeros(resp_length,1);
        for j=1:resp_length
            resp(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
    end
    n_inps = length(tmp);
        
    for kk=1:test_n
        % input 1 separately
        tmp = find(my_filtered_inputs(cone1,:)>=my_bins(i) & my_filtered_inputs(cone1,:)<my_bins(i+1));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone2,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone2,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp1=zeros(resp_length,1);
        for j=1:resp_length
            resp1(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std1(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        % input 2 separately
        tmp = find(my_filtered_inputs(cone2,:)<my_bins(end-i+1) & my_filtered_inputs(cone2,:)>=my_bins(end-i));
        if null_second_cone
            tmp1 = find(my_filtered_inputs(cone1,:)<my_bins(ceil(length(my_bins)/2)) & my_filtered_inputs(cone1,:)>=my_bins(ceil(length(my_bins)/2)-1));
            tmp = intersect(tmp, tmp1);
        end
        tmp(tmp<(pre_spike+1)) = [];
        tmp(tmp>=(length(spike_rate_tmp)+pre_spike-resp_length)) = [];
        if same_n
            tmp = sort(tmp(randperm(length(tmp),min(n_inps, length(tmp)))));
        end
        resp2=zeros(resp_length,1);
        for j=1:resp_length
            resp2(j) = mean(spike_rate_tmp(tmp-pre_spike+j-1)); % step back by pre_spike frames!
            resp_std2(j) = std(spike_rate_tmp(tmp-pre_spike+j-1))/sqrt(length(tmp));
        end
        
        ttt(kk)=resp1(16)+resp2(16);
        ttt_ste(kk)=mean(resp_std1(16)+resp_std2(16));
    end
    
    subplot(3,4,cnt)
    plot(ttt, 'linewidth',2)
    hold on
    plot(ttt + ttt_ste, '--b')
    plot(ttt - ttt_ste, '--b')
    plot(ones(test_n,1)*resp(16), 'linewidth',2)
    plot(ones(test_n,1)*(resp(16) + resp_std(16)),'--r')
    plot(ones(test_n,1)*(resp(16) - resp_std(16)),'--r')
    title(['cone ', int2str(cone1),' vs cone ', int2str(cone2), ', n = ', int2str(n_inps)])
    drawnow
    cnt=cnt+1;
    
end


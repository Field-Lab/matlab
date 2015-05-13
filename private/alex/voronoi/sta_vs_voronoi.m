%% load stuff
datarun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data008-from-d08_11/data008-from-d08_11');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = load_neurons(datarun);
datarun = set_polarities(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data011-from-d08_11/data011-from-d08_11');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);


wn_movie_name = 'BW-2-6-0.48-11111-300x300-60.35.xml';
[inputs, refresh, duration] = get_wn_movie_ath(datarun, wn_movie_name);

vormap = load('/Volumes/Data/2011-12-13-2/Visual/2011-12-13-2_f04_vorcones/map-0000.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2011-12-13-2/d08-11-norefit/data009-from-d08_11/data009-from-d08_11');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-11111-937x1-60.35.xml');

[inputs_fake, refresh_fake, duration_fake] = get_wn_movie_ath(vorrun, 'BW-1-2-0.48-22222-100x1-60.35.xml');

datarunID = 71;
visionID = datarun.cell_ids(datarunID);


%% recalculate sta


offset = 0;
sta_length=4;


starun = datarun;
starun = load_neurons(starun);


sta = squeeze(starun.stas.stas{datarunID});
figure
colormap gray
imagesc(sta(:,:,5))

aind = 150:250;
bind = 100:180;

% aind = 200:290;
% bind = 120:200;
cnt = 1;
myinds = [];
for i = aind
    for j = bind
        myinds(cnt) = sub2ind([300,300], i,j);
        cnt=cnt+1;
    end
end

a = reshape(sta(:,:,5), 300*300,1);
a = a(myinds);
a = reshape(a,length(bind),length(aind))';
figure
colormap gray
imagesc(a)



   
[type_name, my_type] = find_cell_type(starun, visionID);
   
spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/refresh); % spikes in frames
figure
cnt = 1;
spike_rate=zeros(duration,1);
my_sta=zeros(length(myinds),sta_length);
interv = 500;
for i=0:interv:interv*5
    spike_tmp = spikes(spikes>i &spikes<=i+interv);
    spike_tmp(spike_tmp<sta_length-offset)=[];
    tic
    while ~isempty(spike_tmp)
        [c, ia, ic] = unique(spike_tmp);
        spike_rate(spike_tmp(ia))=spike_rate(spike_tmp(ia))+1;
        for j=1:sta_length
            my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(inputs(myinds,spike_tmp(ia)-sta_length+j+offset),2);
        end
        spike_tmp(ia)=[];
    end
    toc
%     my_sta=my_sta/sum(spike_rate);
    tmp_sta = reshape(my_sta,length(bind),length(aind),4);

    subplot(4,6,cnt)
    colormap gray
    imagesc(tmp_sta(:,:,1)')
    drawnow
    title([int2str((i+interv)/(1000/refresh)), ' s'])
    cnt = cnt+1;
end
 


%% voronoi



figure
imagesc(vormap)

tmp_map = vormap(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2);
figure
imagesc(tmp_map)
myinds_vor = unique(tmp_map);
myinds_vor(myinds_vor==0)=[];
% myinds_vor = 1:937;
myinds_vor = myinds_vor';

sta_params.length = 6;
sta_params.offset = 0;


figure

my_sta=zeros(length(myinds_vor),sta_params.length);

spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
spike_rate=zeros(size(inputs_v,2),1);
k=0;
interv = 2700;
cnt1 = 7;
while cnt1<25
    k = k+interv;
    spikes_tmp = spikes(spikes>k & spikes <=k+interv);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
%     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3);
    
    subplot(4,6,cnt1)
    colormap gray
    imagesc(vorsta)
    title([int2str((k+interv)/(1000/refresh_v)), ' s'])
    cnt1 = cnt1+1;
    drawnow
end



%% crap cells


length(datarun.cell_types{1, 12}.cell_ids)

figure
cnt1 = 1;
for kk = 101:125    
    
    visionID = datarun.cell_types{1, 12}.cell_ids(kk);
    datarunID = find(datarun.cell_ids==visionID);
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes;
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
    %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3);
    
    subplot(5,5,cnt1)
    colormap gray
    imagesc(vorsta)
    title([int2str(visionID)])
    
    drawnow
    cnt1=cnt1+1;
end


% plot on parasols

figure
cnt1 = 1;
for visionID = [665 1308 1323 1922 2178 2357 3061 3558 3797 4052 4158 4369 4472 6542 6886 7355 6991 6482]  
    
%     visionID = datarun.cell_ids(datarunID);
    datarunID = find(datarun.cell_ids==visionID);
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes;
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
    %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3);
    
    subplot(4,5,cnt1)
    colormap gray
    imagesc(vorsta)
    title([int2str(visionID)])
    
    drawnow
    cnt1=cnt1+1;
end


% plot on midgets

figure
cnt1 = 1;

frames_to_use = 18000;

for visionID = datarun.cell_types{4}.cell_ids  
    
%     visionID = datarun.cell_ids(datarunID);
    datarunID = find(datarun.cell_ids==visionID);
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes(spikes<frames_to_use);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
    %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3);
    
    subplot(4,5,cnt1)
    colormap gray
    imagesc(vorsta)
    title([int2str(visionID)])
    
    drawnow
    cnt1=cnt1+1;
end



% compare time
first_10=9000;

second_10=54000-9000;

figure
cnt1 = 1;
for visionID = 1322%3736%2357  
    
%     visionID = datarun.cell_ids(datarunID);
    datarunID = find(datarun.cell_ids==visionID);
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    % first 10 minutes
    spikes_tmp = spikes(spikes<first_10);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    nsp = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
    %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3)/nsp;
    
    vor_first_10 = vorsta;
    
    % last 10 minutes
    spikes_tmp = spikes(spikes>second_10);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    nsp = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    vorsta=zeros(600,600,sta_params.length);
    coneX = zeros(max(tmp_map(:)),1);
    coneY = coneX;
    cnt = 1;
    for i=myinds_vor'
        [a, b] = find(tmp_map==tt(i+1));
        % find center of this cone
        coneX(i) = mean(a);
        coneY(i) = mean(b);
        if ~isempty(a)
            for j = 1:length(a)
                vorsta(a(j),b(j),:) = my_sta(cnt,:);
            end
        end
        cnt=cnt+1;
    end
    %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
    vorsta = vorsta(:, :, 3)/nsp;
    
    vor_second_10 = vorsta;
    
    tt = vor_first_10;
%     tt=-tt;
%     tt(tt<0) = 0;
% %     tt = tt+0.5;
% %     tt = -tt+0.5;
%     tt = tt/max(tt(:));
    
    tt1 = vor_second_10;
%     tt1=-tt1;
%     tt1(tt1<0) = 0;
% %     tt1 = tt1+0.5;
% %     tt1 = -tt1+0.5;
%     tt1 = tt1/max(tt1(:));
    
    tt = vor_first_10;
    tt1 = vor_second_10;
    figure
    plot(sort(tt(:))/min(tt(:)))
    hold on
    plot(sort(tt1(:))/min(tt1(:))*1)
    
    ttt = tt/min(tt(:))-tt1/min(tt1(:))*1; 
    ttt1 = ttt;
    ttt1(ttt1<0) = 0;
    
    ttt2 = -ttt;
    ttt2(ttt2<0) = 0;
    
    ttt_orig = tt/min(tt(:));
    
    tmp = zeros(600,600,3);
    tmp(:,:,1) = ttt1;
    tmp(:,:,2) = ttt_orig;
    tmp(:,:,3) = ttt2;
    
    
    figure
    imagesc(tmp)
    title(int2str(visionID))
    
    
%     tmp = zeros(600,600,3);
%     tmp(:,:,1) = tt;
%     tmp(:,:,3) = tt1;
    
    
    
    figure
    colormap gray
    imagesc(tt)
    figure
    colormap gray
    imagesc(tt1)
   
end






% estimate error for off midgets

figure
cnt1 = 1;

frames_to_use = 4500;

for visionID = datarun.cell_types{4}.cell_ids  
    
    
%     visionID = datarun.cell_ids(datarunID);
    datarunID = find(datarun.cell_ids==visionID);
    
    % first 5 minutes
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes(spikes<frames_to_use);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    nsp = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    my_sta=my_sta/nsp;
    
    temp_sta = my_sta(:,3);
    
    % last 5 minutes
    
    my_sta=zeros(length(myinds_vor),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes(spikes>54000-frames_to_use);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    nsp = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(myinds_vor,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    my_sta=my_sta/nsp;
    
    temp_sta2 = my_sta(:,3);

    
    % fake sta
    
    my_sta=zeros(100,sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    spikes_tmp = spikes(spikes<frames_to_use);
    spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
    nsp = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_fake(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    my_sta = my_sta/nsp;
    my_sta = my_sta(:,3);
    thresh = 2*std(my_sta)
    
    figure
    [~, ic] = sort(-temp_sta);
    plot(-temp_sta*1.15)
    hold on
    plot(-temp_sta2)
    line([1 937], [thresh thresh], 'color', 'k');
    
    figure
    plot(-temp_sta2(ic) + temp_sta(ic))
    hold on
    line([1 937], [thresh thresh], 'color', 'k');
    line([1 937], [-thresh -thresh], 'color', 'k');  
    
end


%% plots for presentation and ej

% first do Voronoi - simply 3rd frame, no temporal summation
sta_params.length = 5;
sta_params.offset = 0;
offset = 0;
sta_length=4;
filepath = '/Users/alexth/Desktop/voronoi/2011-12-13-2/d08-11-norefit/';

for datarunID = 1:length(datarun.cell_ids)
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 24 1886 1074]);
    
    my_sta=zeros(size(inputs_v,1),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    k = [sta_params.length-sta_params.offset, 3600,18000,54000];

    for cnt1 = 1:3
        spikes_tmp = spikes(spikes>k(cnt1) & spikes <=k(cnt1+1));
        
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
            for j=1:sta_params.length
                my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                    sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
            end
            spikes_tmp(ia)=[];
        end
        % my_sta=my_sta/nspikes;
        
        % get voronoi map sta
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        vorsta=zeros(600,600,sta_params.length);
        coneX = zeros(max(tmp_map(:)),1);
        coneY = coneX;
        cnt = 1;
        for i=1:937
            [a, b] = find(tmp_map==tt(i+1));
            % find center of this cone
            coneX(i) = mean(a);
            coneY(i) = mean(b);
            if ~isempty(a)
                for j = 1:length(a)
                    vorsta(a(j),b(j),:) = my_sta(cnt,:);
                end
            end
            cnt=cnt+1;
        end
        %     vorsta = vorsta(aind(1)*2:aind(end)*2, bind(1)*2:bind(end)*2, 3);
        vorsta_tmp = vorsta(:, :, 3);
        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(vorsta_tmp)
        title([int2str((k(cnt1+1))/(1000/refresh_v)), ' s'])
        drawnow
    end
    

    val = max(abs(my_sta(:,3)));    
    ind = find(abs(vorsta_tmp)==val, 1);
    [a,b] = ind2sub([600,600], ind);
    a = round(a/2);
    b = round(b/2);
    
    myinds = [];
    my_as = a-50:a+50;
    my_as(my_as<1 | my_as>300) = [];
    my_bs = b-50:b+50;
    my_bs(my_bs<1 | my_bs>300) = [];
    for i = my_as
        for j=my_bs
            myinds = [myinds sub2ind([300, 300],i, j)];
        end
    end
    
    % do single cone WN
       
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
    
    spike_rate=zeros(duration,1);
    my_sta=zeros(length(myinds),sta_length);
    k = [sta_length-offset, 1200, 6000, 18000];

    for cnt1 = 4:6
        spike_tmp = spikes(spikes>k(cnt1-3) &spikes<=k(cnt1-2));
        spike_tmp(spike_tmp<sta_length-offset)=[];
        tic
        while ~isempty(spike_tmp)
            [c, ia, ic] = unique(spike_tmp);
            spike_rate(spike_tmp(ia))=spike_rate(spike_tmp(ia))+1;
            for j=1:sta_length
                my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(inputs(myinds,spike_tmp(ia)-sta_length+j+offset),2);
            end
            spike_tmp(ia)=[];
        end
        toc
        %     my_sta=my_sta/sum(spike_rate);
        tmp_sta = reshape(my_sta,length(my_bs),length(my_as),4);
        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(tmp_sta(:,:,2)')
        drawnow
        title([int2str((k(cnt1-2))/(1000/refresh)), ' s'])
    end

    subplot(2,3,2)
    title(['Vision ID ', int2str(visionID), ', datarun ID ', int2str(datarunID), ', ', folder])
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarunID)]));
    close(fig)

end


%% combined RF maps (mosaics) for cell types

% voronoi
sta_params.length = 3;
sta_params.offset = 0;
clear thrs
clear my_coords

thrs{1} = zeros(1, length(datarun.cell_types{1}.cell_ids))+3.5;
thrs{1}([7 11 12 22]) = 3.7;
thrs{1}(17) = 7;
thrs{1}(18) = 3.1;

thrs{2} = zeros(1, length(datarun.cell_types{2}.cell_ids))+3.5;
thrs{2}(5) = 3.7;
thrs{2}([10 15 22 25]) = 4;
thrs{2}(21) = 4.5;

thrs{3} = zeros(1, length(datarun.cell_types{3}.cell_ids))+4;
thrs{3}([9 28 55 57 64]) = 4.5;
thrs{3}([41 61]) = 5.5;

thrs{4} = zeros(1, length(datarun.cell_types{4}.cell_ids))+4;
thrs{4}([12 67 71]) = 4.5;

for cell_type = 1:4
    cell_type
    
    for my_cell = 1:length(datarun.cell_types{cell_type}.cell_ids)
        
        visionID = datarun.cell_types{cell_type}.cell_ids(my_cell);
        datarunID = find(datarun.cell_ids == visionID);
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
        
        % fake voronoi seed
        
        my_sta=zeros(size(inputs_fake,1),sta_params.length);
        spike_rate=zeros(size(inputs_fake,2),1);
        
        spikes_tmp = spikes;
        spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
        nspikes = length(spikes_tmp);
        
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
            for j=1:sta_params.length
                my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                    sum(inputs_fake(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
            end
            spikes_tmp(ia)=[];
        end
        my_sta=my_sta(:,3)/nspikes;
        
        %      fake_most = 4.5*std(my_sta);
        fake_most = thrs{cell_type}(my_cell)*std(my_sta);
        
        
        
        % real voronoi seed
        my_sta=zeros(size(inputs_v,1),sta_params.length);
        spike_rate=zeros(size(inputs_v,2),1);
        
        spikes_tmp = spikes;
        spikes_tmp(spikes_tmp<sta_params.length-sta_params.offset)=[];
        
        while ~isempty(spikes_tmp)
            [~, ia, ~] = unique(spikes_tmp);
            spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
            for j=1:sta_params.length
                my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                    sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
            end
            spikes_tmp(ia)=[];
        end
        my_sta=my_sta(:,3)/nspikes;
        if datarun.stas.polarities{datarunID}<0
            my_sta = -my_sta;
        end
        
        %     figure
        %     plot(my_sta)
        real_most = find(my_sta>fake_most);
        
        % full map voronoi sta
        
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        vorsta_tmp=zeros(600,600);
        for i=1:937
            [a, b] = find(tmp_map==tt(i+1));
            if ~isempty(a)
                for j = 1:length(a)
                    vorsta_tmp(a(j),b(j),:) = my_sta(i);
                end
            end
        end
        
        
        % partial voronoi map sta
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        vorsta=zeros(600,600);
        for i=real_most'
            [a, b] = find(tmp_map==tt(i+1));
            if ~isempty(a)
                for j = 1:length(a)
                    vorsta(a(j),b(j)) = my_sta(i);
                end
            end
        end
        
        %     figure
        %     colormap gray
        %     imagesc(vorsta)
        [x,y] = find(vorsta);
        
        a=convhull(x,y);
        %      figure
        %      colormap gray
        %      imagesc(vorsta_tmp)
        %      hold on
        %      plot(y(a),x(a), '-r')
        %      drawnow
        
        my_coords{cell_type}{my_cell} = [y(a),x(a)];
        %
        %      [x1,y1] = find(vorsta_tmp);
        %
        %      [IN ON] = inpolygon(x1,y1,...
        %         x,y);
        %     IN = sort([find(IN); find(ON)]);
        %     x1(IN)
        
    end
    
end


figure

for cell_type = 1:4
    subplot(2,2,cell_type)
    axis ij
    
    for my_cell = 1:length(datarun.cell_types{cell_type}.cell_ids)
        hold on
        plot(my_coords{cell_type}{my_cell}(:,1),my_coords{cell_type}{my_cell}(:,2))
    end
    title(datarun.cell_types{cell_type}.name)
end

save('/Users/alexth/Desktop/voronoi/2011-12-13-2/mosaics_coordinates.mat', 'my_coords')


%% make video


% first do Voronoi - simply 3rd frame, no temporal summation
sta_params.length = 3;
sta_params.offset = 0;
offset = 0;
sta_length=2;

for datarunID = 45:311
    
    visionID = datarun.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun, visionID);
    
    if ~strcmp(folder,'crap') && ~strcmp(folder,'duplicates')
        datarunID
        my_sta=zeros(size(inputs_v,1),sta_params.length);
        
        spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
        spike_rate=zeros(size(inputs_v,2),1);
        
        k = [sta_params.length-sta_params.offset, 300:300:54000];
        vorsta_tmp = zeros(600,600,length(k));
        clear nsp
        for cnt1 = 1:length(k)-1
            cnt1
            spikes_tmp = spikes(spikes>k(cnt1) & spikes <=k(cnt1+1));
            nsp(cnt1) = length(spikes_tmp);
            while ~isempty(spikes_tmp)
                [~, ia, ~] = unique(spikes_tmp);
                spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
                for j=1:sta_params.length
                    my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                        sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
                end
                spikes_tmp(ia)=[];
            end
            % my_sta=my_sta/nspikes;
            
            % get voronoi map sta
            tmp_map = vormap;
            tt=0:max(tmp_map(:));
            vorsta=zeros(600,600,sta_params.length);
            cnt = 1;
            for i=1:937
                [a, b] = find(tmp_map==tt(i+1));
                if ~isempty(a)
                    for j = 1:length(a)
                        vorsta(a(j),b(j),:) = my_sta(cnt,:);
                    end
                end
                cnt=cnt+1;
            end
            vorsta_tmp(:,:,cnt1) = vorsta(:, :, 3);
        end
        
        save(['/Users/alexth/Desktop/voronoi/2011-12-13-2/sta_time/voronoi_',int2str(datarunID)],'vorsta_tmp', 'nsp');
    end
end


writerObj = VideoWriter('/Users/alexth/Desktop/video89.avi');
writerObj.FrameRate = 3;
writerObj.Quality=100;
open(writerObj);

figure
for j=1:length(k)-1
    
    colormap('gray')
    tmp = vorsta_tmp(:,:,j);
    imagesc(tmp)
    s = k(j+1)*refresh_v/1000;
    m = floor(s/60);
    if m>0
        s = mod(s,60);
        str  = [int2str(m), ' min  ', int2str(s), ' s'];
    else
        str  = [ int2str(s), ' s'];
    end
    text(295,50, str, 'color',[0.9 0.5 0.5], 'fontsize',28)
    movieFrames(j)=getframe;
    writeVideo(writerObj,movieFrames(j));
end
close(writerObj);



for datarunID = 1:311
    visionID = datarun.cell_ids(datarunID);
    
    % full voronoi    
    my_sta=zeros(size(inputs_v,1),sta_params.length);    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    spikes_tmp = spikes(spikes>sta_params.length-sta_params.offset);
    full_voronoi_sta = zeros(600,600);
    nspikes = length(spikes_tmp);
    while ~isempty(spikes_tmp)
        [~, ia, ~] = unique(spikes_tmp);
        spike_rate(spikes_tmp(ia))=spike_rate(spikes_tmp(ia))+1;
        for j=1:sta_params.length
            my_sta(:,sta_params.length-j+1) = my_sta(:, sta_params.length-j+1) +...
                sum(inputs_v(:,spikes_tmp(ia)-sta_params.length+j+sta_params.offset),2);
        end
        spikes_tmp(ia)=[];
    end
    % my_sta=my_sta/nspikes;
    
    % get voronoi map sta
    tmp_map = vormap;
    tt=0:max(tmp_map(:));
    for i=1:937
        [a, b] = find(tmp_map==tt(i+1));
        if ~isempty(a)
            for j = 1:length(a)
                full_voronoi_sta(a(j),b(j)) = my_sta(i,3);
            end
        end
    end
    
    
    val = max(abs(my_sta(:,3)));
    ind = find(abs(full_voronoi_sta)==val, 1);
    [a,b] = ind2sub([600,600], ind);
    a = round(a/2);
    b = round(b/2);
    
    border = 50;
    myinds = [];
    my_as = a-border:a+border;
    my_as(my_as<1 | my_as>300) = [];
    my_bs = b-30:b+border;
    my_bs(my_bs<1 | my_bs>300) = [];
    for i = my_as
        for j=my_bs
            myinds = [myinds sub2ind([300, 300],i, j)];
        end
    end
    
    % do single cone WN
    
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
    
    spike_rate=zeros(duration,1);
    my_sta=zeros(length(myinds),sta_length);
    ksta = [sta_params.length-sta_params.offset, 100:100:18000];
    sta_tmp = zeros(border*2+1,border*2+1,length(ksta)-1);
    
    for cnt1 = 1:length(ksta)-1
        cnt1
        spike_tmp = spikes(spikes>ksta(cnt1) &spikes<=ksta(cnt1+1));
        spike_tmp(spike_tmp<sta_length-offset)=[];
        nsp(cnt1) = length(spikes_tmp);
        while ~isempty(spike_tmp)
            [c, ia, ic] = unique(spike_tmp);
            spike_rate(spike_tmp(ia))=spike_rate(spike_tmp(ia))+1;
            for j=1:sta_length
                my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(inputs(myinds,spike_tmp(ia)-sta_length+j+offset),2);
            end
            spike_tmp(ia)=[];
        end
        tmp = reshape(my_sta,length(my_bs),length(my_as),2);
        
        sta_tmp(:,:,cnt1) = tmp(:,:,2)';
    end
    
    save(['/Users/alexth/Desktop/voronoi/2011-12-13-2/sta_time/single_cone_',int2str(datarunID)],'sta_tmp', 'nsp');
    
end


writerObj = VideoWriter('/Users/alexth/Desktop/voronoi/2011-12-13-2/sta_videos/video168_fast.avi');
writerObj.FrameRate = 10;
writerObj.Quality=100;
open(writerObj);

f = figure
set(gcf, 'position', [147         590        1262         509])
for j=1:length(k)-1
    
    subplot(1,2,1)
    colormap('gray')
    tmp = sta_tmp(:,:,j);
    imagesc(tmp)
    s = ksta(j+1)*refresh/1000;
    m = floor(s/60);
    if m>0
        s = mod(s,60);
        str  = [int2str(m), ' min  ', int2str(s), ' s'];
    else
        str  = [ int2str(s), ' s'];
    end
    text(25,5, str, 'color',[0.9 0.5 0.5], 'fontsize',28)
    
    subplot(1,2,2)
    colormap('gray')
    tmp = vorsta_tmp(my_as*2,my_bs*2,j);
    imagesc(tmp)
    s = k(j+1)*refresh_v/1000;
    m = floor(s/60);
    if m>0
        s = mod(s,60);
        str  = [int2str(m), ' min  ', int2str(s), ' s'];
    else
        str  = [ int2str(s), ' s'];
    end
    text(25,5, str, 'color',[0.9 0.5 0.5], 'fontsize',28)    
    
    movieFrames(j)=getframe(f);
    writeVideo(writerObj,movieFrames(j));
end
close(writerObj);


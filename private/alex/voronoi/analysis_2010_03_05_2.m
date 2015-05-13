%% rough stability check
datarun = load_data('/Volumes/Analysis/2010-03-05-2/data012/data012');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-03-05-2/data015/data015');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun1, {4}, 'clear', false, 'plot_fits', true, 'fit_color', 'r')


%% load stuff
datarun = load_data('/Volumes/Analysis/2010-03-05-2/d12-15-norefit/data013-from-d12_15/data013-from-d12_15');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2010-03-05-2/d12-15-norefit/data015-from-d12_15/data015-from-d12_15');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

% wn_movie_name = 'RGB-1-8-0.48-33333-320x320.xml';
wn_movie_name = 'BW-1-8-0.48-11111-320x320.xml';
[inputs, refresh, duration] = get_wn_movie_ath(datarun1, wn_movie_name);


vormap = load('/Volumes/Archive/2010-03-05-2/Visual/cone_map-data-014.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2010-03-05-2/d12-15-norefit/data014-from-d12_15/data014-from-d12_15');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);

figure
plot(diff(vorrun.triggers))

[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-1-8-0.48-11111-1517x1.xml');

%% find stable cells

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

cell_type = 2;
bords = 15;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(320,320,3);
tmp(:,:,2) = tmp_map;
for i=1:9%find(ww./pixn<5)
    visionID = datarun.cell_types{cell_type}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    sta = double(datarun.stas.stas{cellInd}(:,:,2,4));
    sta1 = double(datarun1.stas.stas{cellInd}(:,:,1,4));
    
    if datarun.stas.polarities{cellInd} == 1
        [tt, ind] = max(sta(:));
        [tt1, ~] = max(sta1(:));
    else
        [tt, ind] = min(sta(:));
        [tt1, ~] = min(sta1(:));
    end
    [a,b] = ind2sub([320,320],ind);
    
    sta = sta/tt;
    sta1 = sta1/tt1;

    tmp(:,:,1) = sta;
    tmp(:,:,3) = sta1; 
    subplot(3,3,cnt)
    imagesc(tmp)
    axis([b-bords b+bords a-bords a+bords])    
    title(['cell ID ', int2str(cellInd), ', vision ID ', int2str(visionID), ', i=', int2str(i)])
    
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




%% plots for presentation and ej

% first do Voronoi - simply 3rd frame, no temporal summation
sta_params.length = 5;
sta_params.offset = 0;
sta_length = 3;
offset = 0;
filepath = '/Users/alexth/Desktop/voronoi/2010-03-05-2/d12-15-norefit/';
starun = datarun1;

border = 25;

for datarunID = 1:length(datarun1.cell_ids)
    visionID = datarun1.cell_ids(datarunID);
    [folder, my_type] = find_cell_type(datarun1, visionID);
    
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 24 1886 1074]);
    
    my_sta=zeros(size(inputs_v,1),sta_params.length);
    
    spikes=ceil((vorrun.spikes{datarunID}-vorrun.triggers(1))*1000/(refresh_v)); % spikes in frames
    spike_rate=zeros(size(inputs_v,2),1);
    
    k = -9000; % refresh 66ms, 5 min (300s) is about 4500 spike-frames
    interv = 9000;
    cnt1 = 1;
    while cnt1<4
        k = k+interv;
        spikes_tmp = spikes(spikes>k & spikes <=k+interv);
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
        % my_sta=my_sta/nspikes;
        
        % get voronoi map sta
        tmp_map = vormap;
        tt=0:max(tmp_map(:));
        vorsta=zeros(320,320,sta_params.length);
        coneX = zeros(max(tmp_map(:)),1);
        coneY = coneX;
        cnt = 1;
        for i=1:1517
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
        vorsta_tmp = vorsta(:, :, 2);
        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(vorsta_tmp)
        title([int2str((k+interv)/(1000/refresh_v)), ' s'])
        cnt1 = cnt1+1;
        drawnow
    end
    

    val = max(abs(my_sta(:,2)));    
    ind = find(abs(vorsta_tmp)==val, 1);
    [a,b] = ind2sub([320,320], ind);
    
    myinds = [];
    my_as = a-border:a+border;
    my_as(my_as<1 | my_as>320) = [];
    my_bs = b-border:b+border;
    my_bs(my_bs<1 | my_bs>320) = [];
    for i = my_as
        for j=my_bs
            myinds = [myinds sub2ind([320, 320],i, j)];
        end
    end
    
    subplot(2,3,2)
    title(['Vision ID ', int2str(visionID), ', datarun ID ', int2str(datarunID), ', ', folder])
  
    
    % do single cone WN
       
    spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/refresh); % spikes in frames
    
    spike_rate=zeros(duration,1);
    my_sta=zeros(length(myinds),sta_length);
    k = -9000;
    interv = 9000;
    cnt1 = 4;
    while cnt1<7
        k = k+interv;
        spike_tmp = spikes(spikes>k &spikes<=k+interv);
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
        tmp_sta = reshape(my_sta,length(my_bs),length(my_as),3);
        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(tmp_sta(:,:,2)')
        drawnow
        title([int2str((k+interv)/(1000/refresh)), ' s'])
        cnt1 = cnt1+1;
    end

     
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarunID)]));
    close(fig)

end


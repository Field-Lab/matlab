%% load stuff
datarun = load_data('/Volumes/Analysis/2011-10-25-0/data001/data001');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-10-25-0/streamed/data006-0/data006-0');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

figure
plot_rf_summaries(datarun, {4}, 'clear', false, 'label', true, 'label_color', 'b', 'plot_fits', true, 'fit_color', 'b')
plot_rf_summaries(datarun1, {4}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


%% load stuff
datarun = load_data('/Volumes/Analysis/2011-10-25-0/d01-06-norefit/data002-from-d01_06/data002-from-d01_06');
datarun = load_params(datarun,'verbose',1);
datarun = load_sta(datarun);
datarun = set_polarities(datarun);
datarun = load_neurons(datarun);

datarun1 = load_data('/Volumes/Analysis/2011-10-25-0/d01-06-norefit/data006-from-d01_06/data006-from-d01_06');
datarun1 = load_params(datarun1,'verbose',1);
datarun1 = load_sta(datarun1);
datarun1 = set_polarities(datarun1);
datarun1 = load_neurons(datarun1);

wn_movie_name = 'BW-2-5-0.48-22222-300x300-60.35.xml';
[inputs, refresh, duration] = get_wn_movie_ath(datarun, wn_movie_name);


vormap = load('/Volumes/Data/2011-10-25-0/Visual/2011-20-25-0_f01_vorcones/map-0000.txt');
figure
imagesc(vormap)


vorrun = load_data('/Volumes/Analysis/2011-10-25-0/d01-06-norefit/data003-from-d01_06/data003-from-d01_06');
vorrun = load_params(vorrun,'verbose',1);
vorrun = load_neurons(vorrun);
[inputs_v, refresh_v, duration_v] = get_wn_movie_ath(vorrun, 'BW-gaussian-1-5-0.16-11111-1603x1-60.35.xml');


%% find stable cells

tmp_map = vormap;
for i=1:max(tmp_map(:))  
    tmp_map(tmp_map==i) = 0.3+(rand(1)-0.5)/8;
end

cell_type = 4;
offm = datarun.cell_types{4}.cell_ids;
bords = 30;
figure
set(gcf, 'position', [58 102 1182 996])
cnt=1;
tmp = zeros(600,600,3);
tmp(:,:,2) = tmp_map;

for i=18:27
    visionID = datarun.cell_types{4}.cell_ids(i);
    cellInd = find(datarun.cell_ids == visionID);
    sta = imresize(double(datarun.stas.stas{cellInd}(:,:,1,4)),600/datarun.stimulus.field_height,'nearest');
    sta1 = imresize(double(datarun1.stas.stas{cellInd}(:,:,1,4)),600/datarun1.stimulus.field_height,'nearest');
    
    if datarun.stas.polarities{cellInd} == 1
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
    title(['cell ID ', int2str(cellInd), ', vision ID ', int2str(visionID), ', i=', int2str(i)])

    cnt = cnt+1;
end



%% 

% first do Voronoi - simply 3rd frame, no temporal summation
sta_params.length = 3;
sta_params.offset = 0;
offset = 0;
sta_length=4;
filepath = '/Users/alexth/Desktop/voronoi/2011-10-25-0/d01-06-norefit/';

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
    
    % 31600 frames is 2700s;
    k = [sta_params.length-sta_params.offset, 1500, 7200, 31641];

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
        for i=1:1603
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
        vorsta_tmp = vorsta(:, :, 2);
        
        subplot(2,3,cnt1)
        colormap gray
        imagesc(vorsta_tmp)
        title([int2str((k(cnt1+1))/(1000/refresh_v)), ' s'])
        drawnow
    end
    

    val = max(abs(my_sta(:,2)));    
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
    
    subplot(2,3,2)
    title(['Vision ID ', int2str(visionID), ', datarun ID ', int2str(datarunID), ', ', folder])

    
    % do single cone WN
       
    spikes=ceil((datarun.spikes{datarunID}-datarun.triggers(1))*1000/refresh); % spikes in frames
    
    spike_rate=zeros(duration,1);
    my_sta=zeros(length(myinds),sta_length);
    k = [sta_length-offset, 1500, 7200, 21741];

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

    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(datarunID)]));
    close(fig)

end


%% load NSEM
mvpath='/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/';
my_movie=zeros(320,320,3600);
cnt=1;
for i=1:30
    load([mvpath, 'movie_chunk_', int2str(i), '.mat']);
    movie=movie(81:240,:,:);
    movie=imresize(movie,2,'nearest');
    my_movie(:,:,cnt:cnt+119)=movie;
    cnt=cnt+120;
end
my_movie=my_movie/255-0.5;
clear mvpath movie

%% load dataruns
% 04,08 - NDF4 (09)
% 11,14 - NDF3 (12)
% 15,18 - NDF2 (16)
% 19,22 - NDF1 (20)
% 23,27 - NDF0 (24)

data='008';

for i=1
    switch data
        case '004'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            ndf=4;
            nmdata='009';
        case '008'
            wn_movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            ndf=4;
            nmdata='009';
        case '011'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            ndf=3;
            nmdata='012';
        case '014'
            wn_movie_name = 'BW-8-8-0.48-11111-40x40.xml';
            ndf=3;
            nmdata='012';
        case '015'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            ndf=2;
            nmdata='016';
        case '018'  
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            ndf=2;
            nmdata='016';
        case '019' 
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            ndf=1;
            nmdata='020';
        case '022' 
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            ndf=1;
            nmdata='020';
        case '023' 
            wn_movie_name = 'BW-16-4-0.48-11111-20x20.xml';
            ndf=0;
            nmdata='024';
        case '027' 
            wn_movie_name = 'BW-4-2-0.48-11111-80x80.xml';
            ndf=0;
            nmdata='024';
    end  
    starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
    starun = load_params(starun,'verbose',1);
    starun = set_polarities(starun);
    starun = load_neurons(starun);
    starun = load_sta(starun,'load_sta','all','keep_java_sta',true);

    nmrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',nmdata,'-from-d05-d27/data',nmdata,'-from-d05-d27']);
    nmrun = load_params(nmrun,'verbose',1);
    nmrun = load_neurons(nmrun);
end
clear data nmdata

%% process white noise movie
[inputs, refresh, duration] = get_wn_movie_ath(starun, wn_movie_name); 
clear wn_movie_name

%% preprocess NSEM (scale, average)
scaled_movie=zeros(size(inputs,1), 3600);
stix_size=320/sqrt(size(inputs,1));
cnt=1;
for i=1:stix_size:320
    for j=1:stix_size:320
        tmp=reshape(my_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end
clear stix_size


%% find mean time courses for cell types
tc=zeros(length(starun.cell_types),length(starun.vision.timecourses(1).g));
for j=1:size(tc,1)
    cnt=1;
    for i=starun.cell_types{j}.cell_ids
        tmp=find(starun.cell_ids==i);
        if ~isempty(starun.vision.timecourses(tmp).g)
            if sum(starun.vision.timecourses(tmp).g == starun.vision.timecourses(tmp).b)
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).g';
            elseif  max(starun.vision.timecourses(tmp).g)<max(starun.vision.timecourses(tmp).b) % for blue cells, RGB run
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).b';
            end
            cnt=cnt+1;
        end
    end
    tc(j,:)=tc(j,:)'/cnt;
end

%% calculate LF and NL for cells


figure

for my_type = 1:5
    sigstix_threshold = 4;
    
    downsize=1;
    offset = 0;
    
    sta_range=3:20;
    
    r2=zeros(length(starun.cell_types{my_type}.cell_ids),length(sta_range));
    inds_n=zeros(size(r2,1),1);
    cell_cnt=1;
    for cellID = starun.cell_types{my_type}.cell_ids
        cell_cnt
        datarunID=find(starun.cell_ids==cellID);
        
        % sta of a chosen cell
        sta = double(squeeze(starun.stas.stas{datarunID}));
        
        if ndims(sta)==4 % ONLY FOR GREEN CELLS!
            sta=squeeze(sta(:,:,2,:));
        end
        
        % get weights by convolving sta of a chosen cell with mean time course
        temp_sta = reshape(sta, size(sta,1)*size(sta,2),size(sta,3));
        weights=sum(repmat(tc(my_type,:),size(inputs,1),1) .* temp_sta, 2);
        weights = weights/max(weights);
        
        
        % threshold weights
        weight_threshold=robust_std(weights)*sigstix_threshold;
        inds=find(weights>weight_threshold);
        
        my_inputs = inputs(inds, :); % to save memory and time
        inp= zeros(length(inds),size(inputs,2)*downsize);
        % down to a frame: instead of refresh of 66.6192, use refresh/8 = 8.3274
        for i=1:downsize
            inp(:,i:downsize:end) = my_inputs;
        end
        my_inputs=inp;
        
        inds_n(cell_cnt)=length(inds);
        
        sta_length_cnt=1;
        for sta_length = sta_range
            
            spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/(refresh/downsize)); % spikes in ms
            spikes(spikes<sta_length-offset)=[];
            
            spike_rate=zeros(duration*downsize,1);
            my_sta=zeros(length(inds),sta_length);
            
            while ~isempty(spikes)
                [c, ia, ic] = unique(spikes);
                spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
                for j=1:sta_length
                    my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(my_inputs(:,spikes(ia)-sta_length+j+offset),2);
                end
                spikes(ia)=[];
            end
            my_sta=my_sta/sum(spike_rate);
            
            
            
            % convolve raw input with linear filter
            %         my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2));
            my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
            for i=1:length(inds)
                my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
            end
            
            gen_sig=sum(my_filtered_inputs);
            gen_sig=gen_sig/max(gen_sig);
            
            spike_rate = spike_rate(sta_length:end);
            
            bins = (max(gen_sig)-min(gen_sig))/100;
            my_bins = min(gen_sig):bins:max(gen_sig)-bins;
            
            % calculate nonlinearity
            my_nl=zeros(size(my_bins));
            cnt=1;
            for i=my_bins
                tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
                tmp(tmp<1)=[];
                tmp(tmp>length(spike_rate))=[];
                my_nl(cnt)=mean(spike_rate(tmp)); % find mean firing rate after this certain value
                cnt=cnt+1;
            end
            
            % calculate predicted FR
            cnt=1;
            predicted_rate=zeros(size(spike_rate));
            for i=min(gen_sig):bins:max(gen_sig)-bins
                tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
                tmp(tmp<1)=[];
                tmp(tmp>length(spike_rate))=[];
                predicted_rate(tmp)=my_nl(cnt);
                cnt=cnt+1;
            end
            sse = sum((predicted_rate-spike_rate).^2);
            sst = sum((spike_rate-mean(spike_rate)).^2);
            r2(cell_cnt, sta_length_cnt) = 1 - sse/sst;
            sta_length_cnt=sta_length_cnt+1;
        end
        max(max(r2(cell_cnt, :, :)))
        cell_cnt=cell_cnt+1;
    end
    
    subplot(2,3,my_type)
    
    plot(sta_range,r2', '-*')
    hold on
    [t,p]=max(nanmean(r2));
    plot(sta_range,nanmean(r2),'-*k','linewidth',5)
    %     axis([1 12 0 0.6])
    title([starun.cell_types{my_type}.name, ', max r2 ', num2str(t), ' at ', int2str(p+sta_range(1))])
end


%% use LF to calculate GS for NSEM

downsize=1;
offset = 0;
int_method='spline';
sta_length=15;
sigstix_threshold = 4;

nm_triggers=nmrun.triggers;
movie_refresh=refresh/8;%1/120*1000; %ms
movie_duration = 31000; % ms,single repetition
nm_repeats = 20;

figure
for my_type=1:5
    
    r2=zeros(length(starun.cell_types{my_type}.cell_ids),2);
    cell_cnt=1;
    for cellID = starun.cell_types{my_type}.cell_ids
        cell_cnt
        datarunID=find(starun.cell_ids==cellID);
        
        %****** PART 1 calculate linear filter from STA ******
        % sta of a chosen cell
        sta = double(squeeze(starun.stas.stas{datarunID}));
        
        if ndims(sta)==4 % ONLY FOR GREEN CELLS!
            sta=squeeze(sta(:,:,2,:));
        end
        
        % get weights by convolving sta of a chosen cell with mean time course
        temp_sta = reshape(sta, size(sta,1)*size(sta,2),size(sta,3));
        weights=sum(repmat(tc(my_type,:),size(sta,1)*size(sta,2),1) .* temp_sta, 2);
        weights = weights/max(weights);
        
        % threshold weights
        weight_threshold=robust_std(weights)*sigstix_threshold;
        inds=find(weights>weight_threshold);
        
        my_inputs = inputs(inds, :); % to save memory and time
        inp= zeros(length(inds),size(inputs,2)*downsize);
        % down to a frame: instead of refresh of 66.6192, use refresh/8 = 8.3274
        for i=1:downsize
            inp(:,i:downsize:end) = my_inputs;
        end
        my_inputs=inp;
        
        inds_n(cell_cnt)=length(inds);
        
        spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/(refresh/downsize)); % spikes in ms
        spikes(spikes<sta_length-offset)=[];
        
        spike_rate=zeros(duration*downsize,1);
        my_sta=zeros(length(inds),sta_length);
        
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            for j=1:sta_length
                my_sta(:,sta_length-j+1)=my_sta(:,sta_length-j+1)+sum(my_inputs(:,spikes(ia)-sta_length+j+offset),2);
            end
            spikes(ia)=[];
        end
        my_sta=my_sta/sum(spike_rate);
        
        % convolve raw input with linear filter
        my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
        for i=1:length(inds)
            my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
        end
        
        gen_sig=sum(my_filtered_inputs);
        gen_sig=gen_sig/max(gen_sig);
        
        spike_rate = spike_rate(sta_length:end);
        
        bins = (max(gen_sig)-min(gen_sig))/100;
        my_bins = min(gen_sig):bins:max(gen_sig)-bins;
        
        % calculate nonlinearity
        my_nl=zeros(size(my_bins));
        cnt=1;
        for i=my_bins
            tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
            tmp(tmp<1)=[];
            tmp(tmp>length(spike_rate))=[];
            my_nl(cnt)=mean(spike_rate(tmp)); % find mean firing rate after this certain value
            cnt=cnt+1;
        end
        
        % calculate predicted FR
        cnt=1;
        predicted_rate=zeros(size(spike_rate));
        for i=min(gen_sig):bins:max(gen_sig)-bins
            tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
            tmp(tmp<1)=[];
            tmp(tmp>length(spike_rate))=[];
            predicted_rate(tmp)=my_nl(cnt);
            cnt=cnt+1;
        end
        sse = sum((predicted_rate-spike_rate).^2);
        sst = sum((spike_rate-mean(spike_rate)).^2);
        r2(cell_cnt,1) = 1 - sse/sst;
        
        
        %****** PART 2 calculate gen signal from NSEM using STA ******
        
        spikes=ceil((nmrun.spikes{datarunID}-nmrun.triggers(1))*1000/movie_refresh); % spikes in ms
        spikes(spikes<sta_length-offset)=[];
        
        spike_rate=zeros(movie_duration*nm_repeats,1);
        
        while ~isempty(spikes)
            [c, ia, ic] = unique(spikes);
            spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
            spikes(ia)=[];
        end
        
        
        my_begs=round([0; nm_triggers(find(diff(nm_triggers)>0.9)+1)]*1000/movie_refresh);
        acc_spike_rate=0;
        
        for i=1:length(my_begs)-1
            acc_spike_rate = acc_spike_rate + spike_rate(my_begs(i)+1:my_begs(i)+3600);
        end
        acc_spike_rate=acc_spike_rate/(length(my_begs)-1);
        
        % interpolate linear filter to match movie frames
        if downsize==1
            tt=[];
            for i=1:size(my_sta,1)
                tt=[tt; interp1(1:8:8*size(my_sta,2), my_sta(i,:),1:8*size(my_sta,2),int_method)]; % spline, linear
            end
            my_sta=tt;
        end
        % convolve raw input with linear filter
        my_inputs=scaled_movie(inds, :);
        my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
        for i=1:length(inds)
            my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
        end
        
        gen_sig=sum(my_filtered_inputs);
        gen_sig=gen_sig/max(gen_sig);
        
        acc_spike_rate = acc_spike_rate(sta_length:end);
        
        bins = (max(gen_sig)-min(gen_sig))/100;
        my_bins = min(gen_sig):bins:max(gen_sig)-bins;
        
        % calculate nonlinearity
        my_nl=zeros(size(my_bins));
        cnt=1;
        for i=my_bins
            tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
            tmp(tmp<1)=[];
            tmp(tmp>length(acc_spike_rate))=[];
            my_nl(cnt)=mean(acc_spike_rate(tmp)); % find mean firing rate after this certain value
            cnt=cnt+1;
        end
        
        % calculate predicted FR
        cnt=1;
        predicted_rate=zeros(size(acc_spike_rate));
        for i=min(gen_sig):bins:max(gen_sig)-bins
            tmp=find(gen_sig>=i & gen_sig<i+bins)-offset; % find instances when GS had certain value
            tmp(tmp<1)=[];
            tmp(tmp>length(acc_spike_rate))=[];
            predicted_rate(tmp)=my_nl(cnt);
            cnt=cnt+1;
        end
        sse = sum((predicted_rate-acc_spike_rate).^2);
        sst = sum((acc_spike_rate-mean(acc_spike_rate)).^2);
        r2(cell_cnt,2) = 1 - sse/sst;
        
        %     figure
        %     plot(predicted_rate)
        %     hold on
        %     plot(acc_spike_rate)
        
        cell_cnt=cell_cnt+1;
    end
    
    subplot(2,3,my_type)
    plot(r2)
    legend('WN', 'NSEM')
    title(starun.cell_types{my_type}.name)
end

%% compare GS for 2 NDFs

scale = true;
mean_rgb = false;

data='011';
[first_nd.psr, first_nd.asr, first_nd.gs, first_nd.r2, first_nd.ind, first_nd.cellID] = get_LN_ndf(data, scale, mean_rgb);
save(['/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/not_scaled/', 'data_',data,'.mat'], 'first_nd');

data='019';
[second_nd.psr, second_nd.asr, second_nd.gs, second_nd.r2, second_nd.ind, second_nd.cellID] = get_LN_ndf(data, scale, mean_rgb);
save(['/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/not_scaled/', 'data_019.mat'], 'second_nd');


nds={'NDF2','NDF1'};
%%%%%%%% plot it
data='011';
starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
starun = load_params(starun,'verbose',1);
filepath = '/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/not_scaled/ND2_1/';

cell_cnt=1;
for cellID = starun.cell_ids    
    
    datarunID=find(starun.cell_ids==cellID);
    [~, my_type] = find_cell_type(starun, cellID);

    folder=starun.cell_types{my_type}.name;
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 244 1637 854]);
    subplot(2,1,1)
    plot(first_nd.gs(cell_cnt,:))
    hold on
    plot(second_nd.gs(cell_cnt,:))
    title(['GS cell ', int2str(cellID), ', r2: ',nds{1},' = ', num2str(first_nd.r2(cell_cnt,2)), ...
        ', ',nds{2},' = ', num2str(second_nd.r2(cell_cnt,2))])
    legend(nds)
    for i=120:120:3600
        line([i,i],[min([first_nd.gs(cell_cnt,:) second_nd.gs(cell_cnt,:)]), max([first_nd.gs(cell_cnt,:) second_nd.gs(cell_cnt,:)])],'color','k')
    end
    axis([0 3600 -Inf Inf])
    
    subplot(2,1,2)
    plot(first_nd.asr(:,cell_cnt))
    hold on
    plot(second_nd.asr(:,cell_cnt))
    title(['Data cell ', int2str(cellID), ',  ', folder])
    legend(nds)
    for i=120:120:3600
        line([i,i],[0, max([first_nd.asr(:,cell_cnt); second_nd.asr(:,cell_cnt)])],'color','k')
    end
    axis([0 3600 -Inf Inf])
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(cellID)]));
    close(fig)
    
    cell_cnt=cell_cnt+1;
end



%%
mean_rgb = false;
even_bin = true;

data='008';
data1 = fit_LN_ndf(data, mean_rgb, even_bin);

data='0111';
data2 = fit_LN_ndf(data, mean_rgb, even_bin);

scales = zeros(length(data1.sta),2);
for cell_cnt=1:length(data1.sta)

    binnie1 = data1.bins{cell_cnt};
    binnie2 = data2.bins{cell_cnt};
    nl1 = data1.nl{cell_cnt};
    nl2 = data2.nl{cell_cnt};
    nl1 = nl1(2:99);
    nl2 = nl2(2:99);
    binnie1 = binnie1(2:99);
    binnie2 = binnie2(2:99);
    
    t = unique([find(isnan(nl1)) find(isnan(nl2))]);
    nl1(t) = [];
    nl2(t) = [];
    binnie1(t) = [];
    binnie2(t) = [];
    nl1 = nl1(end:-1:1);
    nl2 = nl2(end:-1:1);
    
    x = [binnie1; binnie2];
    x = x-max(x(:));
    y = [nl1 nl2]';
    
    [p, resnorm, residual] = fit_NL_to_NCDF(x,y);
    scales(cell_cnt,:) = p(1:end-2);
    
%     subplot(4,4,cell_cnt)
%     for i = 1:2
%         xt = x(i,:) .* p(i);
%         h(i) = plot(xt, y(i,:));
%         hold on
%     end    
%     for i = 1:2
%         xt = x(i,:);
%         h(i) = plot(xt, y(i,:));
%         hold on
%     end
%     title(num2str(scale_factor'))
    
end

data1.scale = scales(:,1);
data2.scale = scales(:,2);

data1 = fit_NSEM_ndf_raw(data1, even_bin);
data2 = fit_NSEM_ndf_raw(data2, even_bin);

scales_nm = zeros(length(data1.sta),2);
for cell_cnt=1:length(data1.sta)

    binnie1 = data1.nm.bins{cell_cnt}(1:end-1);
    binnie2 = data2.nm.bins{cell_cnt}(1:end-1);
    nl1 = data1.nm.nl{cell_cnt};
    nl2 = data2.nm.nl{cell_cnt};
    nl1 = nl1(2:99);
    nl2 = nl2(2:99);
    binnie1 = binnie1(2:99);
    binnie2 = binnie2(2:99);
    
    t = unique([find(isnan(nl1)) find(isnan(nl2))]);
    nl1(t) = [];
    nl2(t) = [];
    binnie1(t) = [];
    binnie2(t) = [];
    nl1 = nl1(end:-1:1);
    nl2 = nl2(end:-1:1);
    
    x = [binnie1; binnie2];
    x = x-max(x(:));
    y = [nl1'; nl2'];
    
    [p, resnorm, residual] = fit_NL_to_NCDF(x,y);
    scales_nm(cell_cnt,:) = p(1:end-2);
end

data1.nm.scale = scales_nm(:,1);
data2.nm.scale = scales_nm(:,2);

data1 = fit_WNrep_ndf_raw(data1, even_bin);

figure
plot(data1.wnrep.r2)
legend('long run','repeats')

figure
subplot(2,1,1)
tt = unique(data1.asr{15});
tmp = data1.asr{15}/tt(2);
tmp = conv(tmp,kern,'same')/50*tt(2);
plot(tmp,'r')

hold on
tt = unique(data1.psr{15});
tmp = data1.psr{15}/tt(2);
tmp = conv(tmp,kern,'same')/50*tt(2);
plot(tmp,'k')
axis([0 450 0 1.5])
title('long run')

subplot(2,1,2)
tt = unique(data1.wnrep.asr{15});
tmp = data1.wnrep.asr{15}/tt(2);
tmp = conv(tmp,kern,'same')/50*tt(2);
plot(tmp,'r')

hold on
tt = unique(data1.wnrep.psr{15});
tmp = data1.wnrep.psr{15}/tt(2);
tmp = conv(tmp,kern,'same')/50*tt(2);
plot(tmp*5,'k')
axis([0 450 0 1.5])
legend('data','model')
title('repeats')

save('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_008_011_chopped.mat','data1','data2')

%%
clear

load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_011_015.mat','data1','data2')

load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_015_019_chopped.mat','data1','data2')

load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_015_022.mat','data1','data2')

load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_011_015_chopped.mat','data1','data2')


load('/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/data_008_011_chopped.mat','data1','data2')

data='011';
starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
starun = load_params(starun,'verbose',1);
filepath = '/Users/alexth/Desktop/Light_adaptation/movie_GS/2015-03-09-2/d05-27-norefit/nd4_3_chopped/';
nds={'NDF4', 'NDF3'};


sr=120;
sig=50;
st=10000/sr*6.2/(60*sig);
time_list=-3.1:st:3.1;
kern=zeros(1,length(time_list));
for i=1:length(time_list)    
    kern(i)=250/sig*exp((1-time_list(i)^2)/2);
end


cell_cnt=1;
for cellID = starun.cell_ids    
    
    datarunID=find(starun.cell_ids==cellID);
    [my_name, my_type] = find_cell_type(starun, cellID);

    folder=my_name;
    if ~exist([filepath,folder],'dir')
        mkdir([filepath,folder]);
    end
    
    fig=figure('PaperPosition',[0 0 12 8],'PaperSize',[12 8]);
    set(fig,'color','white','position',[1 30 1710 1068]);
    
    
    subplot(5,1,1)    
    plot(1:1/120:30, data1.nm.gs{cell_cnt} * data1.nm.scale(cell_cnt),'r')
    hold on
    plot(1:1/120:30, data2.nm.gs{cell_cnt} * data2.nm.scale(cell_cnt),'k')
%     legend(nds)
    title(['GS, scaled by nonlinearity, ', num2str(data1.nm.scale(cell_cnt))...
        ,',  ', num2str(data2.nm.scale(cell_cnt))])
    axis tight
    tmp=get(gca,'Ylim');
    for i=1:1:30
        line([i,i],[tmp(1), tmp(2)],'color',[1 1 1]*0.5)
    end
    
        
    subplot(5,1,2)    
    a=std(data1.nm.gs{cell_cnt});
    b=std(data2.nm.gs{cell_cnt});
    plot(1:1/120:30, data1.nm.gs{cell_cnt} / a,'r')
    hold on
    plot(1:1/120:30, data2.nm.gs{cell_cnt} / b,'k')
%     legend(nds)
    title(['GS, scaled by std, ', num2str(a),',  ', num2str(b)])
    axis tight
    tmp=get(gca,'Ylim');
    for i=1:1:30
        line([i,i],[tmp(1), tmp(2)],'color',[1 1 1]*0.5)
    end
    
    
    subplot(5,1,3)
    
    tt = unique(data1.nm.psr{cell_cnt});
    tmp = data1.nm.psr{cell_cnt}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(1:1/120:30,tmp,'r')
    hold on
    tt = unique(data2.nm.psr{cell_cnt});
    tmp = data2.nm.psr{cell_cnt}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(1:1/120:30,tmp,'k')    
%     subplot(5,1,3)
%     plot(1:1/120:30,data1.nm.psr{cell_cnt},'r')
%     hold on
%     plot(1:1/120:30,data2.nm.psr{cell_cnt},'k')
%     legend(nds)
    title(['Predicted Spike Rate, r2 = ', num2str(data1.r2(cell_cnt,2)),...
        ',  ', num2str(data2.r2(cell_cnt,2))]);
    axis tight
    tmp=get(gca,'Ylim');
    for i=1:1:30
        line([i,i],[tmp(1), tmp(2)],'color',[1 1 1]*0.5)
    end
    
    subplot(5,1,4)
    tt = unique(data1.nm.asr{cell_cnt});
    tmp = data1.nm.asr{cell_cnt}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(1:1/120:30,tmp,'r')
    hold on
    tt = unique(data2.nm.asr{cell_cnt});
    tmp = data2.nm.asr{cell_cnt}/tt(2);    
    tmp = conv(tmp,kern,'same')/50*tt(2);
    plot(1:1/120:30,tmp,'k')
    legend(nds)
    title('Actual Spike Rate')
    axis tight
    tmp=get(gca,'Ylim');
    for i=1:1:30
        line([i,i],[tmp(1), tmp(2)],'color',[1 1 1]*0.5)
    end
    
    h=subplot(5,1,5);
    hold on
    spikes=data1.nm.raw_spikes{cell_cnt};
    trigs = data1.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        splitSpikes{j-5}=tmp*1000;
    end    

    time_trial=31000;
    fr=zeros(1,time_trial);
    splitRasters=[];
    for j=1:14
            splitRasters=[splitRasters splitSpikes{j}'+time_trial*(j-1)];             
            tmp=convolved(splitSpikes{j}',50,time_trial);
            fr=fr+tmp(((size(tmp,2)-time_trial)/2+1):end-((size(tmp,2)-time_trial)/2))/20;        
    end
    rasterplot_options(splitRasters,14,time_trial,0,'r', h)  
    axis([1000 30000 0 21])
    
    spikes=data2.nm.raw_spikes{cell_cnt};
    trigs = data2.nm.trigs;
    myTrigs=[0 find(diff(trigs)>0.9)'];
    splitSpikes=cell(14,1);
    for j=6:19
        tmp=spikes(spikes>=trigs(myTrigs(j)+1) & spikes<trigs(myTrigs(j+1)))...
            - trigs(myTrigs(j)+1);
        splitSpikes{j-5}=tmp*1000;
    end    
    fr=zeros(1,time_trial);
    splitRasters2=[];
    for j=1:14
            splitRasters2=[splitRasters2 splitSpikes{j}'+time_trial*(j-1)];             
            tmp=convolved(splitSpikes{j}',50,time_trial);
            fr=fr+tmp(((size(tmp,2)-time_trial)/2+1):end-((size(tmp,2)-time_trial)/2))/20;        
    end
    rasterplot_options(splitRasters2,14,time_trial,1.5*14+2,'k', h)
    axis([1000 30000 0  1.5*14*2+3])
    set(h,'xtick', (5:5:30)*1000, 'xticklabel',{'5', '10', '15', '20', '25', '30'}) 
    for k=1:30
        line([1,1]+k*1000,[0,1.5*19*2+3],'color',[1 1 1]*0.5,'linewidth',0.3)
    end
    
    title(['Cell ', int2str(cellID), '   (',int2str(datarunID),'),   ', my_name, ',   trials 6-19'])
    
    print(fig,'-dpdf',sprintf('%s%s%s.pdf',[filepath,folder,'/'],['cell_',int2str(cellID)]));
    close(fig)
    
    cell_cnt=cell_cnt+1;
end



figure
plot(data1.gs{1})
hold on
plot(data2.gs{1})



figure
a=data1.sta{15}
figure
plot(mean(a),'linewidth',2)
hold on
a=data2.sta{15}
plot(mean(a),'linewidth',2)
legend('NDF3', 'NDF2')

a = data1.locs{15,2}
a = data2.locs{15,2}

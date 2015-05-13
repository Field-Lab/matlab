function [total_STA, total_GS, total_ASR, total_bins, total_NL] = fit_NSEM_ndf(data)

downsize=1;
offset = 0;
int_method='spline';
sta_length=15;
sigstix_threshold = 4;
movie_duration = 31000; % ms,single repetition
nm_repeats = 20;


mvpath='/Volumes/Data/stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/';
my_movie=zeros(320,320,3600);
cnt=1;
for i=1:30
    load([mvpath, 'movie_chunk_', int2str(i), '.mat']);
    movie=movie(81:240,:,:);
    movie=imresize(movie,2,'nearest');
    my_movie(:,:,cnt:cnt+119)=movie;
    cnt=cnt+120;
end
figure
hist(my_movie(:),50)
my_movie=my_movie/255-0.5;
clear mvpath movie

for i=1
    switch data.name
        case '004'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            frames=8;
            ndf=4;
            nmdata='009';
        case '008'
            wn_movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            frames=8;
            ndf=4;
            nmdata='009';
        case '011'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            frames=8;
            ndf=3;
            nmdata='012';
        case '014'
            wn_movie_name = 'BW-8-8-0.48-11111-40x40.xml';
            frames=8;
            ndf=3;
            nmdata='012';
        case '015'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            frames=6;
            ndf=2;
            nmdata='016';
        case '018'
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            frames=6;
            ndf=2;
            nmdata='016';
        case '019'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            frames=6;
            ndf=1;
            nmdata='020';
        case '022'
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            frames=6;
            ndf=1;
            nmdata='020';
        case '023'
            wn_movie_name = 'BW-16-4-0.48-11111-20x20.xml';
            frames=4;
            ndf=0;
            nmdata='024';
        case '027'
            wn_movie_name = 'BW-4-2-0.48-11111-80x80.xml';
            frames=2;
            ndf=0;
            nmdata='024';
    end
  
    nmrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',nmdata,'-from-d05-d27/data',nmdata,'-from-d05-d27']);
    nmrun = load_params(nmrun,'verbose',1);
    nmrun = load_neurons(nmrun);
end
clear nmdata

space_size = data.aux(3)*data.aux(4);

scaled_movie=zeros(space_size, 3600);
stix_size=320/sqrt(space_size);
cnt=1;
for i=1:stix_size:320
    for j=1:stix_size:320
        tmp=reshape(my_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end
clear stix_size


nm_triggers=nmrun.triggers;
movie_refresh=data.aux(1)/frames;%1/120*1000; %ms

predicted_spike_rate_acc = [];
spike_rate_acc = [];
gen_sig_acc = [];
r2=zeros(length(nmrun.cell_ids),2);

cell_cnt=1;
for cellID = nmrun.cell_ids
    
    datarunID=find(nmrun.cell_ids==cellID);    
    
    %****** PART 1 fit WN ******   
    
    gen_sig = data.gs{cell_cnt}; * data.scale(cell_cnt);
    
    spike_rate = data.asr{cell_cnt}; 
    
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
    
    
    %****** PART 2 fit NSEM using WN STA ******
    
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
    sta = data.sta{cell_cnt} * data.scale(cell_cnt);
    if downsize==1
        tt=[];
        for i=1:size(sta,1)
            tt=[tt; interp1(1:8:8*size(sta,2), sta(i,:),1:8*size(sta,2),int_method)]; % spline, linear
        end
        sta=tt;
    end
    % convolve raw input with linear filter
    inds = data.locs{cell_cnt,2};
    my_inputs=scaled_movie(inds, :);
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),sta(i,:),'valid');
    end
    
    gen_sig=sum(my_filtered_inputs);
    
    gen_sig_acc = [gen_sig_acc; gen_sig];
    
    acc_spike_rate = acc_spike_rate(size(sta,2):end);
    
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
    
    predicted_spike_rate_acc = [predicted_spike_rate_acc predicted_rate];
    spike_rate_acc = [spike_rate_acc acc_spike_rate];
    
    cell_cnt=cell_cnt+1;
end

figure
plot(r2)


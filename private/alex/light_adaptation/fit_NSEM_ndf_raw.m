function data = fit_NSEM_ndf_raw(data, even_bin)

data.aux.int_method='spline';
data.aux.movie_duration = 31000; % ms,single repetition
data.aux.nm_repeats = 20;


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
my_movie=my_movie/255-0.5;
my_movie=my_movie-mean(my_movie(:));
clear mvpath movie

nmrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data.aux.nmdata,'-from-d05-d27/data',data.aux.nmdata,'-from-d05-d27']);
nmrun = load_params(nmrun,'verbose',1);
nmrun = load_neurons(nmrun);

clear nmdata

space_size = data.aux.height*data.aux.width;

scaled_movie=zeros(space_size, 3600);
stix_size=320/sqrt(space_size);
cnt=1;
for i=1:stix_size:320
    for j=1:stix_size:320
%         tmp=reshape(my_movie(j:j+stix_size-1,i:i+stix_size-1, :), stix_size*stix_size, 3600);
% this takes care of rotation!!
        tmp=reshape(my_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end
clear stix_size


data.nm.trigs=nmrun.triggers;
movie_refresh=data.aux.refresh/data.aux.frames;%1/120*1000; %ms

r2=zeros(length(nmrun.cell_ids),2);


cell_cnt=1;
for cellID = nmrun.cell_ids
    
    datarunID=find(nmrun.cell_ids==cellID);    
    
    %****** PART 1 fit WN ******   
    
%     gen_sig = data.gs{cell_cnt} * data.scale(cell_cnt);
    gen_sig = data.gs{cell_cnt};
    spike_rate = data.asr{cell_cnt};
    
    % nonlinearity    
    if even_bin
        n_gs=floor(length(gen_sig)/100);
        tmp=sort(gen_sig);
        my_bins=[min(gen_sig) tmp(n_gs:n_gs:n_gs*100)];
        my_bins(end)=max(gen_sig)+1;
    else
        n_gs = (max(gen_sig)-min(gen_sig))/100;        
        my_bins = min(gen_sig):n_gs:(max(gen_sig)+0.01);
    end
    my_nl=zeros(size(my_bins,2)-1,1);  
    my_std=my_nl;
    for i=1:100
        tmp=find(gen_sig>=my_bins(i) & gen_sig<my_bins(i+1)); % find instances when GS had certain value
        my_std(i)=std(spike_rate(tmp));
        my_nl(i)=mean(spike_rate(tmp)); % find mean firing rate after this certain value
    end

    % calculate predicted FR
    cnt=1;
    predicted_rate=zeros(size(spike_rate));
    for i=1:100
        tmp=find(gen_sig>=my_bins(i) & gen_sig<my_bins(i+1));
        predicted_rate(tmp)=my_nl(i);
        cnt=cnt+1;
    end
    sse = sum((predicted_rate-spike_rate).^2);
    sst = sum((spike_rate-mean(spike_rate)).^2);
    r2(cell_cnt,1) = 1 - sse/sst;
    
    
    %****** PART 2 fit NSEM using WN STA ******
    
    data.nm.raw_spikes{cell_cnt} = nmrun.spikes{datarunID};
    
    spikes=ceil((nmrun.spikes{datarunID}-nmrun.triggers(1))*1000/movie_refresh); % spikes in frames
    spikes(spikes<1)=[];    
    spike_rate=zeros(data.aux.movie_duration*data.aux.nm_repeats,1);    
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end    
    
    my_begs=round([0; nmrun.triggers(find(diff(nmrun.triggers)>0.9)+1)]*1000/movie_refresh);
    acc_spike_rate=0;
    for i=6:length(my_begs)-1
        acc_spike_rate = acc_spike_rate + spike_rate(my_begs(i)+1:my_begs(i)+3600);
    end
    acc_spike_rate=acc_spike_rate/(length(6:length(my_begs)-1));
    
    % interpolate linear filter to match movie frames
%     sta = data.sta{cell_cnt} * data.scale(cell_cnt);
    sta = data.sta{cell_cnt};
    if data.downsize==1
        tt=[];
        for i=1:size(sta,1)
            tt=[tt; interp1(1:data.aux.frames:data.aux.frames*size(sta,2), sta(i,:), 1:data.aux.frames*size(sta,2),data.aux.int_method)]; % spline, linear
        end
        sta=tt;
    end
    if size(sta,2)<120
        sta(:,end:end+(120-size(sta,2)))=0;
    end
    % convolve raw input with linear filter
    inds = data.locs{cell_cnt,2};
    my_inputs=scaled_movie(inds, :);
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),sta(i,:),'valid');
    end
    
    gen_sig=sum(my_filtered_inputs);    
    data.nm.gs{cell_cnt} = gen_sig;
    
    acc_spike_rate = acc_spike_rate(size(sta,2):end);
    data.nm.asr{cell_cnt} = acc_spike_rate;
    
%     tt = unique(acc_spike_rate);
% tmp = acc_spike_rate/tt(2);
% tmp = conv(tmp,kern,'same')/50*tt(2);
% figure
% plot(tmp,'k')

    
    
    % nonlinearity
    if even_bin
        n_gs=floor(length(gen_sig)/100);
        tmp=sort(gen_sig);
        my_bins=[min(gen_sig) tmp(n_gs:n_gs:n_gs*100)];
        my_bins(end)=max(gen_sig)+1;
    else
        n_gs = (max(gen_sig)-min(gen_sig))/100;        
        my_bins = min(gen_sig):n_gs:(max(gen_sig)+0.01);
    end
    my_nl=zeros(size(my_bins,2)-1,1);
    my_std=my_nl;
    for i=1:100
        tmp=find(gen_sig>=my_bins(i) & gen_sig<my_bins(i+1)); % find instances when GS had certain value
        my_std(i)=std(acc_spike_rate(tmp));
        my_nl(i)=mean(acc_spike_rate(tmp)); % find mean firing rate after this certain value
    end
    
    
    data.nm.bins{cell_cnt} = my_bins;  
    data.nm.nl{cell_cnt} = my_nl;
    
    % calculate predicted FR    
    cnt=1;
    predicted_rate=zeros(size(acc_spike_rate));
    for i=1:100
        tmp=find(gen_sig>=my_bins(i) & gen_sig<my_bins(i+1));
        predicted_rate(tmp)=my_nl(i);
        cnt=cnt+1;
    end
    sse = sum((predicted_rate-acc_spike_rate).^2);
    sst = sum((acc_spike_rate-mean(acc_spike_rate)).^2);
    r2(cell_cnt,2) = 1 - sse/sst;
    
    data.nm.psr{cell_cnt} = predicted_rate;

    cell_cnt=cell_cnt+1;
end

data.r2=r2;

function data = fit_WNrep_ndf_raw(data, even_bin)

data.aux.int_method='spline';
data.wnrep.n_of_repeats = 20; 
data.wnrep.duration = 30000; % ms

wnrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data.aux.repeats,'-from-d05-d27/data',data.aux.repeats,'-from-d05-d27']);
wnrun = load_params(wnrun,'verbose',1);
wnrun = load_neurons(wnrun);

[inputs, refresh, duration] = get_wn_movie_ath(wnrun, data.aux.wn_repeats_name);
% 3588 frames 
wn_frames=round(refresh/(1000/120));
inputs_per_rep=floor(3588/wn_frames);

inputs = inputs(:,1:inputs_per_rep);
inputs = reshape(inputs, sqrt(size(inputs,1)),sqrt(size(inputs,1)),inputs_per_rep);
% figure
% imagesc(inputs(:,:,1))
inputs = imresize(inputs, 320/size(inputs,1), 'nearest');
% 
% a=reshape(scaled_movie, 40, 40, 448);
% figure
% imagesc(a(:,:,1))

space_size = data.aux.height*data.aux.width;
scaled_movie=zeros(space_size, inputs_per_rep);
stix_size=320/sqrt(space_size);
cnt=1;
for i=1:stix_size:320
    for j=1:stix_size:320
        tmp=reshape(inputs(j:j+stix_size-1, i:i+stix_size-1,:), stix_size*stix_size, inputs_per_rep);
        scaled_movie(cnt,:)=mean(tmp);
        cnt=cnt+1;
    end
end
clear stix_size

data.wnrep.trigs=wnrun.triggers;
data.wnrep.refresh=refresh;%ms

r2=zeros(length(wnrun.cell_ids),2);


% cell_cnt=15; cellID=529
cell_cnt=1;
for cellID = wnrun.cell_ids
    
    datarunID=find(wnrun.cell_ids==cellID);    
    
    %****** PART 1 fit WN ******   
    
%     gen_sig = data.gs{cell_cnt} * data.scale(cell_cnt);
    gen_sig = data.gs{cell_cnt};
    spike_rate = data.asr{cell_cnt};
    my_nl = data.nl{cell_cnt};
    my_bins = data.bins{cell_cnt};

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
    
    data.psr{cell_cnt} = predicted_rate;
    
    
    %****** PART 2 fit NSEM using WN STA ******
    
    data.wnrep.raw_spikes{cell_cnt} = wnrun.spikes{datarunID};

    spikes=ceil((wnrun.spikes{datarunID}-wnrun.triggers(1))*1000/refresh); % spikes in frames
    spikes(spikes<1)=[];
    max_duration = data.wnrep.n_of_repeats * ceil(data.wnrep.duration/refresh);
    spikes(spikes>max_duration) = [];
    spike_rate=zeros(max_duration,1);    
    while ~isempty(spikes)
        [c, ia, ic] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        spikes(ia)=[];
    end    
    
    my_begs=round([0; wnrun.triggers(find(diff(wnrun.triggers)>0.84)+1)]*1000/refresh);
    acc_spike_rate=0;
    for i=6:length(my_begs)-1
        acc_spike_rate = acc_spike_rate + spike_rate(my_begs(i)+1:my_begs(i)+inputs_per_rep);
    end
    acc_spike_rate=acc_spike_rate/(length(6:length(my_begs)-1));
    
    % interpolate linear filter to match repeats frames
%     sta = data.sta{cell_cnt} * data.scale(cell_cnt);
    sta = data.sta{cell_cnt};
    if data.downsize==1
        tt=[];
        for i=1:size(sta,1)
            tt=[tt; interp1(1:data.aux.frames:data.aux.frames*size(sta,2), sta(i,:), 1:wn_frames:wn_frames*size(sta,2),data.aux.int_method)]; % spline, linear
        end
        sta=tt;
    end
    if size(sta,2)<120/data.aux.frames
        sta(:,end:end+(120/data.aux.frames-size(sta,2)))=0;
    end
    % convolve raw input with linear filter
    inds = data.locs{cell_cnt,2};
    my_inputs=scaled_movie(inds, :);
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),sta(i,:),'valid');
    end
    
    gen_sig=sum(my_filtered_inputs);    
    data.wnrep.gs{cell_cnt} = gen_sig;
    
    acc_spike_rate = acc_spike_rate(size(sta,2):end);
    data.wnrep.asr{cell_cnt} = acc_spike_rate;
    
    
   
    % calculate nonlinearity for comparison purposes
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
    for i=1:100
        tmp=find(gen_sig>=my_bins(i) & gen_sig<my_bins(i+1)); % find instances when GS had certain value
        my_nl(i)=mean(spike_rate(tmp)); % find mean firing rate after this certain value
    end
    
    data.wnrep.bins{cell_cnt} = my_bins;
    data.wnrep.nl{cell_cnt} = my_nl;
    
    
    % take nonlinearity from long run
    my_nl = data.nl{cell_cnt};
    my_bins = data.bins{cell_cnt};
    
    if min(gen_sig)<my_bins(1) || max(gen_sig)>my_bins(end)
        disp(['cell ', int2str(cellID), '  has GS out of NL range!'])
    end
    
   
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
    
    data.wnrep.psr{cell_cnt} = predicted_rate;

    cell_cnt=cell_cnt+1;
end

data.wnrep.r2=r2;

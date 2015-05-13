function data_st = fit_LN_ndf(data, mean_rgb, even_bin)

data_st.name = data;
data_st.downsize=1;
data_st.offset = 0;
data_st.sta_length=15;
data_st.sigstix_threshold = 4;

for i=1
    switch data
        case '004'
            data_st.aux.wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            data_st.aux.frames=8;
            data_st.aux.ndf=4;
            data_st.aux.nmdata='009';
            data_st.aux.repeats='010';
            data_st.aux.wn_repeats_name = 'BW-10-8-0.48-11111-32x32.xml';
        case '008'
            data_st.aux.wn_movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            data_st.aux.frames=8;
            data_st.aux.ndf=4;
            data_st.aux.nmdata='009';
            data_st.aux.repeats='010';
            data_st.aux.wn_repeats_name = 'BW-10-8-0.48-11111-32x32.xml';
        case '011'
            data_st.aux.wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            data_st.aux.frames=8;
            data_st.aux.ndf=3;
            data_st.aux.nmdata='012';
            data_st.aux.repeats='013';
            data_st.aux.wn_repeats_name = 'BW-8-8-0.48-11111-40x40.xml';
        case '014'
            data_st.aux.wn_movie_name = 'BW-8-8-0.48-11111-40x40.xml';
            data_st.aux.frames=8;
            data_st.aux.ndf=3;
            data_st.aux.nmdata='012';
            data_st.aux.repeats='013';
            data_st.aux.wn_repeats_name = 'BW-8-8-0.48-11111-40x40.xml';
        case '015'
            data_st.aux.wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            data_st.aux.frames=6;
            data_st.aux.ndf=2;
            data_st.aux.nmdata='016';
            data_st.aux.repeats='017';
            data_st.aux.wn_repeats_name = 'BW-8-6-0.48-11111-40x40.xml';
        case '018'
            data_st.aux.wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            data_st.aux.frames=6;
            data_st.aux.ndf=2;
            data_st.aux.nmdata='016';
            data_st.aux.repeats='017';
            data_st.aux.wn_repeats_name = 'BW-8-6-0.48-11111-40x40.xml';
        case '019'
            data_st.aux.wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            data_st.aux.frames=6;
            data_st.aux.ndf=1;
            data_st.aux.nmdata='020';
            data_st.aux.repeats='021';
            data_st.aux.wn_repeats_name = 'BW-8-6-0.48-11111-40x40.xml';
        case '022'
            data_st.aux.wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            data_st.aux.frames=6;
            data_st.aux.ndf=1;
            data_st.aux.nmdata='020';
            data_st.aux.repeats='021';
            data_st.aux.wn_repeats_name = 'BW-8-6-0.48-11111-40x40.xml';
        case '023'
            data_st.aux.wn_movie_name = 'BW-16-4-0.48-11111-20x20.xml';
            data_st.aux.frames=4;
            data_st.aux.ndf=0;
            data_st.aux.nmdata='024';            
            data_st.aux.repeats='025';
            data_st.aux.wn_repeats_name = 'BW-4-2-0.48-11111-80x80.xml';
        case '027'
            data_st.aux.wn_movie_name = 'BW-4-2-0.48-11111-80x80.xml';
            data_st.aux.frames=2;
            data_st.aux.ndf=0;
            data_st.aux.nmdata='024';            
            data_st.aux.repeats='025';
            data_st.aux.wn_repeats_name = 'BW-4-2-0.48-11111-80x80.xml';
    end
    starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
    starun = load_params(starun,'verbose',1);
    starun = set_polarities(starun);
    starun = load_neurons(starun);
    starun = load_sta(starun,'load_sta','all','keep_java_sta',true);
    
end

[inputs, refresh, duration] = get_wn_movie_ath(starun, data_st.aux.wn_movie_name);

data_st.aux.refresh = refresh;
data_st.aux.duration = duration;


% figure
tc=zeros(length(starun.cell_types),length(starun.vision.timecourses(1).g));
for j=1:size(tc,1)
    cnt=1;
   % subplot(3,4,j)
    for i=starun.cell_types{j}.cell_ids
        tmp=find(starun.cell_ids==i);
        if ~isempty(starun.vision.timecourses(tmp).g)
            if mean_rgb % for rod-dominated RGB runs, take average of green and blue channels
                tc(j,:)=tc(j,:)+mean([starun.vision.timecourses(tmp).b'; starun.vision.timecourses(tmp).g']);
            elseif  max(starun.vision.timecourses(tmp).g)<max(starun.vision.timecourses(tmp).b) % for blue cells, RGB run                
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).b';
            else % BW run, all cells, or RGB run, green cells
                %if sum(starun.vision.timecourses(tmp).g == starun.vision.timecourses(tmp).b)
                tc(j,:)=tc(j,:)+starun.vision.timecourses(tmp).g';
            end
%              plot(starun.vision.timecourses(tmp).g','g')
%              hold on
%              plot(starun.vision.timecourses(tmp).b','b')
            cnt=cnt+1;
        end
    end
    tc(j,:)=tc(j,:)'/cnt;
  
end

cell_cnt=1;
for cellID = starun.cell_ids
    datarunID=find(starun.cell_ids==cellID);
    
    [~, my_type] = find_cell_type(starun, cellID);
    
    %****** PART 1 calculate linear filter from STA ******
    % sta of a chosen cell
    sta = double(squeeze(starun.stas.stas{datarunID}));
        
    if ndims(sta)==4
        if mean_rgb
            sta=squeeze(mean(sta,3));
        elseif max(starun.vision.timecourses(datarunID).g)<max(starun.vision.timecourses(datarunID).b)
            sta=squeeze(sta(:,:,3,:)); % take blue channel
        else
            sta=squeeze(sta(:,:,2,:)); % take green channel
        end
    end
    
    % get weights by convolving sta of a chosen cell with mean time course
    temp_sta = reshape(sta, size(sta,1)*size(sta,2),size(sta,3));
    weights=sum(repmat(tc(my_type,:),size(sta,1)*size(sta,2),1) .* temp_sta, 2);
    weights = weights/max(weights);
    
    % threshold weights
    weight_threshold=robust_std(weights)*data_st.sigstix_threshold;
    inds=find(weights>weight_threshold);
    cntstx = 0;
    while length(inds)<2
        cntstx=cntstx+0.5;
        weight_threshold=robust_std(weights)*(data_st.sigstix_threshold-cntstx);
        inds=find(weights>weight_threshold);
    end
    
    data_st.locs{cell_cnt,1} = data_st.sigstix_threshold-cntstx;
    data_st.locs{cell_cnt,2} = inds;
    
    my_inputs = inputs(inds, :); % to save memory and time
%     inp= zeros(length(inds),size(inputs,2)*data_st.downsize);
%     % down to a frame: instead of refresh of 66.6192, use refresh/8 = 8.3274
%     for i=1:data_st.downsize
%         inp(:,i:data_st.downsize:end) = my_inputs;
%     end
%     my_inputs=inp;
        
    spikes=ceil((starun.spikes{datarunID}-starun.triggers(1))*1000/(refresh/data_st.downsize)); % spikes in frames
%     spikes=floor((starun.spikes{datarunID}-starun.triggers(1))*1000/(refresh/data_st.downsize)); % spikes in frames
%     data_st.offset = 2
    spikes(spikes<data_st.sta_length-data_st.offset)=[];
    
    spike_rate=zeros(duration*data_st.downsize,1);
    my_sta=zeros(length(inds),data_st.sta_length);
    
    while ~isempty(spikes)
        [~, ia, ~] = unique(spikes);
        spike_rate(spikes(ia))=spike_rate(spikes(ia))+1;
        for j=1:data_st.sta_length
            my_sta(:,data_st.sta_length-j+1)=my_sta(:,data_st.sta_length-j+1)+sum(my_inputs(:,spikes(ia)-data_st.sta_length+j+data_st.offset),2);
        end
        spikes(ia)=[];
    end
    my_sta=my_sta/sum(spike_rate);
    
     
    data_st.sta{cell_cnt} = my_sta;
    
    % convolve raw input with linear filter
    my_filtered_inputs=zeros(size(my_inputs,1),size(my_inputs,2)-size(my_sta,2)+1);
    for i=1:length(inds)
        my_filtered_inputs(i,:)=conv(my_inputs(i,:),my_sta(i,:),'valid');
    end
         
    gen_sig=sum(my_filtered_inputs);    
    data_st.gs{cell_cnt} = gen_sig; 
    
    spike_rate = spike_rate(data_st.sta_length:end);
    data_st.asr{cell_cnt} = spike_rate;
    
    
    % calculate nonlinearity    
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

    data_st.nl{cell_cnt} = my_nl;
    data_st.bins{cell_cnt} = my_bins;
    cell_cnt=cell_cnt+1;
end

sta = double(squeeze(starun.stas.stas{1}));
data_st.aux.height = size(sta,1);
data_st.aux.width = size(sta,2);
    

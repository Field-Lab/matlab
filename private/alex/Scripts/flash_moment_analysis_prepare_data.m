clear
date='20130220'

path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);

trials=82;
res=zeros(46+28,9000,9);
for unit=1:length(flashes) 
    
   % FULL FIELD FLASH
    load([path2data,flashes(unit).name])
    
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,4501);
    black_flash_tmp=zeros(trials,4501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end    
   
    cnt=1;
    for i=[1 11:9:81]
        white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    
    
    res(unit,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
    
end   
   
date='20130220_1'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);

trials=81;

for unit=1:length(flashes) 
    
   % FULL FIELD FLASH
    load([path2data,flashes(unit).name])
    
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,4501);
    black_flash_tmp=zeros(trials,4501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end    
   
    cnt=1;
    for i=1:9:81
        white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    
    
    res(unit+46,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
    
end   

date='20130225'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);

trials=81;

for unit=1:length(flashes) 
    
   % FULL FIELD FLASH
    load([path2data,flashes(unit).name])
    
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,4501);
    black_flash_tmp=zeros(trials,4501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end    
   
    cnt=1;
    for i=1:9:81
        white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    
    
    res(unit+74,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
    
end   

date='20130226'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
flashes=dir([path2data,'*FFFlash*']);

trials=81;

for unit=1:length(flashes) 
    
   % FULL FIELD FLASH
    load([path2data,flashes(unit).name])
    
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,4501);
    black_flash_tmp=zeros(trials,4501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:4501)=white_flash_tmp(trial,1:4501)+conv(flips(st)-500:flips(st)+4000);
            black_flash_tmp(trial,1:4501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+4000);
        end
    end    
   
    cnt=1;
    for i=1:9:81
        white_flash(cnt,1:4501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:4501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;
    
    
    res(unit+121,:,:)=[white_flash(:,1:4500)'; black_flash(:,1:4500)'];
    
end   


timing_flash=zeros(1,9000);
timing_flash(501:1984+501)=1;
timing_flash(5001:1984+5001)=-1;
flash=res;








clear res
date='20130220'
        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);

trials=82;
res=zeros(46+28,5000,9);
for unit=1:length(moments) 
    
    load([path2data,moments(unit).name])
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,2501);
    black_flash_tmp=zeros(trials,2501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
            black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
        end
    end
    
   
    cnt=1;
    for i=[1 11:9:81]
        white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;    
    
    res(unit,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
    
end   
   
date='20130220_1'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);

trials=81;

for unit=1:length(moments) 
    
    load([path2data,moments(unit).name])
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,2501);
    black_flash_tmp=zeros(trials,2501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
            black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
        end
    end    
    
    cnt=1;
    for i=1:9:80
        white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;    
    
    res(unit+46,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
    
end   
 
date='20130225'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);

trials=81;

for unit=1:length(moments) 
    
    load([path2data,moments(unit).name])
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,2501);
    black_flash_tmp=zeros(trials,2501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
            black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
        end
    end    
    
    cnt=1;
    for i=1:9:80
        white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;    
    
    res(unit+74,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
    
end   
 
date='20130226'        
path2data=['S:\data\alexandra\MEA_data\',date,'\easy_formatted_units\'];
moments=dir([path2data,'*Moment*']);

trials=81;

for unit=1:length(moments) 
    
    load([path2data,moments(unit).name])
    clear white_flash black_flash
    white_flash_tmp=zeros(trials,2501);
    black_flash_tmp=zeros(trials,2501);
    for trial=1:trials
        spikes=cell2mat(spike_info.spike_times(trial));
        flips=cell2mat(spike_info.flip_times(trial));
        flips=flips(:,1);
        conv=convolved(spikes,40,flips(end));
        conv=conv(121:end-120);
        for st=1:4:20
            white_flash_tmp(trial,1:2501)=white_flash_tmp(trial,:)+conv(flips(st)-500:flips(st)+2000);
            black_flash_tmp(trial,1:2501)=black_flash_tmp(trial,:)+conv(flips(st+2)-500:flips(st+2)+2000);
        end
    end    
    
    cnt=1;
    for i=1:9:80
        white_flash(cnt,1:2501)=mean(white_flash_tmp(i:i+8,:));
        black_flash(cnt,1:2501)=mean(black_flash_tmp(i:i+8,:));
        cnt=cnt+1;
    end
    white_flash=white_flash/5;
    black_flash=black_flash/5;    
    
    res(unit+121,:,:)=[white_flash(:,1:2500)'; black_flash(:,1:2500)'];
    
end   

timing_moments=zeros(1,5000);
timing_moments(501:17)=1;
timing_moments(2501:17+2501)=-1;
moments=res;







date='20130220'

filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);
meanFilter=zeros(46+28,500,9);
for unit=1:length(filters)
    load([filterspath,filters(unit).name])    
    LF=HighFilters(:,[1:9 11:82])';
    cnt=1;
    for i=1:9:80
        meanFilter(unit,:,cnt)=mean(LF(i:i+8,:));
        cnt=cnt+1;
    end
end

date='20130220_1'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

for unit=1:length(filters)
    load([filterspath,filters(unit).name])    
    LF=HighFilters';
    cnt=1;
    for i=1:9:80
        meanFilter(unit+46,:,cnt)=mean(LF(i:i+8,:));
        cnt=cnt+1;
    end
end

date='20130225'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

for unit=1:length(filters)
    load([filterspath,filters(unit).name])    
    LF=HighFilters';
    cnt=1;
    for i=1:9:80
        meanFilter(unit+74,:,cnt)=mean(LF(i:i+8,:));
        cnt=cnt+1;
    end
end

date='20130226'
filterspath=['S:\data\alexandra\MEA_data\',date,'\FFFlicker_LF\'];
filters=dir([filterspath,'*FFFlicker_linear_filter.mat']);

for unit=1:length(filters)
    load([filterspath,filters(unit).name])    
    LF=HighFilters';
    cnt=1;
    for i=1:9:80
        meanFilter(unit+121,:,cnt)=mean(LF(i:i+8,:));
        cnt=cnt+1;
    end
end



for unit=1:size(meanFilter,1)
    if mean(meanFilter(unit,130:140,4),2)>0
        a=reshape(max(meanFilter(unit,:,:)),9,1);
    else
        a=reshape(min(meanFilter(unit,:,:)),9,1);
    end
    tmp=reshape(meanFilter(unit,:,:),500,9)'./repmat(abs(a),1,500);
    meanFilter(unit,:,:)=tmp';
end

snr=std(meanFilter(:,60:300,4),0,2)./std(meanFilter(:,1:20,4),0,2);
onOff=mean(meanFilter(:,130:140,4),2);



save('S:\data\alexandra\MEA_data\flash_analysis\flash_moment.mat','flash','timing_flash','moments','timing_moments','meanFilter','snr','onOff')


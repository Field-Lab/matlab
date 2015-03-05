clear
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

for i=1:9
    subplot(3,3,i)
    plot(meanFilter(snr>30&onOff>0,1:300,i)')
end



for i=1:9
    subplot(3,3,i)
    plot(meanFilter(snr>35,1:300,i)')
    axis([0 300 -1 1])
end




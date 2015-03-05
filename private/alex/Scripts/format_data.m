function pulled=format_data(date,typ,experiment_code)

path2data=['S:\user\alexandra\MEA_data\',date,'\'];
fits=dir([path2data,'fit_res_HC\*fit_res_HC.mat']);
pulled=nan(length(fits),96,13);
for unit=1:length(fits)
    load([path2data,'fit_res_HC\',fits(unit).name]);
    name=fits(unit).name(1:end-15);
    unitNumber=str2num(name(end-3:end));
    channelNumber=name(13:14);
    if channelNumber(2)=='_'
        channelNumber=str2num(channelNumber(1));
    else
        channelNumber=str2num(channelNumber);
    end
    quality=name(end-12:end-10);
    if quality(1)=='_'
        quality=str2num(quality(2:3));
    else
        quality=str2num(quality);
    end
    
    cellType=sum(sign(common_res_fit(1,:)));
    if cellType>0
        cellType=1;%on
    else
        cellType=-1;%off
    end
    load([path2data,'easy_formatted_units\',name,'___spike_info_FFFlicker'])
    if typ==12
        cnt=1;
        clear nds
        for i=1:length(spike_info.name_info)
            if ~isempty(regexp(spike_info.name_info{i},'HCseq'))
                tmp=regexp(spike_info.name_info{i},'ND');
                tmp=spike_info.name_info{i}(tmp+2);
                if length(tmp)>1
                    nds(cnt)=str2num(tmp(1))+str2num(tmp(2));
                else
                    nds(cnt)=str2num(tmp);
                end
                if cnt>1&&nds(cnt)>nds(cnt-1)
                    nds(cnt)=[];
                    break
                end
                cnt=cnt+1;
            end
        end 
        startPoint=(8-nds(1))*12+1;
        timeSt=startPoint:startPoint+length(unique(nds))*12-1;
        totTime=60;
    elseif typ==6
        cnt=1;
        clear nds
        for i=1:length(spike_info.name_info)
            if ~isempty(regexp(spike_info.name_info{i},'ffflicker'))
                tmp=regexp(spike_info.name_info{i},'ND');
                tmp=spike_info.name_info{i}(tmp+2);
                if length(tmp)>1
                    nds(cnt)=str2num(tmp(1))+str2num(tmp(2));
                else
                    nds(cnt)=str2num(tmp);
                end
                if cnt>1&&nds(cnt)>nds(cnt-1)
                    nds(cnt)=[];
                    break
                end
                cnt=cnt+1;
            end
        end 
        startPoint=(8-nds(1))*12+1;
        timeSt=startPoint:2:startPoint+length(unique(nds))*12-1;        
        totTime=30;
    else
        cnt=1;
        clear nds
        for i=1:length(spike_info.name_info)
            if ~isempty(regexp(spike_info.name_info{i},'ffflicker'))
                tmp=regexp(spike_info.name_info{i},'ND');
                tmp=spike_info.name_info{i}(tmp+2);
                if length(tmp)>1
                    nds(cnt)=str2num(tmp(1))+str2num(tmp(2));
                else
                    nds(cnt)=str2num(tmp);
                end
                if cnt>1&&nds(cnt)>nds(cnt-1)
                    nds(cnt)=[];
                    break
                end
                cnt=cnt+1;
            end
        end 
        startPoint=(8-nds(1))*12+1;
        timeSt=sort([startPoint+1:12:startPoint+length(unique(nds))*12-1 ...
        startPoint+3:12:startPoint+length(unique(nds))*12-1 ... 
        startPoint+9:12:startPoint+length(unique(nds))*12-1 ...
        startPoint+11:12:startPoint+length(unique(nds))*12-1]);
        totTime=30;
    end
    load([path2data,'FFFlicker_LF\',name,'_FFFlicker_linear_filter.mat'])    

    pulled(unit,timeSt,1:3)=common_res_fit(1:3,1:length(timeSt))';
    pulled(unit,timeSt,4)=cellType;
    pulled(unit,timeSt,5)=abs(common_res_fit(1,1:length(timeSt)))./std(HighFilters(1:50,1:length(timeSt)));
    pulled(unit,timeSt,6)=HighSpikeCount(1:length(timeSt))/totTime;
    pulled(unit,timeSt,7)=totTime;
    pulled(unit,timeSt,8)=nds;
    pulled(unit,timeSt,9)=typ;
    pulled(unit,timeSt,10)=unitNumber;
    pulled(unit,timeSt,11)=channelNumber;
    pulled(unit,timeSt,12)=quality;
    pulled(unit,timeSt,13)=experiment_code;
end         

    



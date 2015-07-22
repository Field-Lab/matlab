spikes=dir('/Volumes/Analysis/nora/NSEM/BlockedSpikes/2012-08-09-3/WN_mapPRJ/');
n=length(spikes);
count=0;
for i=1:n
    try
        if strcmp(spikes(i).name(20:22),'Par')
            count=count+1;
            cid_cell{count}=str2num(spikes(i).name(24:(end-4)));
        elseif strcmp(spikes(i).name(19:21),'Par')
            count=count+1;
            cid_cell{count}=str2num(spikes(i).name(23:(end-4)));
        end
    catch
    end
end

WN_cid=cell2mat(cid_cell);
clear i n spikes cid_cell cid count

spikes=dir('/Volumes/Analysis/nora/NSEM/BlockedSpikes/2012-08-09-3/NSEM_mapPRJ/');
n=length(spikes);
count=0;
for i=1:n
    try
        if strcmp(spikes(i).name(20:22),'Par')
            count=count+1;
            cid_cell{count}=str2num(spikes(i).name(24:(end-4)));
        elseif strcmp(spikes(i).name(19:21),'Par')
            count=count+1;
            cid_cell{count}=str2num(spikes(i).name(23:(end-4)));
        end
    catch
    end
end

NSEM_cid=cell2mat(cid_cell);

%%
for i=1:69
    a=WN_cid(WN_cid==NSEM_cid(i));
    if isempty(a) 
        disp(i); 
    end 
end
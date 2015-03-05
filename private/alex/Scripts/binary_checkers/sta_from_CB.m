function res=sta_from_CB(codeWord, pathway, sta)

load('/Users/alexth/Desktop/Scripts/single_cone/s')


unitsPath=[pathway,'units/'];
units=dir([unitsPath,'*.mat']);
load([unitsPath,units(1).name]);

tmp=dir([pathway,'matlabprot/*',codeWord,'_checker_h_seq_*.mat']);

load([pathway,'matlabprot/',tmp.name]);
ind = cellfun(@(x)( ~isempty(x) ), regexp(unit{1, 2}(:, 1), codeWord));
trial=find(ind);

hekapath=[pathway,'HEKA/'];
heka=dir([hekapath,'*',codeWord,'*.phys']);
protocol=read_header_field_heka(hekapath, heka.name, 'Stimulus Protocol');
voidTime=protocol(3,1);% time when the actual presentation began
tmp=regexp(heka.name,'check_bin');
seed=str2double(heka.name(tmp+9));

defaultStream = RandStream.getGlobalStream;
defaultStream.State = s(seed).State;
% generate sequence (may be done once for all runs using the same seed)
nOfFrames=12000; % should be >= than the max length of matlab protocol using this seed
frame=zeros(nOfFrames,1600);

for i=1:nOfFrames
    frame(i,:)=randi([0 1],1,1600);
end


flips=round(FlipTimeStamps*1000);
comb_spikes=[];
break_units=[];
filter_length=500;
omitTime=4000+filter_length+1;% time to omit in ms (stabilizing time)
nUnits=length(units);
nPixels=1600;

for unitN=1:nUnits
    load([unitsPath,units(unitN).name]);
    
    spikes=round(unit{1, 2}{trial, 2}-voidTime);
    spikes(spikes<omitTime)=[];    
    spikes(spikes>flips(end))=[];
    
    comb_spikes=[comb_spikes spikes];
    break_units=[break_units ones(size(spikes))*unitN];    
end
fprintf([codeWord,', total spikes ',int2str(length(comb_spikes)),'\n'])
[sorted_spikes, init_order]=sort(comb_spikes);
break_units=break_units(init_order);

allFramesDurations=diff(flips);

nFrames=17;
if sta==0
    sta=zeros(filter_length,nPixels,nUnits);
end

for i=1
    
    tmp=sorted_spikes(i);
    
    last_frame=find(flips>=tmp,1)-1;
    twin_spikes=find(sorted_spikes>flips(last_frame)&sorted_spikes<=flips(last_frame+1));
    master_units=break_units(twin_spikes);
    
    
    endTimes=sorted_spikes(twin_spikes)-flips(last_frame);
    first_frame=last_frame-nFrames+1;
    framesDurations=allFramesDurations(first_frame:last_frame);
    cumFrames=[0 cumsum(framesDurations)];
    
    actualFrames=zeros(cumFrames(end),1600);
    tmp=frame(first_frame:last_frame,:);
    for j=1:nFrames
        t=repmat(tmp(j,:),framesDurations(j),1);
        actualFrames(cumFrames(j)+1:cumFrames(j+1),:)=t;
    end
    old_last_frame=last_frame;
    
    for j=1:length(twin_spikes)
        sta(:,:,master_units(j))=sta(:,:,master_units(j))+actualFrames((end-filter_length+1:end)-endTimes(j),:);
    end
    i=i+length(twin_spikes);
end

cntSup=ones(1,nUnits);
sta_temp=zeros(500,1600,nUnits,100);
k=1;
ak=tic;
while i<length(sorted_spikes)
    tmp=sorted_spikes(i);

    last_frame=find(flips>=tmp,1)-1;
    twin_spikes=find(sorted_spikes>flips(last_frame)&sorted_spikes<=flips(last_frame+1));
    master_units=break_units(twin_spikes);
    endTimes=sorted_spikes(twin_spikes)-flips(last_frame);
    
    frame_shift=last_frame-old_last_frame;
    first_frame=first_frame+frame_shift;
    
    actualFrames(1:cumFrames(frame_shift+1),:)=[];

    framesDurations=allFramesDurations(first_frame:last_frame);
    cumFrames=[0 cumsum(framesDurations)];  
    
    tmp=frame(old_last_frame+1:last_frame,:);
    for j=1:frame_shift
        t=repmat(tmp(j,:),framesDurations(end-frame_shift+j),1);
        actualFrames=[actualFrames; t];
    end
    
    old_last_frame=last_frame;

    for j=1:length(twin_spikes)
        sta_temp(:,:,master_units(j),cntSup(master_units(j)))=actualFrames((end-filter_length+1:end)-endTimes(j),:);
        cntSup(master_units(j))=cntSup(master_units(j))+1;   
    end
    if max(cntSup)>92
        sta=sta+sum(sta_temp,4);
        cntSup=ones(1,nUnits);
        sta_temp=zeros(500,1600,nUnits,100);
        if i>5000*k            
            fprintf(['processed ', int2str(i),' spikes in ',num2str(toc(ak)),'s\n'])
            ak=tic;
            k=k+1;
        end
    end
        

    i=i+length(twin_spikes);
end

res=sta;
toc(ak);

% 
% % v1
% while i<length(sorted_spikes)
%     tmp=sorted_spikes(i);
% 
%     last_frame=find(flips>=tmp,1)-1;
%     twin_spikes=find(sorted_spikes>flips(last_frame)&sorted_spikes<=flips(last_frame+1));
%     master_units=break_units(twin_spikes);
%     endTimes=sorted_spikes(twin_spikes)-flips(last_frame);
%     
%     frame_shift=last_frame-old_last_frame;
%     first_frame=first_frame+frame_shift;
%     
%     actualFrames(1:cumFrames(frame_shift+1),:)=[];
% 
%     framesDurations=allFramesDurations(first_frame:last_frame);
%     cumFrames=[0 cumsum(framesDurations)];  
%     
%     tmp=frame(old_last_frame+1:last_frame,:);
%     for j=1:frame_shift
%         t=repmat(tmp(j,:),framesDurations(end-frame_shift+j),1);
%         actualFrames=[actualFrames; t];
%     end
%     
%     old_last_frame=last_frame;
%     
%     the most time-consuming part!
%     tic
%     for j=1:length(twin_spikes)
%         sta(:,:,master_units(j))=sta(:,:,master_units(j))+actualFrames((end-filter_length+1:end)-endTimes(j),:);
%     end
%     toc
%     i=i+length(twin_spikes);
% end
% 
% res=sta;

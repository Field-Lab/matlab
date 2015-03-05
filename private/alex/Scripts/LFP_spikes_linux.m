%% get filtered versions of LFP and spike times
clear all
cd('F:\templates_code')
date='20120628';
mainPath='F:\20120628\MEA_bin\channels\channel';
heka=dir('F:\20120628\HEKA\');
names=dir([mainPath,'1']);

clear file_list
% quick
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'quick', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=72500;
prefix='quick';
% inds='901_1800';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end



clear file_list
% ffflicker
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'ffflicker', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=246000;
prefix='ffflicker';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end



clear file_list
% chirp
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'chirp', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=22500;
prefix='chirp';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end


clear file_list
% checkers
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'check', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=486000;
prefix='check';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end


clear file_list
% NatMov1
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'NatMov1', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=25500;
prefix='nm1';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end


clear file_list
% NatMov2
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'NatMov2', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=23200;
prefix='nm2';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end


clear file_list
% NatMov3
cnt=1;
for i=1:length(heka)-2
    if ~isempty(regexp(heka(i+2).name,'NatMov3', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;
    end
end
length_expected=63500;
prefix='nm3';
inds='all';
for ch=1:60
    tic
    get_splitted(ch,file_list,date,length_expected,prefix,inds);
    toc
end


%do: 0715, 0710, 0411
clear
date='20120411';
chs=4:60;chs(15)=[];
for ch=chs
    tic
    get_spikes_v2(ch,date);
    toc
end



%%%%%%% GET READY %%%%%%%%


% GET READY
clear
cd('F:\20120628\processing')

protocol=read_header_field_heka('F:\20120628\HEKA','20120628_C1#0002_hlc2_quick_BG1_FWND5ND3.phys', 'Stimulus Protocol');
stim=round(protocol(2:21,1));

% find files with quicks
names=dir('F:\20120628\HEKA\*.phys');
cnt=1;
clear file_list
for i=1:177
    if ~isempty(regexp(names(i).name,'quick', 'once'))
        file_list(cnt)=i;
        cnt=cnt+1;        
    end    
end

NoF=length(file_list);
OnLFP=0;
OffLFPminima=0;
OffLFPmaxima=0;
for ch=1:60
    ch
    tic
    if ch~=15
        %lfp
        load(['F:\20120628\processing\preprocessed\lfp\quick_all_CH',int2str(ch),'_0628'])
        whiteFlash=0;
        blackFlash=0;
        for st=1:4:20
            whiteFlash=whiteFlash+lfp(stim(st):stim(st)+3999,:);
            blackFlash=blackFlash+lfp(stim(st+2):stim(st+2)+3999,:);
        end
        clear lfp        
        OnLFP=OnLFP+reshape([min(whiteFlash(1:1984,:)); min(blackFlash(1984:end,:))],NoF*2,1)/59;        
        OffLFPminima=OffLFPminima+reshape([min(whiteFlash(1984:end,:)); min(blackFlash(1:1984,:))],NoF*2,1)/59;        
        OffLFPmaxima=OffLFPmaxima+reshape([max(whiteFlash(1984:end,:)); max(blackFlash(1:1984,:))],NoF*2,1)/59;     
    end
    toc
end

% plotting
figure
subplot(3,1,1)
plot(-OnLFP)
for i=10.5:10:106
    line([i,i],[min(-OnLFP),max(-OnLFP)],'color','r');
end
axis tight
subplot(3,1,2)
plot(-OffLFPminima)
for i=10.5:10:106
    line([i,i],[min(-OffLFPminima),max(-OffLFPminima)],'color','r');
end
axis tight
subplot(3,1,3)
plot(OffLFPmaxima)
for i=10.5:10:106
    line([i,i],[min(OffLFPmaxima),max(OffLFPmaxima)],'color','r');
end
axis tight

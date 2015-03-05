function get_file(flag)

global hekaNames fileName pathName savePath date ch experimenter thrs...
    spikes spikeTimes timeUnit

if isempty(pathName)
    pathName='/mnt/muench_data/data/';
end
persistent hekaPath key_seq wfs


if flag==1
    sep=filesep;
    [fileName,pathName] = uigetfile([pathName,'*.mat']);
    
    a=find(pathName==sep,2,'last');
    savePath=pathName(1:a(1));
    hekaPath=([pathName(1:a(1)),'HEKA',sep]);
    a=find(savePath==sep,2,'last');
    date=pathName(a+1:a+8);
    
    % which format?
    if ~isempty(regexp(fileName, 'waveforms_ch', 'once'))
        key_seq='ch';
        tmp = regexpi(fileName,'ch(\d+)','tokens');       
        wfs=0;
    elseif ~isempty(regexp(fileName, 'CH', 'once'))   % my format
        key_seq='CH';
        tmp = regexpi(fileName,'CH(\d+)','tokens');
        wfs=1;
    else
        display('Unknown format!')
        return
    end    
    ch=int2str(eval(char(tmp{1})));
    
    prompt = {'Path to save','Path to HEKA','Date:','Channel:','Experimenter','Time Unit (1-s,1000-ms)','HEKA'};
    dlg_title = 'File info';
    num_lines = 1;
    def = {savePath,hekaPath,date,ch,'A','1','all'};
    options.Resize='on';
    answer = inputdlg(prompt,dlg_title,num_lines,def,options);
    savePath=answer{1};
    hekaPath=answer{2};
    date=answer{3};
    ch=str2num(answer{4});
    experimenter=answer{5};
    timeUnit=answer{6};
    takeHeka=answer{7};    
    
    if strcmp(takeHeka,'all')
        startHeka=1;
    else
        startHeka=str2num(takeHeka);
    end
    heka=dir([hekaPath,'*.phys']);
    hekaNames=cell(length(heka)-startHeka+1,1);
    for i=startHeka:length(heka)
        hekaNames{i-startHeka+1}=heka(i).name;
    end
    clear heka  
    
else
    ch=ch+1;    
end

if ch==15
    display('Channel 15 is omitted!')
    ch=ch+1;
end
if ch>60
    display('You have only 60 channels!')
    ch=1;
end

% load file from appropriate channel
[tmp,b,c]=regexpi(fileName,[key_seq,'(\d+)'],'tokens','start','end');
new_ch=eval(char(tmp{1}));
if new_ch~=ch % channel changed!
    if wfs
        new_ch=int2str(ch);
    else
        new_ch=int2str(ch);
        if ch<10
            new_ch=['0', new_ch];
        end
    end
else
    new_ch=char(tmp{1});
end

fileName=[fileName(1:b+1),new_ch,fileName(c+1:end)];

% tmp=fileName(a+3);
% if tmp=='_'
%     fileName=[fileName(1:a+1) int2str(ch) fileName(a+3:end)];
% else
%     fileName=[fileName(1:a+1) int2str(ch) fileName(a+4:end)];
% end

load([pathName,fileName]); % load spikes

if wfs % convert format
    suppl_fileName=[fileName(1:end-4), '_spikeTimes.mat']; % load spike times and thresholds
    load([pathName,suppl_fileName]); % load spike times
else
    load([pathName,'spikeTimes_ch',new_ch,'.mat']); % load spike times
    thrs=load([pathName,'thrs_ch',new_ch,'.mat']); % load thresholds
    thrs=thrs.threshold;
    convert_wfs(waveform,spiketime);
end


set_thr(1)

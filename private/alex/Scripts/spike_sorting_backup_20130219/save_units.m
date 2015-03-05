function res=save_units(date,savePath,spikeTimes,thr,thrs,ch,quality,experimenter,hekaNames,sorting)
global rgc getSlice1 getSlice2 notesUI autoPathUI unitUI src
% date - date of experiment
% savePath - path to save units
% spikeTimes - cell array of time stamps of all spikes, splitted by stim files
% thrs - nx1 matrix of thresholds applied for each of the stim files
% ch - channel number (integer)
% quality - quality of sorting. should be queried after tracing (make a
% text field). subj. estimation
% experimenter - one letter
% hekaNames - heka files names
% sorting - axes displayed (get values from the fields)
quality=str2num(quality);

if size(quality,2)<size(rgc,2)
    display('INPUT: quality!')
    return
end
sep=filesep;
notes=get(notesUI,'String');

timeUnit=str2num(get(unitUI,'String'));

if ~get(autoPathUI,'Value')
    savePath=uigetdir;
end

for i=1:2
    switch sorting(i)
        case 1
            par='pc1';
        case 2
            par='pc2';
        case 3
            a=get(getSlice1,'String');
            par=['slice 1: point ',a];
        case 4
            a=get(getSlice2,'String');
            par=['slice 2: point ',a];
        case 5
            par='time';
        case 6
            par='minima';
        case 7
            par='maxima';
        case 8
            par='template'
    end
    if i==1
         sort_parameters=par;
    else
        sort_parameters=[sort_parameters ' against ' par];
    end
end

if ~exist([savePath,'units'],'dir')
    mkdir([savePath,'units'])
end


res=cell(size(rgc,2),1);
if size(spikeTimes,2)<size(spikeTimes,1)
    spikeTimes=spikeTimes';
end
ss=size(spikeTimes,2);
linTim=[];
for i=1:ss
    linTim=[linTim [spikeTimes{i}; ones(1,length(spikeTimes{i}))*i]];
end

% --- which units are already there? ---
so_far=dir([savePath,'units',sep,'*_unit_*.mat']);
last_unit=0;
for i = 1:length(so_far)
    tmp = regexpi(so_far(i).name,'_unit_(\d+)','tokens');
    last_unit = max(last_unit,eval(char(tmp{1})));
end
new_units = last_unit+1:1:last_unit+size(rgc,2);

% --- has this channel been sorted before? ---
last_sort=0;
for i = 1:length(so_far)
    tmp = regexpi(so_far(i).name,strcat(int2str(ch),'_sort(\d+)_'),'tokens');
    if isempty(tmp)
        tmp=0;
    else
        tmp=str2num(char(tmp{1}));
    end
    last_sort = max(last_sort,tmp);
end
sort = last_sort+1;

unitNames=[];

for j=1:size(rgc,2)
    lin_cl=linTim(:,rgc{j});
    unit_cur=['000', int2str(new_units(j))];
    fileName=[experimenter,date,'_CH',int2str(ch),'_sort',int2str(sort),'_',int2str(quality(j)),'_unit_',unit_cur(end-3:end)];
    unit=cell(1,2);
    unit{1,2}=cell(ss,2);
    unit{1,1}{1,1}='threshold';
    unit{1,1}{1,2}=thr;
    unit{1,1}{2,1}='sorting';
    unit{1,1}{2,2}=sort_parameters;
    unit{1,1}{3,1}='channel';
    unit{1,1}{3,2}=ch;
    unit{1,1}{4,1}='number of clusters';
    unit{1,1}{4,2}=size(rgc,2);
    unit{1,1}{5,1}='time unit (1=sec, 1000=ms)';
    unit{1,1}{5,2}=timeUnit;
    unit{1,1}{6,1}='Quality';
    unit{1,1}{6,2}=quality(j);
    unit{1,1}{7,1}='All thresholds';
    unit{1,1}{7,2}=thrs;
    unit{1,1}{8,1}='Notes';
    unit{1,1}{8,2}=notes;
    for i=1:ss
        unit{1,2}{i,1}=hekaNames{i};
        unit{1,2}{i,2}=lin_cl(1,lin_cl(2,:)==i);
    end
    
    save([savePath,'units',sep,fileName],'unit');
    res{j,1}=fileName;
    unitNames=[unitNames,'_',unit_cur];
end
sorting_info=[savePath,'sorting_info',sep,'CH',int2str(ch),sep];
if ~exist(sorting_info,'dir')
    mkdir(sorting_info);
end


saveas(src,[sorting_info,experimenter,date,'_CH',int2str(ch),'_units',unitNames,'_gui.fig']);
h=findobj('Name','Spike Minima Histograms');
saveas(h,[sorting_info,experimenter,date,'_CH',int2str(ch),'_units',unitNames,'_minima.fig']);
h=findobj('Name','ISI <5ms; % - less than 2ms');
saveas(h,[sorting_info,experimenter,date,'_CH',int2str(ch),'_units',unitNames,'_ISI.fig']);
h=findobj('Name','Units mean');
saveas(h,[sorting_info,experimenter,date,'_CH',int2str(ch),'_units',unitNames,'_mean_STD.fig']);

figure(src)
strSaved=[int2str(size(rgc,2)),' unit(s) saved in ',savePath];

textSavedUnitsUI=uicontrol('style','text','string',strSaved,...
    'Units','normalized','position',[0.13 0.018 0.51 0.017],'BackgroundColor',[0.8 0.2 0.2]);


% str=sprintf('Saved %d unit(s) in %sunits', size(rgc,2),savePath);
% 
% choice = questdlg(['Good job.    ',str,' .     Would you like to work more?'], ...
%  'work more!','Yes, let me choose file','Yes, get next channel','Maybe','Maybe');
% 
% switch choice
%     case 'Yes, let me choose file'
%         get_file(1);
%     case 'Yes, get next channel'
%         get_file(2),
%     case 'Maybe'
%         a=ceil(rand(1,1)*3);
%         if a<3
%             get_file(ceil(rand(1,1)*3));
%         else
%             msgbox('happy hamster!!!')
%             pause(1)
%             close all
%             clear all
%             clear functions
%             clear java
%             clear classes
%             pack
%         end
% end
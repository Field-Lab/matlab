%% Fig.O
clear

a = dir('/Users/alexth/Desktop/old_stuff/my_data/20130302/units/*.mat');

for uns = 1:length(a)
%     load('/Users/alexth/Desktop/old_stuff/my_data/20130302_1/units/A20130302_CH52_sort1_100_unit_0030.mat')
    load(['/Users/alexth/Desktop/old_stuff/my_data/20130302/units/', a(uns).name])
    
    hekapath='/Users/alexth/Desktop/old_stuff/my_data/20130302/HEKA/';
    heka=dir([hekapath,'*.phys']);
    
    
    codeWord='NatMov';
    % select heka files
    file_list=[];
    for i=1:length(heka)
        if ~isempty(regexp(heka(i).name,codeWord, 'once'))
            file_list=[file_list i];
        end
    end
    nonEmpty = find(~cellfun('isempty', unit{1, 2}(:,1)),1,'last');
    stimFiles = find(cellfun(@(x)( ~isempty(x) ), regexp(unit{1, 2}(1:nonEmpty,1), codeWord)));
    nds = find(cellfun(@(x)( ~isempty(x) ), regexp(unit{1, 2}(1:nonEmpty,1), 'nd%nd')));
    for j=1:length(nds)
        [~,nd] = regexp(unit{1, 2}{nds(j,1),1},'(?<=nd%nd_)\d','tokens','match');
        nds(j,2)=str2num(nd{1});
    end
    ndFullList=zeros(size(stimFiles));
    
    
    cnt = 8; % number corresponds to NDF
    nm = cell(8,18);
    for i=1:18:144
        for j = 1:18
            spikes=round(cell2mat(unit{1,2}(file_list(i+j-1),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms
            nm{cnt,j} = spikes;
        end
        cnt=cnt-1;
    end
    
    
    % plot raster
    cnt = 1;
    figure
    for i=4:8
        subplot(5,1,cnt)
        hold on
        for j=1:15
            tmp = nm{i,j};
            %             tmp(tmp>7000) = [];
            %             tmp(tmp<2500) = [];
            plot(tmp,zeros(length(tmp),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
        end
        axis off
        cnt = cnt+1;
        
    end
    
end

load('/Users/alexth/Desktop/old_stuff/my_data/analysis/summary_all','bf','wf', 'names')
cell_ID = 387;

codeWord='quick';
nonEmpty = find(~cellfun('isempty', unit{1, 2}(:,1)),1,'last');
stimFiles = find(cellfun(@(x)( ~isempty(x) ), regexp(unit{1, 2}(1:nonEmpty,1), codeWord)));
nds = find(cellfun(@(x)( ~isempty(x) ), regexp(unit{1, 2}(1:nonEmpty,1), 'nd%nd')));
for j=1:length(nds)
    [~,nd] = regexp(unit{1, 2}{nds(j,1),1},'(?<=nd%nd_)\d','tokens','match');
    nds(j,2)=str2num(nd{1});
end
ndFullList=zeros(size(stimFiles));
for j=1:length(stimFiles)
    ndFullList(j)=nds(find(nds(:,1)<stimFiles(j),1,'last'),2);
    if j>1 && ndFullList(j)~=1 && ndFullList(j-1)==1
        ndFullList=ndFullList(1:j-1);
        break;
    end
end
stimFiles=stimFiles(1:length(ndFullList));

protocol=read_header_field_heka(hekapath, heka(stimFiles(10)).name, 'Stimulus Protocol');
protocol=round(protocol(3:2:end-2,1));

black_flash_spike_data = cell(5,3,8);
white_flash_spike_data = cell(5,3,8);

for j=1:length(stimFiles)
        spikes=round(unit{1,2}{stimFiles(j),2}*1000/unit{1,1}{5,2}); %now in ms
        for st=1:2:10
            black_flash_spike_data{st/2+0.5,mod(j-1,3)+1,ndFullList(j)} = spikes(spikes>protocol(st+1,1)-499 & spikes<protocol(st+1,1)+4000) - (protocol(st+1,1)-499);
            white_flash_spike_data{st/2+0.5,mod(j-1,3)+1,ndFullList(j)} = spikes(spikes>protocol(st,1)-499 & spikes<protocol(st,1)+4000) - (protocol(st,1)-499);
            
        end
end

% plot raster
figure
cnt = 1;
for i=4:7
    
    tmp = black_flash_spike_data(:,:,i);
    tmp = tmp(:);
    
    subplot(4,2,cnt)
    hold on
    for j=1:15
        plot(tmp{j},zeros(length(tmp{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    
    tmp = white_flash_spike_data(:,:,i);
    tmp = tmp(:);
    axis off
    
    subplot(4,2,cnt+1)
    hold on
    for j=1:15
        plot(tmp{j},zeros(length(tmp{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    axis off    
    cnt = cnt+2;
    
end

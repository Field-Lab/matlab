%% get responses
clear
mainpath='/Users/alexth/Desktop/old_stuff/my_data/';
date = '20130225';
filename = 'A20130225_CH16_sort1_90_unit_0019';

codeWord='quick';
baseline=1;

units=dir([mainpath,date,'/units/*.mat']);
load([mainpath,date,'/units/',filename]);

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

if sum(ndFullList==9)>9
    ndFullList(1)=[];
    stimFiles(1)=[];
end


heka=dir([mainpath,date,'/HEKA/*phys']);
protocol=read_header_field_heka([mainpath,date,'/HEKA/'], heka(stimFiles(1)).name, 'Stimulus Protocol');
protocol=round(protocol(2:2:end-2,1));

black_flash = nan(5500,5,9,9);
white_flash = nan(5500,5,9,9);

black_flash_spike_data = cell(5,9,9);
white_flash_spike_data = cell(5,9,9);

for j=1:length(stimFiles)
        spikes=round(unit{1,2}{stimFiles(j),2}*1000/unit{1,1}{5,2}); %now in ms
        
        conv=convolved(spikes,40,75000);
        conv=conv(121:end-120);
        for st=1:2:10
            black_flash(:,st/2+0.5,mod(j-1,9)+1,ndFullList(j))=conv(protocol(st+1,1)-1499:protocol(st+1,1)+4000)';
            white_flash(:,st/2+0.5,mod(j-1,9)+1,ndFullList(j))=conv(protocol(st,1)-1499:protocol(st,1)+4000)';
            
            black_flash_spike_data{st/2+0.5,mod(j-1,9)+1,ndFullList(j)} = spikes(spikes>protocol(st+1,1)-1499 & spikes<protocol(st+1,1)+4000) - (protocol(st+1,1)-1499);
            white_flash_spike_data{st/2+0.5,mod(j-1,9)+1,ndFullList(j)} = spikes(spikes>protocol(st,1)-1499 & spikes<protocol(st,1)+4000) - (protocol(st,1)-1499);
            
        end
end


bf=reshape(black_flash,5500,5*9,9);
mean_in_nd = squeeze(mean(bf,2));
figure
for i=1:8
    subplot(8,1,i)
    plot(mean_in_nd(:,i))
end


wf=reshape(white_flash,5500,5*9,9);
mean_in_nd = squeeze(mean(wf,2));
figure
for i=1:8
    subplot(8,1,i)
    plot(mean_in_nd(:,i))
end

% plot raster
cnt = 1;
figure
for i=4:8
    
    tmp = black_flash_spike_data(:,:,i);
    tmp = tmp(:);
    
    subplot(5,2,cnt)
    hold on
    for j=20:45
        plot(tmp{j},zeros(length(tmp{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    
    tmp = white_flash_spike_data(:,:,i);
    tmp = tmp(:);
    axis off
    
    subplot(5,2,cnt+1)
    hold on
    for j=20:45
        plot(tmp{j},zeros(length(tmp{j}),1)+0.8/13*j-0.8/13, '.k', 'markersize',0.1)
    end
    cnt = cnt+2;
    axis off
    
end


%% Full Field Flicker Filter PARTIAL
clear
cd('/mnt/muench_data/user/alexandra/scripts')
date='20130301';
codeWord='HL10';
howMany=7;

mainpath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/'];
hekapath=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/HEKA/'];
heka=dir([hekapath,'*.phys']);

if ~exist([mainpath 'LinearFilters/'],'dir')    
    mkdir([mainpath, 'LinearFilters']);
end
path2save=[mainpath 'LinearFilters/'];

% make list of files to put into easy formatted unit (choosing by the name of the stimulus in heka files)
file_list=[];
for i=1:length(heka)
    if ~isempty(regexp(heka(i).name,codeWord, 'once'))
        file_list=[file_list i];
    end
end

timingPath=[mainpath 'timing/'];
timings=dir([timingPath,'*.mat']);

%get list of units to be reformatted
basic_format_units_list=dir([mainpath, 'units/*.mat']);

% get protocol
protocol=read_header_field_heka(hekapath, heka(file_list(1)).name, 'Stimulus Protocol');
protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
protocol=round(protocol(2:end,[1,4]));
protocol_accum=zeros(size(protocol,1),2,length(file_list));
abused_accum=zeros(size(protocol,1),length(file_list));
for i=1:length(file_list)
        protocol=read_header_field_heka(hekapath, heka(file_list(i)).name, 'Stimulus Protocol');
        protocol(protocol(:,3)==0&protocol(:,4)==0,4)=-1;
        protocol=round(protocol(2:end,[1,4]));
        protocol_accum(:,:,i)=protocol;
        load([timingPath, timings(file_list(i)).name]);
        abused=abused(2:end,5)-abused(2,5);
        abused_accum(:,i)=abused(3:end)*1000;
end


flicker=zeros(protocol_accum(end,1,1),length(file_list));
startingPoints=zeros(1,length(file_list));
maxLength=zeros(1,length(file_list));
changeContrast=zeros(howMany,length(file_list));
for i=1:length(file_list)
    flips=protocol_accum(1:end-2,1,i);
    pxls=(protocol_accum(1:end-2,2,i)-30)/30;
    startingPoints(i)=flips(1);
    flips=flips-flips(1)+1;    
    
    tmp=zeros(flips(end),1);
    pxls(pxls==0)=100;
    tmp(flips(2:end),1)=pxls(2:end);
    % fill stimulus down to ms with brightness values where spikes
    % happen
    droppedMax=max(diff(flips));
    for j=1:droppedMax+1
        subFrame=flips(2:end)-j;
        subFrame(subFrame<1)=[];
        subFrame(find(tmp(subFrame,1)))=[];
        tmp(subFrame,1)=tmp(subFrame+1,1);
    end
    tmp(tmp(:,1)==100,1)=0;
    flicker(1:length(tmp),i)=tmp;
    maxLength(i)=length(tmp);
    changeContrast(:,i)=flips(1:600:end);
end
flicker=flicker(1:min(maxLength),:);

save([path2save,date,'_flicker_',codeWord],'flicker','startingPoints','maxLength','changeContrast')


filter_length=500;

names=cell(1,length(basic_format_units_list));
LinearFilterHigh=zeros(filter_length,length(file_list),length(basic_format_units_list));
LinearFilterLow=zeros(filter_length,length(file_list),length(basic_format_units_list));
SpikeCount=zeros(length(file_list),length(basic_format_units_list),2);

for cnt=1:length(basic_format_units_list)
    basic_format_units_list(cnt).name
    load([mainpath,'units/',basic_format_units_list(cnt).name]);
    for i=1:length(file_list)
        spikes=round(cell2mat(unit{1,2}(file_list(i),2))*1000/cell2mat(unit{1,1}(5,2))); %now in ms        
        spikes=spikes-startingPoints(i)+1;
        spikes(spikes<filter_length|spikes>size(flicker,1))=[];
        n=zeros(length(spikes),filter_length);
        for k=1:filter_length
            n(:,k)=flicker(spikes-k+1,i);
        end
        spikesHigh=zeros(1,size(n,1));
        for m=1:2:howMany-1
            spikesHigh(spikes>changeContrast(m,i)&spikes<=changeContrast(m+1,i))=1;
        end
        spikesHigh=logical(spikesHigh);
        LinearFilterHigh(1:filter_length,i,cnt)=sum(n(spikesHigh,:));
        LinearFilterLow(1:filter_length,i,cnt)=sum(n(~spikesHigh,:));
        SpikeCount(i,cnt,1)=sum(spikesHigh);
        SpikeCount(i,cnt,2)=sum(~spikesHigh);
        names{cnt}=basic_format_units_list(cnt).name;
    end
end
save([path2save,date,'_LF_',codeWord],'LinearFilterHigh','LinearFilterLow','SpikeCount','names')






%% Plot

clear
date='20130301';
nd='87654321';
codeWord='HL10';
load(['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/',date,'_LF_',codeWord]);

path2save=['/mnt/muench_data/data/alexandra/MEA_data/',date,'/LinearFilters/'];


figure(1)
set(gcf,'Position', [1 352 1011 625])
[rows,cols]=opt_subplots(length(nd));

for i=1:size(LinearFilter,3)
    tmp=LinearFilterLow(:,:,i);
    minn=min(tmp(:));
    maxx=max(tmp(:));
    for j=1:2:size(tmp,2)
        k=subplot(rows,cols,floor(j/2)+1);
        plot(tmp(:,j:j+1))
        line([0 500],[0 0],'color','k')
        title(['ND',nd(floor(j/2)+1)])
        axis([0 500 minn maxx])
    end
    subplot('Position',[0.5 0.97 0.00001 0.00001])
    set(gca,'XTick',0,'YTick',0,'XTickLabel','','YTickLabel','')
    title([names{i}(1:end-4),'_',codeWord],'FontSize',12,'FontWeight','bold','Interpreter','None')

    saveas(gcf,[path2save,names{i}(1:end-4),'_',codeWord,'.png'])
end
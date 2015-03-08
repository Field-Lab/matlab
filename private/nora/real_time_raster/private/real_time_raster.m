function real_time_raster(rawdatafiledir,thresholding_file,raster_starts,raster_length)
%
% real_time_raster(rawdatafiledir,thresholding_file,raster_starts,raster_length)
% rawdatafiledir is the location of the raw data file
% thresholding_file is the file to read the thresholding data from
% raster_starts: the trigger numbers where the rasters start. must be
%           integer and nonzero
% raster_length: in seconds
%
% In a continuous run, raster_starts= 
% time we want to start rasters*120/100
%
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar

sampling_rate=20000;
load(thresholding_file);
channel_count=length(channels.thresh);
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);

% get all triggers up to the beginning of the last raster
disp('Finding triggers...')
tic
trigger_times=zeros(raster_starts(end),1); % vector to store the trigger times
trigger_increment=1.6652e+04; % about where the next trigger is supposed to be! 100 frames later
trigger_count=1;
sample_start=0;
sample_end=sampling_rate/2;
trigger_data=rdf.getData(0,sample_start,sample_end-sample_start);
error_count=0;

while trigger_count<=raster_starts(end)
    
    % find the first trigger in the loaded data
    temp=find(trigger_data~=0,1);
    
    % set sample_start and sample_end for the next loop 
    if isempty(temp)
        sample_start=sample_end;
        sample_end=sample_end+sampling_rate;
        error_count=error_count+1;
    else
        trigger_times(trigger_count)=temp+sample_start;
        sample_start=trigger_times(trigger_count)+trigger_increment;
        sample_end=sample_start+5;
        trigger_count=trigger_count+1;
    end
    
    % 
    trigger_data=rdf.getData(0,sample_start, sample_end-sample_start);
    
end

if any(diff(trigger_times)<=0)
    save('trigger.mat','trigger_times');
    error('Triggers not in order???')
end
disp([num2str(error_count) ' triggers were not in the expected place']);

toc

disp('Thesholding...')

tic
% get spikes times
spike=cell(channel_count,length(raster_starts));
warn=0;

% one raster at a time
for i=1:length(raster_starts)
    
    % select channels and subtract threshold
    for j=1:channel_count
        selected_data=-double(rdf.getData(channels.number{j},trigger_times(raster_starts(i)),raster_length*sampling_rate))-channels.thresh{j};
        
        % make binary
        selected_data(selected_data>0)=1;
        selected_data(selected_data<0)=0;
        temp=diff(selected_data);
        
        % store spike times
        spike{j,i}=find(temp==1);
        if isempty(spike{j,i})
            figure(2)
            plot(-double(rdf.getData(channels.number{j},trigger_times(raster_starts(i)),raster_length*sampling_rate))-channels.thresh{j});
            figure(3)
            plot(-double(rdf.getData(channels.number{j},sampling_rate,sampling_rate))-channels.thresh{j});
            spike{j,i}=0;
            warn=1;
        end
    end
    
end

if warn
    warning('No spikes found in one or more rasters. May have picked a bad channel')
end
rdf.close()
toc

% plot rasters
scrsz = get(0,'ScreenSize');
h=figure(1);
set(h,'Position',[1 scrsz(4)/1.5 scrsz(3)/1.5 scrsz(4)/1.5]);
for j=1:channel_count
    subplot(ceil(channel_count/2),2,j)
    if j==2*ceil(channel_count/2)-1
        ylabel('Trial')
        xlabel('Time (seconds)')
        set(gca,'XTick',(0:raster_length)*sampling_rate,'XTickLabel',0:raster_length)
    else
        set(gca,'XTick',[],'YTick',[])
    end
    hold on
    for i=1:length(raster_starts)
        X=[1,1]'*spike{j,i}';
        Y=[i-0.4,i+0.4]'*ones(length(X),1)';
        l=line(X,Y);
        set(l,'Color','black')
    end
    hold off
    ylim([0 length(raster_starts)+1])
end
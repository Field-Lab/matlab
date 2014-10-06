function stream_raster(rawdatafiledir,thresholding_file,raster_starts, raster_length, previous_streaming_length)
%
% real_time_raster(rawdatafiledir,thresholding_file,raster_starts,raster_length,previous_streaming_length)
%
% INPUTS
% rawdatafiledir can be either the location of the raw data files or one specific
%           file
% thresholding_file is the file to read the thresholding data from
% raster_starts: the trigger numbers where the rasters start. must be
%           integer and nonzero
% raster_length: in seconds
% previous_streaming_length: underestimate! approximately how long streaming has been
%           running. This prevents the code from hitting the end of the data file. If the
%           data is done streaming, put a really long time for previous
%           streaming time and it will do it all at once
%
% EXAMPLE
% stream_raster('/Volumes/stream-bertha/Data/9999-99-99-9/data014.bin','thresh.mat',1:10:200,2,20)
% creates a raster for 9999-99-99-9/data014 using the thresholds in thresh.mat,
% with the raster starts every 10 seconds from 1 to 200, and 2 second
% rasters,assuming the streaming has been running for at least 20 seconds
%
% TRIGGERS
% In a continuous run, raster_starts= 
% time we want to start rasters*120/100
%
% COMMON ISSUES
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar
% Run stream_thresh FIRST!

% make sure code doesn't get ahead of streaming
tic

% Allocations and stuff
sampling_rate=20000;
load(thresholding_file);
channel_count=length(channels.thresh);
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
spike=cell(channel_count,length(raster_starts));
warn=0;
trigger_times=zeros(raster_starts(end),1); % vector to store the trigger times
trigger_increment=1.6652e+04; % about where the next trigger is supposed to be! 100 frames later
trigger_count=1;
sample_start=0;
sample_end=sampling_rate/2;
trigger_data=rdf.getData(0,sample_start,sample_end-sample_start);
error_count=0;
raster_count=1;

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
    
    
    if trigger_count==raster_starts(raster_count)+1

       % pause(raster_length);
        % select channels and subtract threshold
        for j=1:channel_count
            selected_data=-double(rdf.getData(channels.number{j},trigger_times(raster_starts(raster_count)),raster_length*sampling_rate))-channels.thresh{j};
                       
            % make binary
            selected_data(selected_data>0)=1;
            selected_data(selected_data<0)=0;
            temp=diff(selected_data);
            
            % store spike times
            spike{j,raster_count}=find(temp==1);
            if isempty(spike{j,raster_count})
                spike{j,raster_count}=0;
                warn=1;
            end
        end
        if warn
            warning('No spikes found in one or more rasters. May have picked a bad channel or threshold')
            warn=0;
        end
        
        % plot rasters
        scrsz = get(0,'ScreenSize');
        h=figure(1);
        set(h,'Position',[1 scrsz(4)/1.5 scrsz(3)/2.1 scrsz(4)/1.3]);
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
            for i=1:raster_count
                X=[1,1]'*spike{j,i}';
                Y=[i-0.4,i+0.4]'*ones(length(X),1)';
                l=line(X,Y);
                set(l,'Color','black')
            end
            hold off
            ylim([0 length(raster_starts)+1])
        end
           raster_count=raster_count+1;
    end  
    
    % prevent it from hitting the end of the data file
    while sample_end/20000>(toc+previous_streaming_length)
        pause(0.1)
    end
    
    trigger_data=rdf.getData(0,sample_start, sample_end-sample_start);
    
    
end
end
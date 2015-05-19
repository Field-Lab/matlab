function spike = stream_raster(rawdatafiledir,thresholding_file, raster_interval, raster_length, number_of_rasters, previous_spikes)
% stream_raster(rawdatafiledir,thresholding_file,raster_interval,raster_length, number_of_rasters)
%
% INPUTS
% rawdatafiledir can be either the location of the raw data files or one specific
%           file
% thresholding_file is the file to read the thresholding data from
% raster_interval: either one number: the number of triggers in between raster starts
%                       or two numbers: the # of triggers between rasters
%                       and the # of triggers to delay
% raster_length: in seconds
% number_of_rasters: just to preallocate the variable, so overestimats
%
% EXAMPLE
% stream_raster('/Volumes/stream-bertha/Data/9999-99-99-9/data014.bin','thresh.mat',10,2, 100)
% creates a raster for 9999-99-99-9/data014 using the thresholds in thresh.mat,
% with 2 sec rasters starting every 10 seconds. If there are more than 100
% rasters, it will slow down because it will have to reallocate the spikes
% variable.
%
% TRIGGERS
% In a continuous run, triggers=seconds*120/100
%
% COMMON ISSUES
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar
% Run stream_thresh FIRST!

load(thresholding_file);
channel_count=length(channels.thresh);
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);

if nargin == 5
    spike=cell(channel_count,number_of_rasters);
    raster_count=0;
else
    raster_count=size(previous_spikes,2);
    number_of_rasters=number_of_rasters + raster_count ;
    spike=cell(channel_count,number_of_rasters);
    spike(:,1:raster_count) = previous_spikes;
    clear previous_spikes
end

if length(raster_interval) == 2
    offset = raster_interval(2)+1;
    raster_interval = raster_interval(1);
else
    offset=1;
end

trigger_channel=0;
sampling_rate=20000;
trigger_increment=1.6652e+04; % about where the next trigger is supposed to be! 100 frames later

% load the initial half second of data to get started
sample_start=0;
sample_end=sampling_rate/2;
trigger_data=rdf.getData(trigger_channel,sample_start,sample_end-sample_start);

% bookkeeping counts, errors and warnings
bad_channel_warn=0;
trigger_count=0;
missed_trigger=0;
end_of_streaming=0;

while ~end_of_streaming
    
    % find the first trigger in the loaded data
    trigger_time=find(trigger_data~=0,1);
 
    % if there's no trigger, load the next second of data and check there
    if isempty(trigger_time)
        sample_start=sample_end
        sample_end=sample_end+sampling_rate;
        missed_trigger=missed_trigger+1;
        
    % if there is a trigger, reset sample start and end and check if
    % trigger is a raster start
    else
        sample_start=trigger_time+trigger_increment+sample_start;
        sample_end=sample_start+5;
        trigger_count=trigger_count+1;
        % if it is a raster start, find the spikes and plot!
        if trigger_count >= offset && ~mod(trigger_count-offset, raster_interval)
            raster_count=raster_count+1;
            if raster_count > number_of_rasters 
                raster_count=raster_count-1; 
                spike=spike(:,1:raster_count);
                disp('Done'); 
                break; 
            end
            
            disp(['Making Raster ' num2str(raster_count)])

            % select channels and subtract threshold
            for j=1:channel_count
                
                % try to load the raster data, or wait for data to come in
                raster_wait=1;
                while raster_wait
                    try
                        rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
                        selected_data=-double(rdf.getData(channels.number{j},trigger_time+sample_start,raster_length*sampling_rate))-channels.thresh{j};
                        raster_wait=0;
                    catch
                        pause(1)
                        raster_wait=raster_wait+1;
                        if raster_wait>raster_length+2
                            raster_wait=0;
                            end_of_streaming=1;
                            disp('End of streaming')
                        end
                    end
                end
                
                if end_of_streaming
                    spike=spike(:,1:raster_count);
                    break; 
                end
                
                % make binary
                selected_data(selected_data>0)=1;
                selected_data(selected_data<0)=0;
                temp=diff(selected_data);
                
                % store spike times
                spike{j,raster_count}=find(temp==1);
                if isempty(spike{j,raster_count})
                    spike{j,raster_count}=0;
                    bad_channel_warn=1;
                end
            end
            if bad_channel_warn
                warning('No spikes found in one or more rasters. May have picked a bad channel or threshold')
                bad_channel_warn=0;
            end
            
            if end_of_streaming
                spike=spike(:,1:raster_count);
                break; 
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
                    l=line(X,Y,'LineWidth',0.25);
                    set(l,'Color','black')
                end
                hold off
                ylim([0 raster_count+1])
                pause(0.05)
            end
        elseif mod(trigger_count,raster_interval)==1
            disp('Waiting for more data')
        end
    end
    
    % load the next bit of data
    wait_count=1;
    while wait_count>0
        % try loading the next bit of data
        try
            rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
            trigger_data=rdf.getData(trigger_channel,sample_start, sample_end-sample_start);
            wait_count=0;
        % if it doesn't load, wait a bit for more data
        catch
            pause(1)
            wait_count=wait_count+1;
            % if you've been waiting awhile, end the program. 
            if wait_count>10
                wait_count=0;
                end_of_streaming=1;
                spike=spike(:,1:raster_count);
                disp('End of streaming')
            end
        end
    end
      
end
end
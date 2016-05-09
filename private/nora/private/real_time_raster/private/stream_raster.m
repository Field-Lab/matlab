function stream_raster(rawdatafiledir,thresholding_file)
% stream_raster(rawdatafiledir,thresholding_file,raster_interval,raster_length, number_of_rasters)
%
% INPUTS
% rawdatafiledir can be either the location of the raw data files or one specific
%           file
% thresholding_file is the file to read the thresholding data from
%
% EXAMPLE
% stream_raster('/Volumes/stream-bertha/Data/9999-99-99-9/data014.bin','thresh.mat')
% creates a raster for 9999-99-99-9/data014 using the thresholds in thresh.mat
% with the time parameters from raster_params
%
% COMMON ISSUES
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar
% Run stream_thresh FIRST!

% Allocations and stuff
raster_params
number_of_rasters = length(raster_start);
load(thresholding_file);
channel_count=length(channels.thresh);
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
spike=cell(channel_count,number_of_rasters);

% Magic Numbers
sampling_rate=20000;
trigger_increment=1.6652e+04; % about where the next trigger is supposed to be! 100 frames later

% load the initial half second of data to get started
sample_start=0;
sample_end=sampling_rate/2;
trigger_data=rdf.getData(1,sample_start,sample_end-sample_start);

% bookkeeping counts, errors and warnings
bad_channel_warn=0;
trigger_count=0;
missed_trigger=0;
raster_count=0;
end_of_streaming=0;

while ~end_of_streaming
    
    % find the first trigger in the loaded data
    trigger_time=find(trigger_data~=0,1);
 
    % if there's no trigger, load the next second of data and check there
    if isempty(trigger_time)
        sample_start=sample_end;
        sample_end=sample_end+sampling_rate;
        missed_trigger=missed_trigger+1;
        
    % if there is a trigger, reset sample start and end and check if
    % trigger is a raster start
    else
        sample_start=trigger_time+trigger_increment+sample_start;
        sample_end=sample_start+5;
        trigger_count=trigger_count+1;
        
        % if it is a raster start, find the spikes and plot!
        if trigger_count==raster_start(raster_count+1)
            raster_count=raster_count+1
            
            disp(['Making Raster ' num2str(raster_count)])

            % select channels and subtract threshold
            for j=1:channel_count
                
                % try to load the raster data, or wait for data to come in
                raster_wait=1;
                while raster_wait
                    try
                        pause(1)
                        rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
                        selected_data=-double(rdf.getData(channels.number{j},trigger_time+sample_start,raster_length*sampling_rate))-channels.thresh{j};
                        raster_wait=0;
                    catch
                        raster_wait=raster_wait+1;
                        if raster_wait>raster_length+2
                            raster_wait=0;
                            end_of_streaming=1;
                            disp('End of streaming')
                        end
                    end
                end
                
                if end_of_streaming; break; end
                
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
            
            if end_of_streaming; break; end

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
                ylim([0 raster_count+1])
            end
       % elseif trigger_count==(raster_start(raster_count)+1)
        %    disp('Waiting for more data')
        end
    end
    
    % load the next bit of data
    wait_count=1;
    while wait_count>0
        % try loading the next bit of data
        try
            rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
            trigger_data=rdf.getData(0,sample_start, sample_end-sample_start);
            wait_count=0;
        % if it doesn't load, wait a bit for more data
        catch
            pause(1)
            wait_count=wait_count+1;
            % if you've been waiting awhile, end the program. 
            if wait_count>10
                wait_count=0;
                end_of_streaming=1;
                disp('End of streaming')
            end
        end
    end
    
    
end
end
function stream_thresh(rawdatafiledir, thresholding_file, number_of_channels)
%
% stream_thresh(rawdatafiledir, thresholding_file, number_of_channels)
%
% INPUTS
% thresholding_file is the file to save the thresholding data to
% number_of_channels is number of channels on the whole array.
%
% EXAMPLE
% stream_thresh(rawdatafiledir, thresholding_file, number_of_channels)
% stream_thresh('/Volumes/stream-bertha/Data/9999-99-99-9','thresh.mat',512)
%
% COMMON ISSUES
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar
% Make sure streaming was going for at least a few seconds before running
% Run stream_raster AFTER this

clearvars -global channels
global channel_count channels
channel_count=0;
sampling_rate=20000;

% let users pick some channels to look at if they want
user_channels=input('Enter any channels you want to look at as an array or hit enter \n');
if ~isempty(user_channels)
    if max(user_channels)>number_of_channels || any(floor(user_channels)~=user_channels)
        user_channels(user_channels>number_of_channels)=[];
        user_channels(floor(user_channels)~=user_channels)=[];
    end
end
n_user_channels=length(user_channels);

% load and plot a second of data to do thresholding and channel selection
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
temp_data=rdf.getData(sampling_rate,sampling_rate);
fig_h=figure(1);
scrsz = get(0,'ScreenSize');
set(fig_h,'Name','Click the channels you want to monitor')
set(fig_h,'Position',[1 scrsz(4)/1.5 scrsz(3)/2.1 scrsz(4)/1.3]);
sample=randsample(number_of_channels,64);
if ~isempty(user_channels)
    for i=1:n_user_channels
        channel=user_channels(i);
        subplot(8,8,i);
        p=plot(-double(temp_data(:,channel+1)),'r');
        set(p,'HitTest','off')
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig, 'UserData', channel)
        ylim([-500 200])
        title(num2str(channel))
    end
    for i=n_user_channels+1:64
        channel=sample(i);
        subplot(8,8,i);
        p=plot(-double(temp_data(:,channel+1)));
        set(p,'HitTest','off')
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig, 'UserData', channel)
        ylim([-500 200])
    end
else
    for i=1:64
        channel=sample(i);
        subplot(8,8,i);
        p=plot(-double(temp_data(:,channel+1)));
        set(p,'HitTest','off')
        set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig, 'UserData', channel)
        ylim([-500 200])
    end
end
uicontrol('String','Done','Style','pushbutton','Callback',@done_selecting);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks

% thresholding a channel
    function createnew_fig(cb,~)
        
        channel_count=channel_count+1;
        % replot the selected channel
        h2=figure(2);
        hh = copyobj(cb,h2);
        set(h2,'Position',[1 scrsz(4)/1.5 scrsz(3)/2.1 scrsz(4)/1.3]);
        set(gcf,'Name','Click to set threshold')
        set(hh,'ButtonDownFcn',[]);
        set(hh, 'Position', get(0, 'DefaultAxesPosition'));
        
        % save the channel number and threshold for later use
        [~,channels.thresh{channel_count}]=ginput(1);
        channels.number{channel_count}=get(hh,'UserData');
        
        close(figure(2))
        
    end

% once we've selected all the channels we want to use
    function done_selecting(~,~,~)
        close(figure(1))
        save(thresholding_file,'channels')
        
        clearvars -global channels channel_count
    end
end
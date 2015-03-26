function real_time_thresholding(rawdatafiledir, thresholding_file, number_of_channels)
%
% real_time_raster(rawdatafiledir, thresholding_file, number_of_channels)
% thresholding_file is the file to save the thresholding data to
% number_of_channels is number of channels on the whole array.
%
% Needs your path to include something like Vision.app/Contents/Resources/Java/Vision.jar

clearvars -global channels
global channel_count channels
channel_count=0;
sampling_rate=20000;

% load and plot a second of data to do thresholding and channel selection
rdf=edu.ucsc.neurobiology.vision.io.RawDataFile(rawdatafiledir);
temp_data=rdf.getData(sampling_rate,sampling_rate);
fig_h=figure(1);
set(fig_h,'Name','Click the channels you want to monitor')
sample=randsample(number_of_channels,64);
for i=1:64
    channel=sample(i);
    subplot(8,8,i);
    p=plot(-double(temp_data(:,channel+1)));
    set(p,'HitTest','off')
    set(gca,'XTick',[],'YTick',[],'ButtonDownFcn',@createnew_fig, 'UserData', channel)
    ylim([-500 200])
end
uicontrol('String','Done','Style','pushbutton','Callback',@done_selecting);

%%%%%%%%%%%%%%%%%%%%%%%%%
% Callbacks

% thresholding a channel
    function createnew_fig(cb,~)
        
        channel_count=channel_count+1;
        
        % replot the selected channel
        hh = copyobj(cb,figure(2));
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
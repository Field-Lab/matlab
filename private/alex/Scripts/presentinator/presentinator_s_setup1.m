% case 'S' % stimulus: [action frames posX um posY um]

command = strrep(command,'\n',sprintf('\n'));
Screen('FillRect',w,bgcolor);
Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
Screen('Flip',w);

bullseyeOn = 0;
if isempty(command(3:end))           % use last command
    command = last_s_command;
end
last_s_command = command;

if isempty(command(3:end))
    command = ' '; % delete command buffer
    return;
end

rgb = regexp(command(2:end),'[ \f\r\t\v]*([^\n])*[\n]*(.*)','tokens');
command = rgb{1}{2};
protocolFileName = rgb{1}{1};
bmppathstr = fileparts(protocolFileName);

if isempty(command)
    command = ' '; % delete command buffer
    % load protocol: protocolFileName
    if (verbose > 1)
        Screen('FillRect',w,bgcolor);
        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
%         Screen(w,'DrawText',['Load ' protocolFileName], 25, 15, infocolor);
        Screen('Flip',w);
    end
    protocol = [];
    
    fid = fopen(protocolFileName,'r');
    if fid == -1
        txt = ['Missing file ' protocolFileName];
        presentinator_error;
        return;
    end
    params = fgetl(fid); % first line
    params = regexp(params,' (\S+)\%(\d+)\%\S+','tokens');
    tline = fgetl(fid); % bmpdir/comment line
    tline = strtok(tline,sprintf('\t'));
    if exist(tline,'dir')
        bmppathstr = tline; % bmpdir
    end
    
    while true
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end % until EOF
        % replace params
        for pid=1:length(params)
            tline = strrep(tline,params{pid}{1},params{pid}{2});
        end
        clear T;
        tline = regexp(tline,'(\S+)\s+(\S+)\s+(\S*)\s*(\S*)','tokens');
        try
            protocol(end+1,:) = cellfun(@eval, tline{1});
        catch % hopefully only T is missing
            try
                totaltime = eval(tline{1}{2});
                for T = 0:1/frameRateFactor:totaltime-1/frameRateFactor
                    try
                        protocol(end+1,:) = cellfun(@eval, tline{1});
                        protocol(end,2) = 1/frameRateFactor;
                    catch
                        protocol(end+1,:) = [-1 1/frameRateFactor -3 0];
                        txt = ['Unknown parameter? in stim file']; %#ok<NBRAK>
                        presentinator_error;
                    end
                end
            catch
                protocol(end+1,:) = [-1 0 -3 0];
                txt = ['Unknown duration: ' tline{1}{2}];
                presentinator_error;
            end
        end
    end
    fclose(fid);
else
    if (verbose > 1)
        Screen('FillRect',w,bgcolor);
        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
        Screen(w,'DrawText',['Skip ' protocolFileName], 25, 15, infocolor);
        Screen('Flip',w);
    end
    %     fprintf(logFid, '>%s', command);
    protocol = textscan(command,'%f%f%f%f', -1, 'delimiter', sprintf(' \t'), 'multipleDelimsAsOne', 1); % 'headerlines',1
    command = ' '; % delete command buffer
    if isempty(protocol)
        txt = ['Missing command: ' rgb{1}{2}];
        presentinator_error;
        return;
    end
    protocol = [protocol{:}];
end

if size(protocol,2) == 2 % old style file
    protocol(:,3) = 0;
    protocol(:,4) = 0;
end
if size(protocol,2) == 5 % new style file
    
end

% empty line can make NaN
protocol = protocol(~isnan(protocol(:,1)),:);
protocol = protocol(isfinite(protocol(:,2)),:);

% first line for init
protocol = [-1 waitFrames(1) 0 0; protocol]; % wait 1 blink at first
if length(waitFrames) > 1
    protocol = [protocol; -1 waitFrames(2) 0 0];
end
% bmp id 1 based
protocol(:,1) = floor(protocol(:,1) + 1);

% duration in blankings - adjust for frameRateFactor
protocol(:,2) = ceil(protocol(:,2)*frameRateFactor);

%% adjust offsets
% offset in pixel
bmplines = find(protocol(:,1) > 0);
protocol(bmplines,3) = protocol(bmplines,3)/um2pixel;
protocol(bmplines,4) = protocol(bmplines,4)/um2pixel;
% external offset
protocol(bmplines,3) = floor(-protocol(bmplines,3)-offset(1));
protocol(bmplines,4) = floor(-protocol(bmplines,4)+offset(2));

%% check user command exist
so = instrfind();
for i=1:size(so,2)
    so_arr{str2num(so(i).Port(4))} = so(i); %#ok<*ST2NM>
    %TODO? cmd{6} = so_arr{str2double(cmd{2}(4))};    
end

bmplines = find(protocol(:,1) == 0 & protocol(:,3) == 4);

for i=1:length(bmplines)
    if protocol(bmplines(i),4) < 1 || ...
            protocol(bmplines(i),4) > length(user_commands)  || ...
            isempty(user_commands{protocol(bmplines(i),4)})
        txt = ['Unknown user command ' num2str(protocol(bmplines(i),4)) ' in stimulus'];
        presentinator_error;
        protocol(bmplines(i),3) = -3;
    end
end

%% polychrome
bmplines = find(protocol(:,1) == 0 & protocol(:,3) == 5);

if ~isempty(bmplines) % Polychrome is in use
    if ~polychrome.exist
        txt = ['TILL Polychrome dll is missing']; %#ok<NBRAK>
        presentinator_error;
        protocol(bmplines,3) = -3;
    else % open and set to bgcolor
        try
            if ~polychrome.open
                [errCode polychrome.pTill] = calllib('TILLPolychrome','TILLPolychrome_Open',polychrome.pTill,0);
            else
                errCode = 0;
            end
            if ~errCode
                polychrome.open = true;
                errCode = calllib('TILLPolychrome','TILLPolychrome_SetRestingWavelength',polychrome.pTill,polychrome.color);
            end
        catch %#ok<*CTCH>
            txt = 'TILL Polychrome crashed';
            presentinator_error;
            protocol(bmplines,3) = -3;
        end
        if errCode
            txt = ['TILL Polychrome errorcode: ' num2str(errCode)];
            presentinator_error;
            protocol(bmplines,3) = -3;
        end
    end
end

%% convert color
% bmplines = find(protocol(:,1) == 0 & protocol(:,3) == 2);
% for i=1:length(bmplines)
%     protocol(bmplines(i),4) = clut(max(min(round(protocol(bmplines(i),4))+1,length(clut)),1));
% end

%% load imgs
if (verbose > 1)
    Screen('FillRect',w,bgcolor);
    Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
%     Screen(w,'DrawText',' Load bitmaps', 500, 15, infocolor);
    Screen('Flip',w);
end


if ~strcmp(pathstr_old,bmppathstr) && ~isempty(pathstr_old)
    % read canvas and sort
    pngcanvas = dir(fullfile(bmppathstr , '*.png'));
    pngcanvas = sort({pngcanvas.name});
    pngcanvas(find(cellfun('isempty',regexp(pngcanvas,'background.png'))))=[];
    % to clear old canvas; the first one is the background image, should not be
    % cleared
    scrCanvas = [scrCanvas(1) zeros(1,size(pngcanvas,2))];
%     currentCanvasIndex = 0; % temporarily: show zero background; in future - should be done through init file before (when no init.txt exist, make "0" as default)
end


%% if new directory
if ~strcmp(pathstr_old,bmppathstr) || isempty(scrOff)
    % read images and sort filenames
    bmpfiles = dir(fullfile(bmppathstr , '*.bmp'));
    bmpfiles = sort({bmpfiles.name});
    bmpfiles(find(~cellfun('isempty',regexp(bmpfiles,'background.bmp'))))=[];
    pngfiles = dir(fullfile(bmppathstr , '*.png'));
    pngfiles = sort({pngfiles.name});
    pngfiles(find(~cellfun('isempty',regexp(pngfiles,'background.png'))))=[];
    bmpfiles = [bmpfiles pngfiles];
    
    % clear old bitmaps
    pathstr_old = bmppathstr;
    scrOff = zeros(1,size(bmpfiles,2));
end

usedbmps = unique(protocol(find(protocol(2:end,1) > 0)+1,1))';
for i=usedbmps
    if i < length(scrOff) && scrOff(i)    % already loaded
        continue;
    end
    if i > length(bmpfiles)
        txt = ['Missing bmp id ' num2str(i-1) ' in ' bmppathstr];
        presentinator_error;
        protocol((protocol(:,1) == i),1) = 0;
        protocol((protocol(:,1) == i),3) = -3;
        continue;
    end

    [imageArray map alpha] = imread(fullfile(bmppathstr,bmpfiles{i}));
    imageInfo = imfinfo(fullfile(bmppathstr,bmpfiles{i}));
    if size(map,1) > 0
        if (isfield(imageInfo,'SimpleTransparencyData'))
            tp = logical(1-imageInfo.SimpleTransparencyData);
            map(tp,:) = repmat(bgcolor/255, size(map(tp,1)), 3);
        end
        imageArray = ind2gray(imageArray, map); % resolve indexed file
    end
    if ~isa(imageArray,'uint8')
        imageArray = (double(imageArray)+1)/256-1;
    end
    if ~isempty(alpha)
        imageArray(:,:,4) = alpha(:,:);
    end
    imageArray = uint8(imageArray);

    scrOff(i) = Screen('MakeTexture',w,imageArray);
    tmpIA = floor(size(imageArray)/2);
    sizeImage(i,1:2) = tmpIA(1,1:2);  

    %% Dont forget to make a correction for zooming later 
    
end

%% read canvas (background images)
usedcanvas = unique(protocol((find(protocol(2:end,3) == 6 & protocol(2:end,1) == 0)+1),4))';
for i=usedcanvas
    if i < length(scrCanvas) && scrCanvas(i+1)    % already loaded
        continue;
    end
    if i > length(pngcanvas)
        txt = ['Missing canvas id ' num2str(i-1) ' in ' bmppathstr];
        presentinator_error;
        protocol((protocol(:,1) == i),1) = 0;
        protocol((protocol(:,1) == i),3) = -3;
        continue;
    end

    [canvasArray map canvasAlpha] = imread(fullfile(bmppathstr,pngcanvas{i}));
    if isa(canvasArray, 'uint16')
        canvasArray(:,:,1:3) = (double(canvasArray(:,:,1:3))+1)/256-1;
    elseif ~isa(canvasArray, 'uint8')
        txt = ['Background file ',pngcanvas{i} ,' is neither UINT8 nor UINT16'];
        presentinator_error;
    end
    if ~isempty(canvasAlpha)
        canvasArray(:,:,4) = canvasAlpha(:,:);
    end
    canvasArray = uint8(canvasArray);    
    scrCanvas(i+1) = Screen('MakeTexture',w,canvasArray);
end

%% read checkers info
% checkers_presented=0;
checkers_protocol=protocol((find(protocol(2:end,3) == 7 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_protocol)
        imageSize=600;
        checkerSize=20;
        color1=0;
        color2=50;
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        time2pres=zeros(1,length(checkers_protocol));
        scr=zeros(2,length(checkers_protocol));
        % time of presentation in frames
        pres_cycles=protocol((find(protocol(2:end,3) == 7 & protocol(2:end,1) == 0)+1),2);        
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 7 & protocol(2:end,1) == 0)+1),2)=1;
    for checkers_presented=1:length(checkers_protocol)
        load(fullfile(bmppathstr,['checkers',int2str(checkers_protocol(checkers_presented))]));
        checkers=checkers+1;
        % time of presentation in s. substract 2 frames to be on the safe side (against time error)
        time2pres(checkers_presented)=pres_cycles(checkers_presented)*frametime-2/60;
        % create templates textures
        imageArray=zeros(checkerSize,checkerSize)+color1;
        imageArray2=zeros(checkerSize,checkerSize)+color2;
        imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
        imageArray2=uint8(cat(3,imageArray2,imageArray2,imageArray2));
        scr(1,checkers_presented)=Screen('MakeTexture',w,imageArray);
        scr(2,checkers_presented)=Screen('MakeTexture',w,imageArray2);

    end
    FlipTimestamp=zeros(max(pres_cycles),checkers_presented);
end
checkers_presented=0;


% 8 - high-low contrast
checkers_hl_protocol=protocol((find(protocol(2:end,3) == 8 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_hl_protocol)
        imageSize=600;
        checkerSize=20;
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        scr=zeros(1,255);
        % make texture of color from 0 to 255
        for i=0:255
            imageArray=zeros(checkerSize,checkerSize)+i;
            imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
            scr(i+1)=Screen('MakeTexture',w,imageArray);
        end
        corr_factor=0.2; % multiply sigma by this value
        sigma=9;
        mu=30;
        sigma_corr_factor=round(255*corr_factor);
        mu_corr_factor=(mu-mu*corr_factor)/255;
        % time of presentation in frames
        pres_cycles=12600-2;        
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 8 & protocol(2:end,1) == 0)+1),2)=1;
        seq=protocol(find(protocol(2:end,3) == 8,1),4);
        load(fullfile(bmppathstr,['HC_seq',int2str(seq+1)]));
        FlipTimestamp=zeros(max(pres_cycles),1);
end
checkers_hl_presented=0;

% 9 - checkers: only high contrast
checkers_hl_protocol=protocol((find(protocol(2:end,3) == 9 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_hl_protocol)
        imageSize=800;
        checkerSize=20;
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        scr=zeros(1,255);
        % make texture of color from 0 to 255
        for i=0:255
            imageArray=zeros(checkerSize,checkerSize)+i;
            imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
            scr(i+1)=Screen('MakeTexture',w,imageArray);
        end
        % time of presentation in frames
        pres_cycles=28800-2;        
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 8 & protocol(2:end,1) == 0)+1),2)=1;
        load(fullfile(bmppathstr,'HC_seq'));
        FlipTimestamp=zeros(max(pres_cycles),1);
end
checkers_h_presented=0;

% 10 - checkers: only high contrast, different sizes
checkers_s_protocol=protocol((find(protocol(2:end,3) == 10 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_s_protocol)
        imageSize=800;
        checkerSize=checkers_s_protocol/3; %in px
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        scr=zeros(1,61);
        % make texture of color from 0 to 255
        for i=0:60
            imageArray=zeros(checkerSize,checkerSize)+i;
            imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
            scr(i+1)=Screen('MakeTexture',w,imageArray);
        end
        % time of presentation in frames
        pres_cycles=18000-2;        
        % replace protocol duration entry with 1 frame as it is passed to rush
%         protocol((find(protocol(2:end,3) == 10 & protocol(2:end,1) == 0)+1),2)=1;
        load(fullfile(bmppathstr,['size_',int2str(checkers_s_protocol),'.mat']));
%         FlipTimestamp=zeros(max(pres_cycles),1);
end
checkers_s_presented=0;

% 11 - checkers: shuffle coord
checkers_h_protocol=protocol((find(protocol(2:end,3) == 11 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_h_protocol)
        imageSize=780;
        checkerSize=60;
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        scr=zeros(1,255);
        % make texture of color from 0 to 255
        for i=0:255
            imageArray=zeros(checkerSize,checkerSize)+i;
            imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
            scr(i+1)=Screen('MakeTexture',w,imageArray);
        end
        % time of presentation in frames
        pres_cycles=36000-2;        
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 8 & protocol(2:end,1) == 0)+1),2)=1;
        load(fullfile(bmppathstr,['seq',int2str(checkers_h_protocol)]));
        load(fullfile(bmppathstr,'coord_shift'));
        FlipTimestamp=zeros(max(pres_cycles),1);
end
checkers_h_presented=0;

% 12 - checkers: single cones
checkers_sg_protocol=protocol((find(protocol(2:end,3) == 12 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_sg_protocol)
        load(fullfile(bmppathstr,'s'));
        stim_duration=floor(protocol(find(protocol(:,3) == 12),2)/60);
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 12 & protocol(2:end,1) == 0)+1),2)=1;
end
checkers_sg_presented=0;

% 14 - checkers: binary
checkers_h_protocol=protocol((find(protocol(2:end,3) == 14 & protocol(2:end,1) == 0)+1),4)';
if ~isempty(checkers_h_protocol)
        imageSize=800;
        checkerSize=20;
        %create coordinates map
        steps=imageSize/checkerSize;
        center=screenRect(3:4)/2;
        coord_checkers=zeros(4,steps*steps);
        cnt_coord=1;
        for i_coord=1:steps
            for j_coord=1:steps
                coord_checkers(1,cnt_coord)=center(1)+(i_coord-steps/2-1)*checkerSize;
                coord_checkers(2,cnt_coord)=center(2)+(j_coord-steps/2-1)*checkerSize;
                cnt_coord=cnt_coord+1;
            end
        end
        coord_checkers(3,:)=coord_checkers(1,:)+checkerSize;
        coord_checkers(4,:)=coord_checkers(2,:)+checkerSize;
        
        scr=zeros(1,2);
        colors=[10,50];
        % make texture of color 10 and 50
        for i=0:1
            imageArray=zeros(checkerSize,checkerSize)+colors(i+1);
            imageArray=uint8(cat(3,imageArray,imageArray,imageArray));
            scr(i+1)=Screen('MakeTexture',w,imageArray);
        end
        % time of presentation in frames
        pres_cycles=18000-2;        
        % replace protocol duration entry with 1 frame as it is passed to rush
        protocol((find(protocol(2:end,3) == 8 & protocol(2:end,1) == 0)+1),2)=1;
        load(fullfile(bmppathstr,'HC_seq'));
        FlipTimestamp=zeros(max(pres_cycles),1);
end
checkers_h_presented=0;



if port % if port == 0 it was keyboard event
     txt = mat2str(protocol);
     pnet(udp,'write', sprintf(['STIM @ %s %s\n' strrep(txt(2:min(length(txt),40000)),';','\n')], ...
         datestr(now), protocolFileName));     
     pnet(udp,'writepacket',ip,port);   % Send buffer as UDP packet
end

% start 
triggeringError = false; % if false trigger error
timingTest = zeros(size(protocol,1),2);

% Rush(rushloop,2);
presentinator_rush; 


Screen('FillRect',w,bgcolor);
Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
Screen('Flip',w);

if triggeringError
    txt = ['No trigger arrived within ' num2str(triggeredStart(4)) 'ms']; %#ok<*UNRCH>
    presentinator_error;
else
    timingTest = timingTest / frametime;
    stepnum = size(timingTest,1);
    tmpTimingTest=timingTest;
    tmpProt=protocol;
    timingTest(:,1) = timingTest(:,1) - protocol(1:stepnum,2);

    errorlines = find(abs(timingTest(2:end,1))>= 1,5,'first');
    errorlines = [errorlines+1; find(abs(timingTest(:,2)) > 1,5,'first')];
    errorlines = sort(unique(errorlines));
    % to prevent stop if timing error
%     if ~isempty(errorlines)
%         txt = ['Timing error! ' ...
%             num2str(reshape(timingTest(errorlines,:)',1, length(errorlines)*2)) 'frames'];
%         presentinator_error;
%     end
    longestProcessing = ceil(100*max(timingTest(:,2)));
end

if strcmp(lastCommandScreen(1:5),'ERROR')
    errorLine = lastCommandScreen;
else
    errorLine = [];
end
% txt = mat2str(protocol);
% txt = sprintf(['s %s\n' strrep(txt(2:min(length(txt),40000)),';','\n')],protocolFileName);
% presentinator_log;

if ~isempty(errorLine)
    lastCommandScreen = errorLine;
end

protocol = [];
command = ' '; % delete command buffer
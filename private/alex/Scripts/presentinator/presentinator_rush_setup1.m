% set-up commands
do_checker=0;
do_checker_SG=0;
do_checker_HL=0;
do_checker_H=0;
do_checker_S=0;
do_checker_Sh=0;
contrast_sine=0;
do_checker_bin=0;
protocolLength = size(protocol,1); % used many times
% array of actions, {1,1} = number of elements
% {:,1} p - parallel, s - serial  {:,2} - port  {:,3} - command
data_cmd = repmat({0, '', ''},20,1); % pre-allocate space
flip = 0;       % flip square/signal HAS TO start with 0
keyIsDown = false; % to break free
pretime = GetSecs();
timingTest = zeros(protocolLength,2); % preallocate timing test with zeros
if triggeredStart(1) && triggeredStart(2) % active trigger
    if parallel_codes(6)
        parallel_states(parallel_codes(6)) = 1; % trigger
        parallel_out(parallel_codes(1),bin2dec(num2str(parallel_states)));
    end
    
    if verbose
        Screen('FillRect',w,bgcolor);
        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
        Screen(w,'DrawText',' Wait for trigger! ', 25, 30, infocolor);
        Screen('Flip',w);
    end
    parallel_out(890,uint32(bitset(uint8(parallel_in(890)),6))); % set to read
    p = parallel_in(triggeredStart(1));
    %disp(num2str([bitand(2^triggeredStart(3),p) p uint8(triggeredStart(2)) triggeredStart(3)  ]));
    while ~bitand(2^(triggeredStart(3)-1),p) && (triggeredStart(4)/1000 > GetSecs() - pretime)
        p = uint8(parallel_in(triggeredStart(1)));
    end
    parallel_out(890,uint32(bitset(uint8(parallel_in(890)),6,0))); % set to write
    if parallel_codes(6)
        parallel_states(parallel_codes(6)) = 0; % trigger received
    end
    
    if ~bitand(2^(triggeredStart(3)-1),p) % not triggered but time-out
        triggeringError = true;
    end
else
    for waitsync=1:waitFrames(3)   % give some time to the sender
        Screen(w,'WaitBlanking');
    end
    if parallel_codes(4) % first/last
        parallel_states(parallel_codes(4)) = 1;
    end
end
abused=protocol;

if ~triggeringError
    wholeFlipTime=tic; % counts whole time elapsed between the flips of the screen
    pretime=tic; % counts time elapsed for all preparations, but not waiting time and flips themselves
    % present - procItem=1 for init first
    procItem = 1; waitItem = 1;
    while procItem < protocolLength+1 % wait for poststim time
        infopos = 15;
        data_cmd{1,1} = 1; % valid actions in this array (1 = no action)
        doit = (procItem <= protocolLength); % if last line don't comp just wait
        while doit % do until duration not zero
            % generate next frame
            if verbose > 2 % debug
                infocolor = txtcolor;
            else
                if bgcolor < 160
                    infocolor = bgcolor+40;
                else
                    infocolor = bgcolor-70;
                end
            end
            if protocol(procItem,1) > 0 % bmp
                destRect = [screenRect(3)/2-sizeImage(protocol(procItem,1),2)-1-protocol(procItem,3)...
                    screenRect(4)/2-sizeImage(protocol(procItem,1),1)-1-protocol(procItem,4)...
                    screenRect(3)/2+sizeImage(protocol(procItem,1),2)-protocol(procItem,3)...
                    screenRect(4)/2+sizeImage(protocol(procItem,1),1)-protocol(procItem,4)];
                Screen('FillRect',w,bgcolor);
                Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
                Screen('DrawTexture',w,scrOff(protocol(procItem,1)),[],destRect);
            elseif protocol(procItem,1) == 0 % bg -1
                switch protocol(procItem,3)
                    case 0 % bg color
                        Screen('FillRect',w,bgcolor);
                        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
                    case 2 % full-field
                        bgcolor=uint8(protocol(procItem,4));
                        Screen('FillRect',w,bgcolor);
                        if isempty(polyMask) && verbose < 3
                            if protocol(procItem,4) < 160
                                infocolor = protocol(procItem,4)+40;
                            else
                                infocolor = protocol(procItem,4)-70;
                            end
                        end
                    case 4 % user
                        Screen('FillRect',w,bgcolor);
                        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
                        cmd = user_commands{protocol(procItem,4)};
                        if verbose > 1 && protocol(procItem,2)
                            Screen(scrOffline,'DrawText', [cmd{4} ' (' cmd{2} ' ' cmd{3} ')'], 25, infopos, infocolor);infopos=infopos+15;
                        end
                        
                        data_cmd{1,1} = data_cmd{1,1} + 1;
                        data_cmd{data_cmd{1,1},1} = cmd{5};
                        data_cmd{data_cmd{1,1},2} = cmd{6};
                        data_cmd{data_cmd{1,1},3} = cmd{7};
                        
                        %% time 0.04 / 24
                        switch cmd{5}(1)
                            case 'p'
                                if cmd{6} == parallel_codes(1) && abs(cmd{3}) < 9
                                    data_cmd{1,1} = data_cmd{1,1} - 1; % dont store
                                    % same parallel port, so cmd{3} is a pin:
                                    if cmd{3} < 0 % -1 0 4 -2 = pin2 is low
                                        parallel_codes(9+cmd{3}) = 0;
                                    else          % -1 0 4 2  = pin2 is high
                                        parallel_codes(9-cmd{3}) = 1;
                                    end
                                end
                            case 's'
                                so = so_arr{str2double(cmd{2}(4))};
                                data_cmd{data_cmd{1,1},2} = so;
                                if ~isempty(so) && so.status(1) == 'o'
                                    if so.BytesAvailable
                                        fread(so, so.BytesAvailable, 'char');
                                    end
                                end
                            case 'x'
                                if verbose > 1
                                    %                             Screen(scrOffline,'DrawText', [' ' cmd{7}], 25, infopos, infocolor);infopos=infopos+15;
                                    Screen(w,'DrawText', [' ' cmd{7}], 25, infopos, infocolor);infopos=infopos+15;
                                end
                        end
                        
                    case 5 % till polychrome
                        Screen('FillRect',w,bgcolor);
                        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
                        if verbose > 1
                            Screen(w,'DrawText', [' Till Polychrome to ' num2str(protocol(procItem,4)) 'nm'], 25, infopos, infocolor);infopos=infopos+15;
                        end
                        data_cmd{1,1} = data_cmd{1,1} + 1;
                        data_cmd{data_cmd{1,1},1} = 't';
                        data_cmd{data_cmd{1,1},2} = protocol(procItem,4); % color
                    case 6
                        currentCanvasIndex=protocol(procItem,4);
                        Screen('FillRect',w,bgcolor);
                        Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
                    case 7 % pseudo random flicker
                        do_checker=1;
                        checkers_presented=checkers_presented+1;
                    case 8 % pseudo random flicker High Low contrast
                        do_checker_HL=1;
                    case 9 % pseudo random checkerboard flicker High contrast
                        do_checker_H=1;
                    case 10 % pseudo random flicker dif. sizes
                        do_checker_S=1;
                    case 11 % pseudo random flicker shuffle coord
                        do_checker_Sh=1;
                    case 12 % pseudo random flicker single cones
                        do_checker_SG=1;
                    case 13 % contrast sine
                        contrast_sine=1;
                    case 14 % binary checkerboard 40x40
                        do_checker_bin=1;

                end
            end
            if ~isempty(polyMask)
                Screen('FillPoly',w,bgcolor,polyMask);
            end
            
            if verbose == 3 % debug
                 flip = ~flip;
            end
            
            if protocol(procItem,2) || procItem == protocolLength
                doit = false;
                if verbose > 1 % write only the last
                    %% time 0.02 / 36
                    txt = [protocolFileName ' step ' num2str(procItem-1) 'of' num2str(protocolLength) ' remains ' num2str(sum(protocol(procItem:end,2))*frametime) 's'];
                    % possible bug %                 Screen(scrOffline,'DrawText',txt, 25, screenRect(4), infocolor);
                end
            else % if duration == 0 do the next line
                abused(procItem,5)=0;
                procItem = procItem + 1;
            end
        end % doit while duration not zero
        
        if parallel_codes(5) % show action after next Blanking when its on screen
            parallel_states(parallel_codes(5)) = waitItem < protocolLength; % posttime not signaled
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
        end

        timingTest(procItem,2) = toc(pretime);
 
        if parallel_codes(3) %first frame on screen
            parallel_states(parallel_codes(3)) = (procItem <= protocolLength && protocol(procItem,2)); % stim
        end
        if parallel_codes(5) % show action after next Blanking when its on screen
            parallel_states(parallel_codes(5)) = (procItem == 1); % pretime not signaled
        end
        
        % if protocol is done, show just canvas
        if procItem>protocolLength
            Screen('FillRect',w,bgcolor);
            Screen('DrawTexture',w,scrCanvas(currentCanvasIndex+1));
        end
        
        % wait the req. time
        while toc(wholeFlipTime)<(protocol(waitItem,2))*frametime
            if keyIsDown
                break;             %#ok<UNRCH>
            end
        end
        timingTest(waitItem,1)=toc(wholeFlipTime);
        wholeFlipTime=tic;
        
        if do_checker_bin
            
            Screen('FillRect',w,uint8(30));
            not_needed=0;
            FlipTimestamp=zeros(1,18500);  % for max 15000 flips (app.6 min with 60hZ)
            load('C:\stimfiles\_singleCone_\s')
            defaultStream = RandStream.getDefaultStream;
            defaultStream.State = s(1).State;
            b=uint8(randi([0 1],1,1600)); 
            Screen('DrawTextures',w,scr(b+1),[],coord_checkers);
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            
            tt=FlipTimestamp(1);
            kkk=2;
            while kkk<=18000 && FlipTimestamp(kkk-1)-tt<300
                b=uint8(randi([0 1],1,1600));
                Screen('DrawTextures',w,scr(b+1),[],coord_checkers)
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w,0,1,0);
                kkk=kkk+1;
            end
            FlipTimestamp=FlipTimestamp(1:kkk-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_bin=0;
            checkers_h_presented=1;
        end

        
        if do_checker
            not_needed=0;
            Screen('DrawTextures',w,scr(checkers(1,:),checkers_presented),[],coord_checkers);
            [ttt not_needed FlipTimestamp(1,checkers_presented)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            iii=2;
            while FlipTimestamp(iii-1,checkers_presented)-FlipTimestamp(1,checkers_presented)<time2pres(checkers_presented)
                Screen('DrawTextures',w,scr(checkers(iii,:),checkers_presented),[],coord_checkers);
                [not_needed not_needed FlipTimestamp(iii,checkers_presented)]=Screen('Flip',w,0,1,0);
                iii=iii+1;
            end
            FlipTimestamp(1:iii-1,checkers_presented)=FlipTimestamp(1:iii-1,checkers_presented)-ttt;
            do_checker=0;
        end
        
        if do_checker_HL
            Screen('FillRect',w,uint8(30));
            not_needed=0;            
            FlipTimestamp=zeros(1,25300);            
            glsl1 = MakeTextureDrawShader(w, 'SeparateAlphaChannel',[mu_corr_factor*ones(1,3),0]);
            Screen('DrawTextures',w,scr(high_contr(1,:)+1),[],coord_checkers);
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            
            tt=FlipTimestamp(1);
            iii=2;kkk=2;
            for j=1:8
                if mod(j,2)==1
                    while kkk<901 && FlipTimestamp(iii-1)-tt<15
                        Screen('DrawTextures',w,scr(high_contr(kkk,:)+1),[],coord_checkers)
                        [not_needed not_needed FlipTimestamp(iii)]=Screen('Flip',w,0,1,0);
                        iii=iii+1;
                        kkk=kkk+1;
                    end
                    tt=FlipTimestamp(iii-1);
                    kkk=1;
                else
                    while kkk<901 && FlipTimestamp(iii-1)-tt<15
                        Screen('DrawTextures',w,scr(high_contr(kkk,:)+1),[],coord_checkers,[],[],[],sigma_corr_factor*ones(1,3),glsl1);
                        [not_needed not_needed FlipTimestamp(iii)]=Screen('Flip',w,0,1,0);
                        iii=iii+1;
                        kkk=kkk+1;
                    end
                    kkk=1;
                    tt=FlipTimestamp(iii-1);
                end
            end
            FlipTimestamp=FlipTimestamp(1:iii-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_HL=0;
            checkers_hl_presented=1;
        end
      
        if do_checker_H
            Screen('FillRect',w,uint8(30));
            not_needed=0;            
            FlipTimestamp=zeros(1,35000);            
            Screen('DrawTextures',w,scr(high_contr(1,:)+1),[],coord_checkers);
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            
            tt=FlipTimestamp(1);
            kkk=2;
            while kkk<=28800 && FlipTimestamp(kkk-1)-tt<480
                Screen('DrawTextures',w,scr(high_contr(kkk,:)+1),[],coord_checkers)
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w,0,1,0);
                kkk=kkk+1;
            end
            FlipTimestamp=FlipTimestamp(1:kkk-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_H=0;
            checkers_h_presented=1;
        end

        if do_checker_S
            Screen('FillRect',w,uint8(30));
            not_needed=0;            
            FlipTimestamp=zeros(1,22000);            
            Screen('DrawTextures',w,scr(high_contr(1,:)+1),[],coord_checkers);
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            
            tt=FlipTimestamp(1);
            kkk=2;
            while kkk<=18000 && FlipTimestamp(kkk-1)-tt<300
                Screen('DrawTextures',w,scr(high_contr(kkk,:)+1),[],coord_checkers)
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w,0,1,0);
                kkk=kkk+1;
            end
            FlipTimestamp=FlipTimestamp(1:kkk-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_S=0;
            checkers_s_presented=1;
        end

        if do_checker_Sh
            Screen('FillRect',w,uint8(30));
            not_needed=0;
            FlipTimestamp=zeros(1,50000);
            coor=coord_checkers+repmat(coord_shift(:,1),2,size(coord_checkers,2));
            Screen('DrawTextures',w,scr(high_contr(1,:)+1),[],coor);
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            
            tt=FlipTimestamp(1);
            kkk=2;kkksh=2;
            while kkk<=36000 && FlipTimestamp(kkk-1)-tt<600
                coor=coord_checkers+repmat(coord_shift(:,kkksh),2,size(coord_checkers,2));            
                Screen('DrawTextures',w,scr(high_contr(kkk,:)+1),[],coor)
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w,0,1,0);
                kkk=kkk+1;
                kkksh=kkksh+1;
                if kkksh>size(coord_shift,2)
                    kkksh=1;
                end
            end
            FlipTimestamp=FlipTimestamp(1:kkk-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_Sh=0;
            checkers_h_presented=1;
        end
        
        if do_checker_SG
            Screen('FillRect',w,uint8(30));
            not_needed=0;            
            FlipTimestamp=zeros(1,15000);  % for max 15000 flips (app.6 min with 60hZ)
            defaultStream = RandStream.getDefaultStream;
            defaultStream.State = s(checkers_sg_protocol).State;
            imageArray=zeros(800,800);
            b=uint8(randi([0 1],1,160000)*60); % colors 0 and 60
            b=repmat(b,2,1);
            b=reshape(b,320000,1);
            b=reshape(b,800,400);
            imageArray(:,1:2:end)=b;
            imageArray(:,2:2:end)=b;
            scr=Screen('MakeTexture',w,imageArray);
            Screen('DrawTexture',w,scr,[])
            [ttt not_needed FlipTimestamp(1)]=Screen('Flip',w,0,1,0);
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            tt=FlipTimestamp(1);
            kkk=2;
            while FlipTimestamp(kkk-1)-tt<stim_duration %max time in s
                % make texture of color 0 and 60
                imageArray=zeros(800,800);
                b=uint8(randi([0 1],1,160000)*60); % colors 0 and 60
                b=repmat(b,2,1);
                b=reshape(b,320000,1);
                b=reshape(b,800,400);
                imageArray(:,1:2:end)=b;
                imageArray(:,2:2:end)=b;
                Screen('Close',scr);
                scr=Screen('MakeTexture',w,imageArray);
                Screen('DrawTexture',w,scr,[])
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w,0,1,0);
                kkk=kkk+1;
            end
            FlipTimestamp=FlipTimestamp(1:kkk-1);
            FlipTimestamp=FlipTimestamp-ttt;
            do_checker_SG=0;
            checkers_sg_presented=1;
        end
        
        if contrast_sine
            load('C:\stimfiles\_hlc2_\sineContr.mat') 
            mySine=uint8(mySine);
            FlipTimestamp=zeros(1,length(mySine));
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
            kkk=1;
            while kkk<=36000 && FlipTimestamp(kkk)-FlipTimestamp(1)<=600000
                Screen('FillRect',w,mySine(kkk));
                [not_needed not_needed FlipTimestamp(kkk)]=Screen('Flip',w);
                kkk=kkk+1;
            end

        end
        % Copy to monitor
        [not_needed not_needed currentFlip]=Screen('Flip',w);
        abused(procItem,5)=currentFlip;
%         Screen('Flip',w);
        
%         posttime = GetSecs();
        if parallel_codes(1) % show action after next Blanking when its on screen
            parallel_out(parallel_codes(1),2*(2*(2*(2*(2*(2*(2*parallel_states(1)+parallel_states(2))+parallel_states(3))+parallel_states(4))+parallel_states(5))+parallel_states(6))+parallel_states(7))+parallel_states(8));
        end
        if parallel_codes(4)
            parallel_states(parallel_codes(4)) = 0;
            %%% code for making parallel_codes(4) work as a vsync indicator
            %         parallel_states(parallel_codes(4)) = 1-parallel_states(parallel_codes(4));
            %%%            
        end
        
        pretime=tic;
        
        waitItem = procItem;
        % Do the data communications
        for dataItem=2:data_cmd{1,1} % first line is the number of actions
            switch data_cmd{dataItem,1}
                case 'p'
                    parallel_out(data_cmd{dataItem,2}, data_cmd{dataItem,3});%for parallel commands, action is stored in 3, port in 2
                case 's'
                    while ~strcmp('idle',data_cmd{dataItem,2}.TransferStatus)
                    end
                    fprintf(data_cmd{dataItem,2}, data_cmd{dataItem,3},'async');%for serial commands, action is stored in 3, port in 2
                case 'sOlympus'
                    while ~strcmp('idle',data_cmd{dataItem,2}.TransferStatus)
                    end
                    fwrite(data_cmd{dataItem,2}, data_cmd{dataItem,3},'int8');%for serial commands, action is stored in 3, port in 2
                case 'x'
                    eval(strcat('h_Sol.',data_cmd{dataItem,2}));%for ActX commands, action is stored in 2, comment in 3
                case 't'
                    calllib('TILLPolychrome','TILLPolychrome_SetRestingWavelength',polychrome.pTill,data_cmd{dataItem,2});
            end
        end

        procItem = procItem + 1;
        % check keyboard activity
        if keyIsDown
            break; %#ok<UNRCH>
        end
    end

end % no trigger error

% set all signals to zero
if parallel_codes(4)
    parallel_states(parallel_codes(4)) = 1;
    parallel_out(parallel_codes(1),bin2dec(num2str(parallel_states)));
    Screen(w,'WaitBlanking');
end
parallel_states = false(1,8);
if parallel_codes(1)
    parallel_out(parallel_codes(1),bin2dec(num2str(parallel_states)));
end
Screen(w,'WaitBlanking');
parallel_states = true(1,8);
if parallel_codes(1)
    parallel_out(parallel_codes(1),bin2dec(num2str(parallel_states)));
end
Screen(w,'WaitBlanking');
parallel_states = false(1,8);
if parallel_codes(1)
    parallel_out(parallel_codes(1),bin2dec(num2str(parallel_states)));
end

timingTest = timingTest(1:end-1,:); % last line is after postrec
if checkers_presented
    if ~isdir(['D:\Users\alexandra\checkers\',date_tmp])
        mkdir('D:\Users\alexandra\checkers\',date_tmp)
    end
    for aaa=1:checkers_presented
        heka_file_name=['000' int2str(heka_file_ID)];
        FlipTimeStamps=FlipTimestamp(:,aaa);
        if ~isempty(find(FlipTimeStamps==0,1))
            FlipTimeStamps=FlipTimeStamps(1:find(FlipTimeStamps==0,1)-1);
        end            
        save(['D:\Users\alexandra\checkers\',date_tmp,'\',date_tmp,'_#',heka_file_name(end-3:end),'_checkerSequ_',int2str(checkers_protocol(aaa)),'_',int2str(aaa),'.mat'],'FlipTimeStamps');
    end
    checkers_presented=0;
end

if checkers_hl_presented
    if ~isdir(['D:\Users\alexandra\checkers_hl\',date_tmp])
        mkdir('D:\Users\alexandra\checkers_hl\',date_tmp)
    end
    heka_file_name=['000' int2str(heka_file_ID)];
    FlipTimeStamps=FlipTimestamp;
    if ~isempty(find(FlipTimeStamps==0,1))
        FlipTimeStamps=FlipTimeStamps(1:find(FlipTimeStamps==0,1)-1);
    end
    save(['D:\Users\alexandra\checkers_hl\',date_tmp,'\',date_tmp,'_#',heka_file_name(end-3:end),'_checker_hl_seq_',int2str(checkers_hl_protocol),'.mat'],'FlipTimeStamps');
    checkers_hl_presented=0;
end


if checkers_h_presented
    if ~isdir(['D:\Users\alexandra\checkers_h\',date_tmp])
        mkdir('D:\Users\alexandra\checkers_h\',date_tmp)
    end
    heka_file_name=['000' int2str(heka_file_ID)];
    FlipTimeStamps=FlipTimestamp;
    if ~isempty(find(FlipTimeStamps==0,1))
        FlipTimeStamps=FlipTimeStamps(1:find(FlipTimeStamps==0,1)-1);
    end
    save(['D:\Users\alexandra\checkers_h\',date_tmp,'\',date_tmp,'_#',heka_file_name(end-3:end),'_checker_h_seq_',int2str(checkers_hl_protocol),'.mat'],'FlipTimeStamps');
    checkers_h_presented=0;
end

if checkers_s_presented
    if ~isdir(['D:\Users\alexandra\checkers_s\',date_tmp])
        mkdir('D:\Users\alexandra\checkers_s\',date_tmp)
    end
    heka_file_name=['000' int2str(heka_file_ID)];
    FlipTimeStamps=FlipTimestamp;
    if ~isempty(find(FlipTimeStamps==0,1))
        FlipTimeStamps=FlipTimeStamps(1:find(FlipTimeStamps==0,1)-1);
    end
    save(['D:\Users\alexandra\checkers_s\',date_tmp,'\',date_tmp,'_#',heka_file_name(end-3:end),'_checker_size_',int2str(checkers_s_protocol),'.mat'],'FlipTimeStamps');
    checkers_s_presented=0;
end

if checkers_sg_presented
    if ~isdir(['D:\Users\alexandra\checkers_sg\',date_tmp])
        mkdir('D:\Users\alexandra\checkers_sg\',date_tmp)
    end
    heka_file_name=['000' int2str(heka_file_ID)];
    FlipTimeStamps=FlipTimestamp;
    if ~isempty(find(FlipTimeStamps==0,1))
        FlipTimeStamps=FlipTimeStamps(1:find(FlipTimeStamps==0,1)-1);
    end
    save(['D:\Users\alexandra\checkers_sg\',date_tmp,'\',date_tmp,'_#',heka_file_name(end-3:end),'_checker_sg_seq_',int2str(checkers_sg_protocol),'.mat'],'FlipTimeStamps');
    checkers_sg_presented=0;
end


if ~isdir(['D:\matlabProtocols\',date_tmp])
    mkdir(['D:\matlabProtocols\',date_tmp])    
end
heka_file_name=['000' int2str(heka_file_ID)];
fileName=['D:\matlabProtocols\',date_tmp,'\',heka_file_name(end-3:end),'.mat'];
cnt=1;
while exist(fileName,'file')
    if cnt==1
        fileName=[fileName(1:end-4),'_',int2str(cnt),'.mat'];
    else
        fileName=[fileName(1:end-5),int2str(cnt),'.mat'];
    end
    cnt=cnt+1;
end
save(fileName,'abused');

if contrast_sine
    heka_file_name=['000' int2str(heka_file_ID)];
    fileName=['D:\matlabProtocols\',date_tmp,'\',heka_file_name(end-3:end),'.mat'];
    cnt=1;
    while exist(fileName,'file')
        if cnt==1
            fileName=[fileName(1:end-4),'_',int2str(cnt),'.mat'];
        else
            fileName=[fileName(1:end-5),int2str(cnt),'.mat'];
        end
        cnt=cnt+1;
    end
    save([fileName(1:end-4),'_sine.mat'],'FlipTimestamp');
    contrast_sine=0;
end
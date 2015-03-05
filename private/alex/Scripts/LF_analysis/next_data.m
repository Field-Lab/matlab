function next_data(flag,LinearFilter,SpikeCount,trialsPerScreen,handles,names,polarity)

global fitParams currentUnit start_trial currentLFs 


for i=1:size(fitParams.subplots,1)
    if ishandle(fitParams.subplots(i,1))
        delete(fitParams.subplots(i,1))
    end
end

if flag==0 % next trials

    start_trial=start_trial+trialsPerScreen;

    if start_trial>size(LinearFilter,2)
        next_data(1,LinearFilter,SpikeCount,trialsPerScreen,handles,names,polarity)
    end
        
elseif flag ==1 % next unit
    
    currentUnit=currentUnit+1;
    if currentUnit>size(LinearFilter,3)
        display('GOING ROUND! Unit 1')
        currentUnit=1;
    end        
    start_trial=1;
    currentLFs=LinearFilter(:,:,currentUnit);
    set(handles.on,'Value',polarity(currentUnit)/2+0.5);
    set(handles.off,'Value',~(polarity(currentUnit)/2+0.5));
    calcParams(currentLFs,SpikeCount(:,currentUnit),handles,trialsPerScreen,start_trial)
    
else % jump to a different unit

    currentUnit=str2num(get(handles.enterGotoUnit,'String'));
    if currentUnit>size(LinearFilter,3)
        display('You don''t have so many units! went to unit 1')
        currentUnit=1;
    end
    start_trial=1;
    currentLFs=LinearFilter(:,:,currentUnit);
    set(handles.on,'Value',polarity(currentUnit)/2+0.5);
    set(handles.off,'Value',~(polarity(currentUnit)/2+0.5));
    calcParams(currentLFs,SpikeCount(:,currentUnit),handles,trialsPerScreen,start_trial)
end




makePlots(currentLFs,trialsPerScreen,start_trial,[],1);

set(handles.info,'string', [names{currentUnit}, '  (',int2str(currentUnit),' out of ',...
    int2str(size(LinearFilter,3)),')   Trials ',int2str(start_trial),' to ',int2str(min(start_trial+trialsPerScreen-1,size(LinearFilter,2))),' out of ',int2str(size(LinearFilter,2))]);
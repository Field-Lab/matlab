function set_thr(flag)

global spikes thr hist_fig ch
persistent line_id

if flag==1 % first creation: automatic threshold
    h=findobj('Name','Spike Minima Histograms');
    if ~isempty(h)
        close(h)
    end
    [vals,ranges]=hist(min(spikes),500);
    [max_val,pos_max]=max(vals);
    thr=ranges(pos_max-find(vals(pos_max:-1:1)<max_val/5,1,'first'));
    hist_fig=figure;
    set(hist_fig,'Name','Spike Minima Histograms','Units', 'normalized','Position',[0.002 0.5171 0.2 0.4075])
elseif flag==0 % reset threshold
    hist_fig=findobj('Name','Spike Minima Histograms');
    if isempty(hist_fig)
        display('Hardly believable, huh?')
    end
    [thr,~]=ginput(1);
    delete(line_id)
elseif flag==3 % forel reset
    hist_fig=findobj('Name','Spike Minima Histograms');
    if isempty(hist_fig)
        display('Hardly believable, huh?')
    end
    delete(line_id)
end

figure(hist_fig)
subplot(2,1,1)
hist(min(spikes),250);
axis tight
line_id=line([thr,thr],get(gca,'YLim'),'color','r','lineWidth',2);
title(['CH',int2str(ch),' minima of all extracted wfs. N=',int2str(size(spikes,2))])
set(gca,'ytick',0,'yTickLabel','','position',[0.1    0.58   0.85    0.33])
subplot(2,1,2)
hist(min(spikes(:,min(spikes)<thr)),200);
axis tight
title(['minima of wfs for pca. N=',int2str(sum(min(spikes)<thr)),'. Threshold: ',int2str(thr)])
set(gca,'ytick',0,'yTickLabel','','position',[0.1   0.1    0.85    0.33])

% reset threshold UI
resetThrUI=uicontrol('style','pushbutton','string','reset threshold');
set(resetThrUI,'Units', 'normalized','Position',[0.72 0.485 0.25 0.05],'fontsize',12, 'Callback','set_thr(0)');

% delete spikes not reaching threshold. CAREFUL! No copy saved, for
% memory's sake. Reload everything if you want your spikes back.
% k='delete_spikes()';
deleteSpikesThrUI=uicontrol('style','pushbutton','string','remove spikes','BackgroundColor','r');
set(deleteSpikesThrUI,'Units', 'normalized','Position',[0.02 0.485 0.3 0.05],'fontsize',12, 'Callback','delete_spikes(thr)');

getFileUI=uicontrol('style','pushbutton','string','get new','BackgroundColor',[0.5 0.5 1]);
set(getFileUI,'Units', 'normalized','Position',[0.05 0.01 0.2 0.05],'fontsize',12,'fontweight','bold','Callback','get_file(1)');


getNextFileUI=uicontrol('style','pushbutton','string','get next','BackgroundColor',[0.5 0.5 1]);
set(getNextFileUI,'Units', 'normalized','Position',[0.45 0.01 0.22 0.05],'fontsize',12,'fontweight','bold','Callback','get_file(2)');

runUI=uicontrol('style','pushbutton','string','RUN','BackgroundColor',[0.2 0.8 0.2]);
set(runUI,'Units', 'normalized','Position',[0.85 0.01 0.1 0.05],'fontsize',12,'fontweight','bold','Callback','myGUI');

set(gcf,'Toolbar','figure')

end
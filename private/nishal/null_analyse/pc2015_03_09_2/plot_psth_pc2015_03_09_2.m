function [h] = plot_psth_pc2015_03_09_2(spkCondColl,Condtime,event_threshold,condDraw)

binSz=1/1200;
len=ceil(Condtime/binSz);
convolve=50;
h=figure('Color','w');
shift=0;
lb=0.2;
ub=0.9;
%event_tol=0.05; % accuracy of events in seconds
col='rkrkrkrkrkr';
for icond=condDraw(end:-1:1);
    
    rec_rast= makeSpikeMat(spkCondColl(icond).spksColl,binSz,len);
    spkCondColl(icond).rec_rast=rec_rast;
    
    [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,rec_rast);

    
    events = (PSTH_rec/max(PSTH_rec(time>lb*max(time) & time<ub*max(time))))>event_threshold & time>lb*max(time) & time<ub*max(time);
    
    colors= 0*spkCondColl(icond).xPoints;
    for ipts=1:length(spkCondColl(icond).xPoints)
        if(~isnan(spkCondColl(icond).xPoints(ipts)))
            if(floor((spkCondColl(icond).xPoints(ipts)/20000)/binSz)~=0)
    colors(ipts) = events(floor((spkCondColl(icond).xPoints(ipts)/20000)/binSz));
            end
        end
    end
    
    %    plot(spkCondColl(icond).xPoints/20000,spkCondColl(icond).yPoints,'r');
   plot(spkCondColl(icond).xPoints/20000,spkCondColl(icond).yPoints+shift,col(icond));
    hold on;
    PSTH_rec_plot =PSTH_rec*max(spkCondColl(icond).yPoints)/max(PSTH_rec);
    plot(time,PSTH_rec_plot+shift,col(icond));
    hold on;  
    spkCondColl(icond).events=events;
    shift=shift+max(spkCondColl(icond).yPoints);
   % shift=0;
    
    plot(time(logical(events)),(shift+1)*ones(sum(events),1),strcat(col(icond),'.'));
    shift=shift+1;
   % shift=0;
  
end
%xlim([0,max(time)]);
ylim([0,shift]);
xlim([4,10]);
set(gca,'ytick',[]);


end
function [spkCondColl,eventOverlap,h] = event_count(spkCondColl,Condtime,event_threshold, event_tol)

binSz=1/1200;
len=ceil(Condtime/binSz);
convolve=50;
h=figure('Color','w');
shift=0;
lb=0.2;
ub=0.9;
%event_tol=0.05; % accuracy of events in seconds
col='rkrkrkrkrkr';
for icond=length(spkCondColl):-1:1
    
    rec_rast= makeSpikeMat(spkCondColl{icond}.spksColl,binSz,len);
    spkCondColl{icond}.rec_rast=rec_rast;
    
    [PSTH_rec,time]=calculate_psth_fcn2(convolve,binSz,len,rec_rast);

    
    events = (PSTH_rec/max(PSTH_rec(time>lb*max(time) & time<ub*max(time))))>event_threshold & time>lb*max(time) & time<ub*max(time);
    
    colors= 0*spkCondColl{icond}.xPoints;
    for ipts=1:length(spkCondColl{icond}.xPoints)
        if(~isnan(spkCondColl{icond}.xPoints(ipts)))
            if(floor((spkCondColl{icond}.xPoints(ipts)/20000)/binSz)~=0)
    colors(ipts) = events(floor((spkCondColl{icond}.xPoints(ipts)/20000)/binSz));
            end
        end
    end
    
    %    plot(spkCondColl{icond}.xPoints/20000,spkCondColl{icond}.yPoints,'r');
    plot(spkCondColl{icond}.xPoints/20000,spkCondColl{icond}.yPoints+shift,col(icond));
    hold on;
    PSTH_rec_plot =PSTH_rec*max(spkCondColl{icond}.yPoints)/max(PSTH_rec);
    plot(time,PSTH_rec_plot+shift,'b');
    hold on;  
    spkCondColl{icond}.events=events;
    shift=shift+max(spkCondColl{icond}.yPoints);
    
    plot(time(logical(events)),(shift+1)*ones(sum(events),1),'g.');
    shift=shift+1;
 
  
end
xlim([0,max(time)]);
ylim([0,shift]);


%% make event start and end windows

for icond=1:length(spkCondColl)
   spkCondColl{icond}.events_list_start=[];
   spkCondColl{icond}.events_list_end=[];
   spkCondColl{icond}.events_list_mean=[];
   spkCondColl{icond}.events_list_length=[];
    rec_rast=spkCondColl{icond}.rec_rast;
    
    ievent=0;
    inside_event=0;
    for ibin=1:length(spkCondColl{icond}.events)
        
        if(spkCondColl{icond}.events(ibin)==1 & inside_event==0)
            ievent=ievent+1;
            inside_event=1;
           spkCondColl{icond}.events_list_start(ievent)=time(ibin);
           spkCondColl{icond}.events_bin_start(ievent)=ibin; 
        end
        
        if(spkCondColl{icond}.events(ibin)==0 & inside_event==1)
           spkCondColl{icond}.events_list_end(ievent)=time(ibin);
            spkCondColl{icond}.events_list_mean(ievent) =  (spkCondColl{icond}.events_list_end(ievent) +spkCondColl{icond}.events_list_start(ievent))/2;
            spkCondColl{icond}.events_list_length(ievent) = (spkCondColl{icond}.events_list_end(ievent)-spkCondColl{icond}.events_list_start(ievent));
            spkCondColl{icond}.events_spks(ievent) = sum(sum(rec_rast(:,time >= spkCondColl{icond}.events_list_start(ievent) & time<= spkCondColl{icond}.events_list_end(ievent))))/(spkCondColl{icond}.events_list_length(ievent)*size(rec_rast,1));
         spkCondColl{icond}.events_bin_end(ievent)=ibin; 
             inside_event=0;
        end
        
    end
end

%% Event overlap index
eventOverlap=zeros(length(spkCondColl),length(spkCondColl));
for icond1=1:length(spkCondColl)
    for icond2=1:length(spkCondColl)
        if(icond1~=icond2)
            for ievent1=1:length(spkCondColl{icond1}.events_list_mean)
                event_time=spkCondColl{icond1}.events_list_mean(ievent1);
                 for ievent2=1:length(spkCondColl{icond2}.events_list_mean)
                     if(abs(spkCondColl{icond2}.events_list_mean(ievent2) - event_time)<event_tol)
                         eventOverlap(icond1,icond2)= eventOverlap(icond1,icond2)+1;
                         break;
                     end
                 end
            end
        end
         eventOverlap(icond1,icond2)= eventOverlap(icond1,icond2)/length(spkCondColl{icond1}.events_list_mean);
    end
end


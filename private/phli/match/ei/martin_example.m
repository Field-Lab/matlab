if 0 %poster
    clf
    %cell_ids=[6181];ch=[413 399 430 190 174];
   
    cell_ids=[5283];ch=[353 372 326 240 216];
   
    cell_ids=[332];ch=[23 17 53 95 139];
   
    index=get_cell_indices(datarun{1},cell_ids);

    points=600;
    start=12300;
   
    points=800;
    start=9550;


    rawFile = edu.ucsc.neurobiology.vision.io.RawDataFile('/Volumes/Twist/Data/Greschner/2009-12-03-2/data003');
    t=rawFile.getData(start, points);
   
    samplingrate=20000;
   
    for i=1:length(ch)
        subplot(length(ch),1,i)
        plot(t(:,ch(i)+1),'k')
        set(gca,'yLim',[-500 500]);
    end
    subplot(length(ch),1,1)
    hold on
    tt=datarun{1}.spikes{index}*samplingrate;
    ttt=find(tt>start & tt<start+points);
    plot(tt(ttt)-start,0,'r.')
   
end



if 0
    clf
   
    cell_id_ac=5283;ch=[353 348 372 326 240 216];
   
    if 1
        scale=4;thr=0.015;
        cell_id_ac=[6304];
        ch=[421 410 427 461 186 169 163];
        tscale=[-200 200; -100 100; -100 100; -100 100; -100 100; -100 100; -100 100; -100 100; -100 100; -100 100; -100 100; ];
    end
   
    ind_ac=get_cell_indices(datarun{1},cell_id_ac);
   
    if 1
        for i=1:length(ch)
            subplot(length(ch)+2,1,i)
            t=datarun{1}.ei.eis{ind_ac}(ch(i)+0,:);
            plot(t(31:150),'k')
            set(gca,'yLim',tscale(i,:));
        end
        subplot(length(ch)+2,1,length(ch)+[1 2]) 
    end
   
    plot_ei(datarun{1},cell_id_ac,'flipud',1,'fliplr',1,'scale',scale,'cutoff',thr,'sub_scale',thr,'neg_color',[0 0 0],'pos_color',[0 0 0]);              
    hold on
    plot(datarun{1}.ei.position(ch,1),datarun{1}.ei.position(ch,2),'r.')
   
    el=[455 386 264 195 125 4 455]
    plot(datarun{1}.ei.position(el,1)*1.05,datarun{1}.ei.position(el,2)*1.05,'k')
   
    set(gca,'XTick',[],'YTick',[]);
end
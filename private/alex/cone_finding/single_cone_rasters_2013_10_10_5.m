
clear

%%%% 2013-10-10-5 INFO %%%%%%

offM=[721 726 2809 3005 3183 4222 4758 5090] % from data002 for data003 for 2013-10-10-5, 117 cones, stimulated only 92?

piece = '2013-10-10-5';
run = 'data003';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);

datarun=load_data(datarun,opt);
    
stimulus=read_stim_lisp_output_ath('2013-10-10-5','s03');

parsed=parse_stim_rgbs_ath(stimulus);

load(['/Volumes/Data/Auxiliary/2013-10-10-5/stimuli/map-0001.txt'])
figure
imagesc(map_0001)
length(unique(map_0001))


a=parsed.rgbs;

k=cell(1,stimulus.numcones); % stimuli list for each cone
m=zeros(size(stimulus.pulses,2)/stimulus.numcones,stimulus.numcones); % when stimuli were applied
for i=1:size(stimulus.pulses,2)
    p=a{i};
    t=find(p(:,1));
    k{t}=[k{t} p(t,1)];
    p=find(m(:,t)==0,1);
    m(p,t)=i;
end





piece = '2013-10-10-5';
run = 'data002';

% define data path
datarunA = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));  
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = set_polarities(datarunA);
datarunA = load_cones(datarunA,2); % load_cones(datarun, 'Analysis');
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  

[cell_list, failed_cells]=map_ei(datarun,datarunA);

load(['/Volumes/Analysis/deprecated/maps/2013-10-10-5/allconesOffmidgets/cone_indices'])

ind=cell2mat(indices);
kmaps=map_0001;
for i=1:size(ind,2)
    bb=round([datarunA.cones.centers(ind(i),1)*2,datarunA.cones.centers(ind(i),2)*2]);
    indInMap(i)=kmaps(bb(2),bb(1));
end
testedInd=find(indInMap);
emptyInd=find(indInMap==0);

figure
plot(datarunA.cones.centers(:,1)*2,datarunA.cones.centers(:,2)*2,'*')
hold on
plot(datarunA.cones.centers(ind(testedInd),1)*2,datarunA.cones.centers(ind(testedInd),2)*2,'*r')
plot(datarunA.cones.centers(ind(emptyInd),1)*2,datarunA.cones.centers(ind(emptyInd),2)*2,'*g')
axis([0 600 0 600])

figure
% cones which were really tested (colorful regions)
imagesc(kmaps) 
hold on
% all cones found online, red stars
plot(datarunA.cones.centers(:,1)*2,datarunA.cones.centers(:,2)*2,'*r') % all cones map
% cones which should have been tested, black stars
plot(datarunA.cones.centers(ind(testedInd),1)*2,datarunA.cones.centers(ind(testedInd),2)*2,'*k')
% cones which were not tested because of trigger settings? white stars
plot(datarunA.cones.centers(ind(emptyInd),1)*2,datarunA.cones.centers(ind(emptyInd),2)*2,'*','color',[1 1 1])







offM=[721 726 2809 3005 3183 4222 4758 5090];
ccnt=1; clear realresp cone_sta
figure
for mycell=offM % loops through all cells of Interest
    cellOfInterest=mycell; % cell ID in datarunA (cone finding)
    
    mycellInA=find(datarunA.cell_ids==cellOfInterest); % cell index in datarunA (cone finding)
    
    
    a=cellfun(@(x) x==cellOfInterest, cell_list,'UniformOutput',false); % find cell IDs in datarun corresponding to cell ID in datarunA
    a( cellfun(@isempty, a) ) = {false};
    a = cell2mat(a);
    cellIndexInDatarun=find(a); % cell index in datarun (cone testing)
    if ~isempty(cellIndexInDatarun)
        
         
        tr=datarun.triggers(1:2:end);
        mmin=0;mmax=0.75;
        
        kmaps=map_0001;
        tmaps=map_0001;
        for i=cellIndexInDatarun % cell Index in datarun (not ID!)
            spikes=datarun.spikes{i};
            
%             figure
            for j=1:stimulus.numcones
                trig_local=m(k{j}==-0.48,j);
                conv=zeros(size(trig_local,1),750);
                cnt=1;
                for f=trig_local'
                    h=spikes-tr(f);
                    hh= h>=mmin & h<mmax;
                    hh=round(h(hh)*1000);
                    tmp=convolved(hh,40,750);
                    tmp=tmp(121:end-120);
                    conv(cnt,:)=tmp;
                    cnt=cnt+1;
                end
%                 subplot(10,10,j)
%                 plot(nanmean(conv(1:10,:)))
%                 hold on
%                 plot(mean(conv(11:20,:)),'r')
%                 plot(mean(conv(21:30,:)),'k')
%                 plot(mean(conv(31:40,:)),'g')
%                 axis tight
                
                meanResp=nanmean(conv);
                meanRespVal(j,ccnt)=sum(meanResp(1:400))/400;
                peakResp(j,ccnt)=max(meanResp(1:400));
                meanAfterRespVal(j,ccnt)=sum(meanResp(551:700))/150;
                realresp(j,ccnt)=nanmean(nanmean(corr(conv')));                
                kmaps(map_0001==j)=meanRespVal(j,ccnt);
                cone_sta(j,ccnt)=datarunA.cones.weights(ind(j),mycellInA);               
                tmaps(map_0001==j)=datarunA.cones.weights(ind(j),mycellInA);
            end
            
        end
        
        subplot(3,3,ccnt)
%         colormap(gray)
%         imagesc(kmaps)

tmp=kmaps/max(kmaps(:));
tmaps=tmaps/max(abs(tmaps(:)));
a=find(tmaps<0);
tmaps2=zeros(600,600);
tmaps2(a)=-tmaps(a);
tmaps(a)=0;
tmp(:,:,2)=tmaps2;
tmp(:,:,3)=tmaps;
image(tmp)

        hold on
        plot_rf_summaries(datarunA, cellOfInterest, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
        title(['cell ', int2str(cellOfInterest)])
        axis([230 380 300 450])
        xlabel('red: SCS, blue/green: STA')
    
    else
        subplot(3,3,ccnt)
        imagesc(map_0001)
        title(['cell ', int2str(cellOfInterest), ' match not found in data003'])
        plot_rf_summaries(datarunA, cellOfInterest, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
        axis([230 380 300 450])
    end
    ccnt=ccnt+1;
end

respStrength=meanRespVal./meanAfterRespVal;

figure
for i=1:8
    subplot(3,3,i)
    plot(sort(realresp(:,i)),'*')
    axis([0 100 0 1])
end

figure
for i=1:8
    subplot(3,3,i)
    plot(realresp(:,i),meanRespVal(:,i),'*')
%     axis([0 100 0 1])
end

figure
for i=1:8
    subplot(3,3,i)
    plot(realresp(:,i),peakResp(:,i),'*')
%     axis([0 100 0 1])
end


figure
for i=1:8
    subplot(3,3,i)
    plot(cone_sta(:,i),peakResp(:,i),'*')
    title(int2str(offM(i)))
    xlabel('SCS: peak response')
    ylabel('STA: cone weight')
%     axis([0 100 0 1])
end


figure
for i=1:8
    subplot(3,3,i)
    hold on
    plot(cone_sta(:,i),meanRespVal(:,i),'*r')
    title(int2str(offM(i)))
    xlabel('SCS: mean response')
    ylabel('STA: cone weight')
%     axis([0 100 0 1])
end



%% recalculate cone weights in sta by stimulated cone areas (only participating cones)
offM=[721 2809 3005 3183 4758 5090];

datarunA = load_sta(datarunA, struct('load_sta',offM));

new_weight=zeros(stimulus.numcones,length(offM));

for cone=1:stimulus.numcones
    [x,y]=find(map_0001==cone);
    x=ceil(x/2);
    y=ceil(y/2);

    for i=1:length(x)
        for rgc=1:6
            a=squeeze(datarunA.stas.stas{datarunA.cell_ids==offM(rgc),1});
            new_weight(cone,rgc)=new_weight(cone,rgc)+sum(a(x(i),y(i),6));
        end
    end
    
end



%% shared cones
offM=[721 2809 3005 3183 4758 5090];
col=[1 3 3 2 1 2];
rgc=1; clear realresp cone_sta
comb_sta=zeros(600,600,3);
comb_scs=comb_sta;
comb_cons=comb_sta;
comb_new_weight=comb_sta;
allresponses1=zeros(750,length(offM),stimulus.numcones);
allresponses2=zeros(750,length(offM),stimulus.numcones);
spikes_rasters=cell(length(offM),stimulus.numcones);

for mycell=offM % loops through all cells of Interest
    cellOfInterest=mycell; % cell ID in datarunA (cone finding)
    
    mycellInA=find(datarunA.cell_ids==cellOfInterest); % cell index in datarunA (cone finding)
    
    
    a=cellfun(@(x) x==cellOfInterest, cell_list,'UniformOutput',false); % find cell IDs in datarun corresponding to cell ID in datarunA
    a( cellfun(@isempty, a) ) = {false};
    a = cell2mat(a);
    cellIndexInDatarun=find(a); % cell index in datarun (cone testing)
    if ~isempty(cellIndexInDatarun)
        
         
        tr=datarun.triggers(1:2:end);
        mmin=0;mmax=0.75;
        
        kmaps=map_0001;
        tmaps=map_0001;
        cmaps=map_0001;
        nmaps=map_0001;
        for i=cellIndexInDatarun % cell Index in datarun (not ID!)
            spikes=datarun.spikes{i};
            
            for cone=1:stimulus.numcones
                trig_local=m(k{cone}==-0.48,cone);
                conv=zeros(size(trig_local,1),750);
                cnt=1;
                t=[];
                for f=trig_local'
                    h=spikes-tr(f);
                    hh= h>=mmin & h<mmax;
                    hh=round(h(hh)*1000);
                    tmp=convolved(hh,40,750);
                    tmp=tmp(121:end-120);
                    conv(cnt,:)=tmp;
                    cnt=cnt+1;
                    
                    t=[t hh'+750*(cnt-1)];
                end
                spikes_rasters{rgc,cone}=t;
                
                cc=triu(corr(conv(:,1:400)'),1);     

                consist(cone,rgc)=sum(cc(~isnan(cc)))/780;
                
                allresponses1(:,rgc,cone)=nanmean(conv(1:2:end,:));
                allresponses2(:,rgc,cone)=nanmean(conv(2:2:end,:));
                
                meanResp=mean(conv);                
                
                meanRespVal(cone,rgc)=mean(meanResp(1:400));
                peakResp(cone,rgc)=max(meanResp(1:400));
                
                cone_sta(cone,rgc)=datarunA.cones.weights(ind(cone),mycellInA);                
            
            end
            
        end
        
        tt=sort(meanRespVal(:,rgc));
        meanRespVal(:,rgc)=meanRespVal(:,rgc)-mean(tt(1:floor(length(tt)/2)));
        
        for cone=1:stimulus.numcones
            kmaps(map_0001==cone)=meanRespVal(cone,rgc);  % peak response
            tmaps(map_0001==cone)=cone_sta(cone,rgc); % sta weights
            cmaps(map_0001==cone)=consist(cone,rgc); % consistency response
            nmaps(map_0001==cone)=-new_weight(cone,rgc); % multiply by -1 since using OFF midgets
        end
        
        comb_sta(:,:,col(rgc))=comb_sta(:,:,col(rgc))+tmaps/max(tmaps(:)); 
        
        comb_scs(:,:,col(rgc))=comb_scs(:,:,col(rgc))+kmaps/max(kmaps(:));
        
        comb_cons(:,:,col(rgc))=comb_cons(:,:,col(rgc))+cmaps/max(cmaps(:)); 
        
        comb_new_weight(:,:,col(rgc))=comb_new_weight(:,:,col(rgc))+nmaps/max(nmaps(:)); 
    end
    rgc=rgc+1;
end

figure
subplot(1,2,1);
comb_sta(comb_sta>1)=1;
comb_sta(comb_sta<0)=0;
image(comb_sta)
hold on
plot_rf_summaries(datarunA, offM, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
axis([230 380 300 450])
title('STA estimation')

subplot(1,2,2);
a=comb_scs;
% a=comb_new_weight;
a(a>1)=1;
a(a<0)=0;
image(a)
hold on
plot_rf_summaries(datarunA, offM, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
axis([230 380 300 450])
title('Single cone stim estimation, mean response')


getResp=uicontrol('style','pushbutton','string','6 rgc rasters');
set(getResp,'Units', 'normalized','Position',[0.7 0.9 0.07 0.02],'fontsize',12, 'Callback','getResponse(1,spikes_rasters, map_0001,offM,consist,cone_sta,meanRespVal)');


getCones=uicontrol('style','pushbutton','string','25 str cones rasters');
set(getCones,'Units', 'normalized','Position',[0.4 0.9 0.07 0.02],'fontsize',12, 'Callback','getResponse(2,spikes_rasters, map_0001,offM,consist,cone_sta,meanRespVal)');





figure
titles={'original STA', 'cone area summed STA','mean SC response','consistent SC responses'};

allweights={comb_sta,comb_new_weight,comb_scs,comb_cons};
for i=1:4
    subplot(2,2,i);
    a=allweights{i};
    a(a>1)=1;
    a(a<0)=0;
    image(a)
    hold on
    plot_rf_summaries(datarunA, offM, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
    axis([230 380 300 450])
    title(titles{i})
end
getResp=uicontrol('style','pushbutton','string','6 rgc rasters');
set(getResp,'Units', 'normalized','Position',[0.47 0.8 0.07 0.02],'fontsize',12, 'Callback','getResponse(1,spikes_rasters, map_0001,offM,consist,cone_sta,meanRespVal)');
getCones=uicontrol('style','pushbutton','string','25 str cones rasters');
set(getCones,'Units', 'normalized','Position',[0.47 0.9 0.07 0.02],'fontsize',12, 'Callback','getResponse(2,spikes_rasters, map_0001,offM,consist,cone_sta,meanRespVal)');



figure

allweights={comb_sta,comb_new_weight,comb_scs,comb_cons};

a=allweights{1}-allweights{2}+0.5;
image(a)
hold on
plot_rf_summaries(datarunA, offM, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
axis([230 380 300 450])
title('original STA - cone area STA')

%%
figure

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};
for i=1:4
    [~,ia,~]=intersect(datarunA.cell_ids,datarunA.cell_types{i}.cell_ids);
    
    a=sort(datarunA.cones.weights(:,ia),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    subplot(2,2,i)
    
    plot(a)
    hold on 
    plot(mean(a,2),'k','linewidth',3)
    title(titles{i})
    axis([0 100 0 1])
end


figure
col='brgk'
for i=1:4
    [~,ia,~]=intersect(datarunA.cell_ids,datarunA.cell_types{i}.cell_ids);
    
    a=sort(datarunA.cones.weights(:,ia),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    b=mean(a,2);
    ind=find(b<0.1,1);
    xdata=1:100/(ind+1):100;
    p(i)=numel(xdata);

    plot(xdata,b(1:numel(xdata)),col(i),'linewidth',3)
    hold on 
end
legend({['ON parasol, ', int2str(p(1))],['OFF parasol ', int2str(p(2))]...
    ['ON midget, ', int2str(p(3))],['OFF midget, ', int2str(p(4))]})
title(int2str(p))
line([0 100],[0.1 0.1],'linestyle','--','color','k')





figure
col='brgk'
for i=1:4
    [~,ia,~]=intersect(datarunA.cell_ids,datarunA.cell_types{i}.cell_ids);
    
    a=sort(datarunA.cones.weights(:,ia),'descend');
%     a=a./repmat(max(a),size(a,1),1);
    
    b=mean(a,2);
    ind=find(b<0.1,1);

    plot(b,[col(i),'--x'],'linewidth',3)
    hold on 
end
legend({['ON parasol, ', int2str(p(1))],['OFF parasol ', int2str(p(2))]...
    ['ON midget, ', int2str(p(3))],['OFF midget, ', int2str(p(4))]})
title(int2str(p))
line([0 100],[0.1 0.1],'linestyle','--','color','k')



figure
for i=1:4
    [~,ia,~]=intersect(datarunA.cell_ids,datarunA.cell_types{i}.cell_ids);
    
    a=sort(datarunA.cones.weights(:,ia),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    b=mean(a,2);
    ind=find(b<0.1,1);
    xdata=1:100/(ind+1):10000;

    plot(xdata(1:numel(b)),b,[col(i), '--*'],'linewidth',2)
    hold on 
end



%%

figure


thresh=0.2;
find(realresp>thresh)


datarunA = load_sta(datarunA,'load_sta',721);
a=datarunA.stas.stas{19};
a=squeeze(a);
figure
colormap(gray)
imagesc(a(:,:,6))

plot_sta(datarunA,721)

kmaps=map_0000;
for i=indInDatarun%408
    spikes=datarun.spikes{i};
    conv=zeros(29,750);
    figure
    for j=1:29
        trig_local=m(k{j}==-0.48,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        subplot(5,6,j)
        plot(conv(j,:)/50)
        kkk=conv(j,:)/50;
        hold on
        
        conv(j,:)=conv(j,:)-conv(j,:);
        trig_local=m(k{j}==-0.288,j);
        for f=trig_local'
            h=spikes-tr(f);
            hh= h>=mmin & h<mmax;
            hh=round(h(hh)*1000);
            tmp=convolved(hh,40,750);
            tmp=tmp(121:end-120);
            conv(j,:)=conv(j,:)+tmp;
        end
        plot(conv(j,:)/50,'r')
        
        kmaps(map_0000==j)=max(kkk(1:300))-mean(kkk(450:650));
    end
  
end



figure
subplot(1,2,1)
imagesc(map_0001)
colormap(gray)
hold on
plot_rf_summaries(datarunA, [721 726 2809 3005 3183 4222 4758 5090], 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')

subplot(1,2,2)
colormap(gray)
imagesc(kmaps)
hold on
plot_rf_summaries(datarunA, cellOfInterest, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


%%%%%%%% datarunA cone plotting %%%%%%%%
% load data
datarunA = load_params(datarunA,struct('verbose',1));  
datarunA = load_sta(datarunA,'load_sta',[]);
datarunA = set_polarities(datarunA);
% if error loading cones - probably multiple cone maps saved, put the number as second argument
datarunA = load_cones(datarunA); % load_cones(datarun, 'Analysis');
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  
figure

figure
imagesc(map_0000)
hold on
plot_rf_summaries(datarunA, {4}, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarunA, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')

offM_d09 = [151 1396 1741 1816 1863 1966 2026 2116 6796 6991];
offP_d09 = [376 1471 4636 4771 5551];

for i=1:size(datarunA.cell_types,2)
   a=datarunA.cell_types{1, i}.cell_ids;
   intersect(a,offP_d09)
end



figure
plot_rf_summaries(datarunA, {2}, 'clear', false, 'label', true, 'label_color', 'k', 'plot_fits', true, 'fit_color', 'k')
hold on
plot_rf_summaries(datarun, {2}, 'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')


figure
imagesc(map_0000)
hold on
plot_rf_summaries(datarunA, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'r', 'plot_fits', true, 'fit_color', 'r')
plot_rf_summaries(datarun, {2}, 'scale',2,'clear', false, 'label', true, 'label_color', 'y', 'plot_fits', true, 'fit_color', 'y')





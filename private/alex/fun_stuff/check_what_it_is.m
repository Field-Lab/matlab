%% dataruns loading

piece = '2013-10-10-5';
run = 'data002';

% define data path
datarunA = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);

opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunA=load_data(datarunA,opt);
datarunA = load_params(datarunA,struct('verbose',1));  
tic
datarunA = load_sta(datarunA,'load_sta','all');
toc
datarunA = set_polarities(datarunA);
datarunA = load_cones(datarunA,2);
datarunA = make_mosaic_struct(datarunA);
datarunA = get_sta_fits_from_vision(datarunA);  
datarunA = make_voronoi_masks(datarunA);




piece = '2012-09-13-2';
run = 'data009';

% define data path
datarunB = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunB=load_data(datarunB,opt);
datarunB = load_params(datarunB,struct('verbose',1));  
datarunB = load_sta(datarunB,'load_sta',[]);
datarunB = set_polarities(datarunB);
datarunB = load_cones(datarunB,'bayes-msf_5.0');
datarunB = make_mosaic_struct(datarunB);
datarunB = get_sta_fits_from_vision(datarunB);  
datarunB = make_voronoi_masks(datarunB);



piece = '2008-08-27-5';
run = 'data003';

% define data path
datarunC = load_data(['/Volumes/Analysis/' piece '/' run '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunC=load_data(datarunC,opt);
datarunC = load_params(datarunC,struct('verbose',1));  
datarunC = load_sta(datarunC,'load_sta',[]);
datarunC = set_polarities(datarunC);
datarunC = load_cones(datarunC,'standard');
datarunC = make_mosaic_struct(datarunC);
datarunC = get_sta_fits_from_vision(datarunC);  
datarunC = make_voronoi_masks(datarunC);

piece = '2010-03-05-2';
run = 'data001';

% define data path
datarunD = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunD=load_data(datarunD,opt);
datarunD = load_params(datarunD,struct('verbose',1));  
datarunD = load_sta(datarunD,'load_sta',[]);
datarunD = set_polarities(datarunD);
datarunD = load_cones(datarunD);
datarunD = make_mosaic_struct(datarunD);
datarunD = get_sta_fits_from_vision(datarunD);  
datarunD = make_voronoi_masks(datarunD);

piece = '2013-08-19-4';
run = 'data001';

% define data path
datarunE = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunE=load_data(datarunE,opt);
datarunE = load_params(datarunE,struct('verbose',1));  
datarunE = load_sta(datarunE,'load_sta',[]);
datarunE = set_polarities(datarunE);
datarunE = load_cones(datarunE);
datarunE = make_mosaic_struct(datarunE);
datarunE = get_sta_fits_from_vision(datarunE);  
datarunE = make_voronoi_masks(datarunE);


piece = '2008-12-12-1';
run = 'data005';

% define data path
datarunF = load_data(['/Volumes/Analysis/' piece '/' run '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunF=load_data(datarunF,opt);
datarunF = load_params(datarunF,struct('verbose',1));  
datarunF = load_sta(datarunF,'load_sta',[]);
datarunF = set_polarities(datarunF);
datarunF = load_cones(datarunF, 'fit');
datarunF = make_mosaic_struct(datarunF);
datarunF = get_sta_fits_from_vision(datarunF);  
datarunF = make_voronoi_masks(datarunF);


piece = '2008-08-26-2';
run = 'data001-s6369-s9551';

% define data path
datarunG = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunG=load_data(datarunG,opt);
datarunG = load_params(datarunG,struct('verbose',1));  
datarunG = load_sta(datarunG,'load_sta',[]);
datarunG = set_polarities(datarunG);
datarunG = load_cones(datarunG,'indiv');
datarunG = make_mosaic_struct(datarunG);
datarunG = get_sta_fits_from_vision(datarunG);  
datarunG = make_voronoi_masks(datarunG);

piece = '2013-08-19-5';
run = 'data001';

% define data path
datarunH = load_data(['/Volumes/Analysis/' piece '/streamed/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarunH=load_data(datarunH,opt);
datarunH = load_params(datarunH,struct('verbose',1));  
datarunH = load_sta(datarunH,'load_sta',[]);
datarunH = set_polarities(datarunH);
datarunH = load_cones(datarunH);
datarunH = make_mosaic_struct(datarunH);
datarunH = get_sta_fits_from_vision(datarunH);  
datarunH = make_voronoi_masks(datarunH);


%% finding thresholds (inflexion points)
threshold=3;

ft = fittype('max(a*x+b,c*x+d)');

for letterCode='ABCDEFGH'
    letterCode
    datarun=eval(['datarun',letterCode]);
    
    optThreshold=zeros(size(datarun.cell_ids));
    for i=1:4
        [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
        
        
        selectRGCCones=datarun.cones.weights(:,selectRGCind);
        selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
        
        ncones=size(selectRGCCones,1);
        ncells=size(selectRGCCones,2);
        
        [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');
        fitted=~isnan(sum(selectRGCCones));
        
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, datarun.cell_types{i}.cell_ids(fitted),...
            'thresh', threshold, 'radius', [0 inf], 'polarity', 1,'contiguity', true);
        
        j=find(~fitted);
        if length(j)>1
            disp('\nWARNING! more than 1 unfitted\n')
            disp(letterCode)
        end
        for tmp=j
            selection=[selection(:,1:tmp-1), false(ncones,1) selection(:,tmp:end)];
        end
        
        totalCones=sum(selection);

        for currentCell=find(fitted)
            
            lb=round(totalCones(currentCell)/3);
            rb=200;
            x=lb:rb;
            if totalCones(currentCell)<4
                continue
            end
            
            a=allRGCsorted(x,currentCell);
            
            b=min(x);
            x=x-b+1;
            
            res=fit(x',a,ft,'startpoint',[-0.005 max(a) -0.001 0.2]);
            
            
            optCone=find(res.a*(x)+res.b-(res.c*(x)+res.d)<0,1)+b-1;
            new_y=max(res.a*(optCone-b+1)+res.b,res.c*(optCone-b+1)+res.d);
            aa=new_y/optCone;
            
            optCone=find(allRGCsorted(1:200,currentCell)'./(1:200)<aa,1);
            if ~isempty(optCone)
                optThreshold(selectRGCind(currentCell))=optCone;
            end   
        end
        
        for currentCell=datarun.cell_types{i}.cell_ids
            mosaic_weights=datarun.cones.weights(:,datarun.cell_ids==currentCell);
            a=optThreshold(datarun.cell_ids==currentCell);
            tmp=robust_std(mosaic_weights);
            eval(['datarun',letterCode,'.RobustSTD(',int2str(find(datarun.cell_ids==currentCell)),')=',num2str(tmp),';']);
            tmp=std(mosaic_weights);
            eval(['datarun',letterCode,'.STD(',int2str(find(datarun.cell_ids==currentCell)),')=',num2str(tmp),';']);
                
            if a>0
                tmp=sort(mosaic_weights,'descend');
                tmp=tmp/max(mosaic_weights);
                eval(['datarun',letterCode,'.LowestInclWeight(',int2str(find(datarun.cell_ids==currentCell)),')=',num2str(tmp(a-1)),';']);
                eval(['datarun',letterCode,'.HighestExclWeight(',int2str(find(datarun.cell_ids==currentCell)),')=',num2str(tmp(a)),';']);                
                
            else
                eval(['datarun',letterCode,'.LowestInclWeight(',int2str(find(datarun.cell_ids==currentCell)),')=0;']);
                eval(['datarun',letterCode,'.HighestExclWeight(',int2str(find(datarun.cell_ids==currentCell)),')=0;']);
                
            end
        end
        
    end
    
    eval(['datarun',letterCode,'.inflexionPoints=optThreshold;']);
end


%%
a=datarun.cones.weights;
figure
hist(a(:,200),500)
robust_std(a(:,200))
std(a(:,200))

x=-1:0.001:1;
sigm=robust_std(a(:,1));
[b,c]=hist(a(:,200),100);
y=ampl*exp(-((x-sigm)/wid).^2);



%%
figure
col='brgk'
for i=1:4
    [~,selectRGCind,~]=intersect(datarunB.cell_ids,datarunB.cell_types{i}.cell_ids);
    
    a=sort(datarunB.cones.weights(:,selectRGCind),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    b=nanmean(a,2);
    sortedByCell=find(b<0.1,1);
    xdata=1:100/(sortedByCell+1):100;
    p(i)=numel(xdata);

    plot(xdata,b(1:numel(xdata)),col(i),'linewidth',3)
    hold on 
end
legend({'ON parasol','OFF parasol','ON midget','OFF midget'})
title(int2str(p))




figure
col='brgk'
for i=1:4
    [~,selectRGCind,~]=intersect(datarunB.cell_ids,datarunB.cell_types{i}.cell_ids);
    
    a=sort(datarunB.cones.weights(:,selectRGCind),'descend');
%     a=a./repmat(max(a),size(a,1),1);
    
    b=nanmean(a,2);
    sortedByCell=find(b<0.1,1);

    plot(b,[col(i),'--x'],'linewidth',3)
    hold on 
end
legend({['ON parasol, ', int2str(p(1))],['OFF parasol ', int2str(p(2))]...
    ['ON midget, ', int2str(p(3))],['OFF midget, ', int2str(p(4))]})
title(int2str(p))
line([0 100],[0.1 0.1],'linestyle','--','color','k')




figure
for i=1:4
    [~,selectRGCind,~]=intersect(datarunB.cell_ids,datarunB.cell_types{i}.cell_ids);
    
    a=sort(datarunB.cones.weights(:,selectRGCind),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    b=nanmean(a,2);

    plot(b,[col(i), '--*'],'linewidth',2)
    hold on 
end
legend({'ON parasol','OFF parasol','ON midget','OFF midget'})
line([0 size(b,1)],[0.1 0.1],'linestyle','--','color','k')
line([0 size(b,1)],[0 0],'color','k')


%% weight vs distances

datarun=datarunG;

centerCones=3;
meanConeDistances=zeros(size(datarun.cones.weights,1),4);
stErrorData=meanConeDistances;meanWeights=meanConeDistances;
for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones));    

    
    for currentCell=find(fitted)
        current_order=sortedByCell(:,currentCell);
        
        x=datarun.cones.centers(current_order,1);
        y=datarun.cones.centers(current_order,2);
        w=allRGCsorted(:,currentCell);
        
        x_center=x(1:centerCones);
        y_center=y(1:centerCones);
        w_center=w(1:centerCones);
        x_com=sum(w_center.*x_center)/sum(w_center);
        y_com=sum(w_center.*y_center)/sum(w_center);
        
        cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2);
    end
    
    cone_distances(:,~fitted)=[];
    meanConeDistances(:,i)=mean(cone_distances,2);
    stErrorData(:,i)=std(cone_distances,0,2)/sqrt(sum(fitted));    
    meanWeights(:,i)=meanRGCsorted;
end

figure
plot(meanConeDistances,'-x')
xlabel('cone number (sorted by weight)')
ylabel('distance to strongest cone (euclidian)')
legend('ON parasol','OFF parasol','ON midget','OFF midget')
title(['center cones ',int2str(centerCones)])

figure
plot(meanWeights,meanConeDistances,'x')
xlabel('cone weight (norm to strongest)')
ylabel('distance to strongest cone (euclidian)')
legend('ON parasol','OFF parasol','ON midget','OFF midget')


figure
% errorbar(meanData,stErrorData,'-x')
plot(meanConeDistances,meanWeights,'.')
ylabel('cone weight (norm to strongest)')
xlabel('distance to strongest cone (euclidian)')
legend('ON parasol','OFF parasol','ON midget','OFF midget')


%% cone ranks and shared cones

datarun=datarunE;

threshold=0.1;
titles={'ON parasol','OFF parasol','ON midget','OFF midget'};
figure(5)
for i=1:4
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
   
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
     ncones=size(selectRGCCones,1);
     ncells=size(selectRGCCones,2);
    
    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');    
    fitted=~isnan(sum(selectRGCCones));
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    figure(6)
    subplot(2,2,i)
    plot(allRGCsorted)
    hold on
    plot(meanRGCsorted,'linewidth',3,'color','k')
    title([titles{i},', N=',int2str(sum(fitted))])
    line([0 ncones],[0 0],'color','k','linewidth',2)
    line([0 ncones],[threshold threshold],'color','r','linewidth',3)
    axis tight
    
    
    sharedConeWeights=zeros(size(selectRGCCones));
 
    sharedConesNumber=zeros(size(selectRGCCones,1),1);
    
    for currentCell=find(fitted)
        
        allOtherCells=true(1,ncells);
        allOtherCells(currentCell)=false;
        sortByCurrentCell=selectRGCCones(sortedByCell(:,currentCell),:);
        
%         tmp=repmat(allRGCsorted(:,currentCell),1,ncells);
%         tmp=sortByCurrentCell>tmp;
        
        tmp=sortByCurrentCell>threshold;
        
        a=sortByCurrentCell(:,allOtherCells&fitted).*tmp(:,allOtherCells&fitted);       

        a=sum(a,2)./sum(tmp(:,allOtherCells&fitted),2);
        a(isnan(a))=0;
        
        sharedConeWeights(:,currentCell)=a;

        sharedConesNumber=sharedConesNumber+sum(tmp(:,allOtherCells),2);
%         if sum(tmp(1,allOtherCells))>0
%             currentCell
%         end
    end
    sharedConesNumber=sharedConesNumber/sum(fitted);
%     meanSharedConeWeights=mean(sharedConeWeights(:,fitted),2);
    a=sharedConeWeights>0;
    meanSharedConeWeights=sum(sharedConeWeights,2)./sum(a,2);

    figure(5)
    subplot(2,2,i)
    plot(meanRGCsorted,'-ok')
    hold on
    plot(meanSharedConeWeights,'-xb')
    plot(sharedConesNumber,'*r')
    legend('mean weights, all','shared weights','share of cells')
    line([0 ncones],[threshold threshold],'linestyle','--','color','k')
    line([0 ncones],[0 0],'color','k')
    title([titles{i},', N=',int2str(sum(fitted))])
%     axis([0 100 0 Inf])
end
    

%% finding cone threshold

datarun=datarunD;

threshold=3;

ft = fittype('max(a*x+b,c*x+d)');

for i=1:4
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
   
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
     ncones=size(selectRGCCones,1);
     ncells=size(selectRGCCones,2);
    
    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');    
    fitted=~isnan(sum(selectRGCCones));
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, datarun.cell_types{i}.cell_ids,...
            'thresh', threshold, 'radius', [0 inf], 'polarity', 1,'contiguity', true);
        
    totalCones=sum(selection);
        
    figure
    cc=[];
    for currentCell=find(fitted) 
        
%         lb=round(totalCones(currentCell)/2);
%         rb=round(totalCones(currentCell)/2*3);
lb=round(totalCones(currentCell)/3);
rb=200;
x=lb:rb;
        if totalCones(currentCell)<4
            continue
        end
        
        a=allRGCsorted(x,currentCell);
        
        b=min(x);
        x=x-b+1;
        
        res=fit(x',a,ft,'startpoint',[-0.005 max(a) -0.001 0.2]); 
        
%         y1=allRGCsorted(1:200,currentCell);
%         res_acc=[];
%         for j=1:10:180
%             res=fit([j:j+9]',y1(j:j+9),'poly1');
%             res_acc=[res_acc res.p1];
%         end
        
%         figure
%         plot(x+b-1,a)
%         hold on
%         plot(x,max(res.a*(x-b+1)+res.b,res.c*(x-b+1)+res.d),'r')
        
        c=find(res.a*(x)+res.b-(res.c*(x)+res.d)<0,1)+b-1;
        new_y=max(res.a*(c-b+1)+res.b,res.c*(c-b+1)+res.d);
        aa=new_y/c;
        
        if currentCell<36
        subplot(6,6,currentCell)
        hold on
        plot(0:200,aa*(0:200),'r')        
        c=find(allRGCsorted(1:200,currentCell)'./(1:200)<aa,1);
        
        plot(allRGCsorted(:,currentCell))
        
%         aa=new_y/c;
%         x=0:0.11;
%         new_line=aa*x;
        

        plot(c,allRGCsorted(c,currentCell),'*g')
        title([int2str(c), '  w=',num2str(allRGCsorted(c,currentCell),3),', tot con ',int2str(totalCones(currentCell))])
        
        plot(x,res.a*(x-b+1)+res.b,'m')
        plot(x,res.c*(x-b+1)+res.d,'m')
        line([lb lb],[-0.1 1],'color','k')
        line([rb rb],[-0.1 1],'color','k')
        axis([0 200 -.1 1])
        cc=[cc aa];
        end
        

    end

end
    


%% cone ranks and shared cones - PULLING
letterCode='ABCDFG';
threshold=0.1;
titles={'ON parasol','OFF parasol','ON midget','OFF midget'};
figure
col='rgbm';
for whichRun=1:6
    datarun=eval(['datarun',letterCode(whichRun)]);
    
    for i=1:4
        [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
        
        
        selectRGCCones=datarun.cones.weights(:,selectRGCind);
        selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
        
        ncones=size(selectRGCCones,1);
        ncells=size(selectRGCCones,2);
        
        [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');
        fitted=~isnan(sum(selectRGCCones));
        meanRGCsorted=mean(allRGCsorted(:,fitted),2);
        
%         xdata=1:99/(ncones-1):100;
xdata=1:ncones;
        hold on
        plot(xdata,meanRGCsorted,'linewidth',2,'color',col(i))
        if whichRun==1&&i==4
            legend('ON parasol','OFF parasol','ON midget','OFF midget')
        end

        axis tight

    end
end
line([0 100],[0 0],'color','k','linewidth',1)
line([0 100],[threshold threshold],'color','r','linewidth',1)


%% cone ranks and distances

datarun=datarunG;

threshold=0.1;
titles={'ON parasol','OFF parasol','ON midget','OFF midget'};
figure(5)
for i=1:4
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
   
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
     ncones=size(selectRGCCones,1);
     ncells=size(selectRGCCones,2);
    
    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');    
    fitted=~isnan(sum(selectRGCCones));
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    figure(6)
    subplot(2,2,i)
    plot(allRGCsorted)
    hold on
    plot(meanRGCsorted,'linewidth',3,'color','k')
    title(titles{i})
    line([0 ncones],[0 0],'color','k','linewidth',2)
    line([0 ncones],[threshold threshold],'color','r','linewidth',3)
    axis tight
    
    
    sharedConeWeights=zeros(size(selectRGCCones));
 
    sharedConesNumber=zeros(size(selectRGCCones,1),1);
    
    for currentCell=find(fitted)
        
        allOtherCells=true(1,ncells);
        allOtherCells(currentCell)=false;
        sortByCurrentCell=selectRGCCones(sortedByCell(:,currentCell),:);
        
%         tmp=repmat(allRGCsorted(:,currentCell),1,ncells);
%         tmp=sortByCurrentCell>tmp;
        
        tmp=sortByCurrentCell>threshold;
        
        a=sortByCurrentCell(:,allOtherCells&fitted).*tmp(:,allOtherCells&fitted);       

        a=sum(a,2)./sum(tmp(:,allOtherCells&fitted),2);
        a(isnan(a))=0;
        
        sharedConeWeights(:,currentCell)=a;

        sharedConesNumber=sharedConesNumber+sum(tmp(:,allOtherCells),2);
%         if sum(tmp(1,allOtherCells))>0
%             currentCell
%         end
    end
    sharedConesNumber=sharedConesNumber/sum(fitted);
%     meanSharedConeWeights=mean(sharedConeWeights(:,fitted),2);
    a=sharedConeWeights>0;
    meanSharedConeWeights=sum(sharedConeWeights,2)./sum(a,2);

    figure(5)
    subplot(2,2,i)
    plot(meanRGCsorted,'-ok')
    hold on
    plot(meanSharedConeWeights,'-xb')
    plot(sharedConesNumber,'*r')
    legend('mean weights, all','shared weights','share of cells')
    line([0 ncones],[threshold threshold],'linestyle','--','color','k')
    line([0 ncones],[0 0],'color','k')
    title([titles{i},', N=',int2str(sum(fitted))])
%     axis([0 100 0 Inf])
end
    

%%

figure
col='brgk'
for i=1:4
    [~,selectRGCind,~]=intersect(datarunA.cell_ids,datarunA.cell_types{i}.cell_ids);
    
    a=sort(datarunA.cones.weights(:,selectRGCind),'descend');
    a=a./repmat(max(a),size(a,1),1);
    
    b=mean(a,2);
    sortedByCell=find(b<0.1,1);
    xdata=1:100/(sortedByCell+1):100;
    p(i)=numel(xdata);

    plot(xdata,b(1:numel(xdata)),col(i),'linewidth',3)
    hold on 
end
legend({'ON parasol','OFF parasol','ON midget','OFF midget'})
title(int2str(p))



%%
clear scpieces allpieces mypieces
cnt=1;
for i=1:size(athstixelSizeLess3,1)
        if ~isempty(athstixelSizeLess3{i})            
            a=athstixelSizeLess3{i};
            if ~isempty(find(a==' ',1))
                a=a(1:find(a==' ',1)-1);
            end
            if ~isstrprop(a(7),'digit')
                a=[a(1:5),'0',a(6:end)];
            end
            if ~isstrprop(a(10),'digit')
                a=[a(1:8),'0',a(9:end)];                
            end
            scpieces(cnt,1:12)=a;
            cnt=cnt+1;
        end            
end

tmp=dir('/Volumes/Analysis/dataset_overview/*.pdf');
for i=1:length(tmp)
    allpieces{i}=tmp(i).name(1:12);    
end

cnt=1;piece_number=[];
for i=1:length(scpieces)
    a=strcmp(allpieces,scpieces(i,:));
    if ~isempty(find(a))
        mypieces{cnt}=allpieces{find(a)};
        piece_number=[piece_number find(a)];
        cnt=cnt+1;
    end
end

tmp=dir('/Volumes/Analysis/dataset_overview/*.pdf');
for i=piece_number
    copyfile(['/Volumes/Analysis/dataset_overview/',tmp(i).name],...
        ['/Users/alexth/Desktop/single_cones_candidates/',tmp(i).name])    
end


%% Model

a=rand(10000,3);
[w, ia]=sort(a(:,3),'descend');
x=a(ia,1)*100;
y=a(ia,2)*100;
figure


for i=1:size(a,1)
    centerX = a(i,1);
    centerY = a(i,2);
    radius = a(i,3)/10;
    rectangle('Position',[centerX - radius/2, centerY - radius/2, radius, radius],...
        'Curvature',[1,1],...
        'FaceColor','r');
    axis([0 100 0 100]);
    hold on
end


centerCones=3;
x_center=x(1:centerCones);
y_center=y(1:centerCones);
w_center=w(1:centerCones);
x_com=sum(w_center.*x_center)/sum(w_center);
y_com=sum(w_center.*y_center)/sum(w_center);

cone_distances=sqrt((x-x_com).^2+(y-y_com).^2);
cone_distances=cone_distances/max(cone_distances);

% reverse proportional
figure
w=ones(size(cone_distances))./cone_distances;
w=w/max(w);
plot(w,cone_distances,'*')
xlabel('weight')
ylabel('distance')

% linear
figure
w=(1-cone_distances);
w=w/max(w);
plot(w,cone_distances,'*')
xlabel('weight')
ylabel('distance')

% random
figure
w=rand(size(cone_distances));
w=w/max(w);
plot(w,cone_distances,'*')
xlabel('weight')
ylabel('distance')


% gaussian
figure
mu=0;
sigma=0.3;
w=mu+sigma.*randn(size(cone_distances,1)*3,1);
w=-w(w<=0);
w=w/max(w);
w=sort(w,'descend');
w=w(1:10000);
[a,ai]=sort(cone_distances);
plot(w,a,'*')
xlabel('weight')
ylabel('distance')





datarun=datarunF;

threshold=0.1;
titles={'ON parasol','OFF parasol','ON midget','OFF midget'};
figure(5)
for i=1:4
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
   
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
     ncones=size(selectRGCCones,1);
     ncells=size(selectRGCCones,2);
    
    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');    
    fitted=~isnan(sum(selectRGCCones));
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    figure(6)
    subplot(2,2,i)
    plot(allRGCsorted)
    hold on
    plot(meanRGCsorted,'linewidth',3,'color','k')
    title([titles{i},', N=',int2str(sum(fitted))])
    line([0 ncones],[0 0],'color','k','linewidth',2)
    line([0 ncones],[threshold threshold],'color','r','linewidth',3)
    axis tight
    
    
    sharedConeWeights=zeros(size(selectRGCCones));
 
    sharedConesNumber=zeros(size(selectRGCCones,1),1);
    
    for currentCell=find(fitted)
        
        allOtherCells=true(1,ncells);
        allOtherCells(currentCell)=false;
        sortByCurrentCell=selectRGCCones(sortedByCell(:,currentCell),:);
        
%         tmp=repmat(allRGCsorted(:,currentCell),1,ncells);
%         tmp=sortByCurrentCell>tmp;
        
        tmp=sortByCurrentCell>threshold;
        
        a=sortByCurrentCell(:,allOtherCells&fitted).*tmp(:,allOtherCells&fitted);       

        a=sum(a,2)./sum(tmp(:,allOtherCells&fitted),2);
        a(isnan(a))=0;
        
        sharedConeWeights(:,currentCell)=a;

        sharedConesNumber=sharedConesNumber+sum(tmp(:,allOtherCells),2);
%         if sum(tmp(1,allOtherCells))>0
%             currentCell
%         end
    end
    sharedConesNumber=sharedConesNumber/sum(fitted);
%     meanSharedConeWeights=mean(sharedConeWeights(:,fitted),2);
    a=sharedConeWeights>0;
    meanSharedConeWeights=sum(sharedConeWeights,2)./sum(a,2);

    figure(5)
    subplot(2,2,i)
    plot(meanRGCsorted,'-ok')
    hold on
    plot(meanSharedConeWeights,'-xb')
    plot(sharedConesNumber,'*r')
    legend('mean weights, all','shared weights','share of cells')
    line([0 ncones],[threshold threshold],'linestyle','--','color','k')
    line([0 ncones],[0 0],'color','k')
    title([titles{i},', N=',int2str(sum(fitted))])
%     axis([0 100 0 Inf])
end
    





datarun=datarunF;
figure
nOfCones=5066;
nOfCells=1;
centerCones=1;
meanConeDistances=zeros(size(datarun.cones.weights,1),4);
stErrorData=meanConeDistances;meanWeights=meanConeDistances;
for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones));    

    
    for currentCell=find(fitted)
        current_order=sortedByCell(:,currentCell);
        
        x=datarun.cones.centers(current_order,1);
        y=datarun.cones.centers(current_order,2);
        w=allRGCsorted(:,currentCell);
        
        x_center=x(1:centerCones);
        y_center=y(1:centerCones);
        w_center=w(1:centerCones);
        x_com=sum(w_center.*x_center)/sum(w_center);
        y_com=sum(w_center.*y_center)/sum(w_center);
        
        cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2);
    end
    
    cone_distances(:,~fitted)=[];
    meanConeDistances(:,i)=mean(cone_distances,2);
    stErrorData(:,i)=std(cone_distances,0,2)/sqrt(sum(fitted));    
    meanWeights(:,i)=meanRGCsorted;
    
    subplot(2,2,i)
    a=allRGCsorted(:,fitted);
%     a(a<0.5)=nan;    
    b=cone_distances;
%     b(isnan(a))=nan;
    plot(b(1:nOfCones,1:end),a(1:nOfCones,1:end),'.k')
    hold on
    plot(nanmean(b(1:nOfCones,:),2),nanmean(a(1:nOfCones,:),2),'k.-','markersize',10,'linewidth',3)
    xlabel('distance')
    ylabel('weight')
    axis([0 75 -0.5 1])
end


figure
plot(meanWeights,meanConeDistances,'x')
xlabel('cone weight (norm to strongest)')
ylabel('distance to strongest cone (euclidian)')
legend('ON parasol','OFF parasol','ON midget','OFF midget')




figure;
a=allRGCsorted(:,fitted);
plot(cone_distances(1:30,:),a(1:30,:),'-x')
xlabel('distance')
ylabel('weight')
[a,ai]=sort(meanWeights(:,1),'descend');
th=find(a<0.2,1);
a=a(1:th);
b=meanConeDistances(ai,1);
b=b(1:th);
b=b/max(b);
figure
plot(a,b,'.')
hold on


a=rand(10000,3);
[w, ia]=sort(a(:,3),'descend');
x=a(ia,1)*100;
y=a(ia,2)*100;
centerCones=3;
x_center=x(1:centerCones);
y_center=y(1:centerCones);
w_center=w(1:centerCones);
x_com=sum(w_center.*x_center)/sum(w_center);
y_com=sum(w_center.*y_center)/sum(w_center);

cone_distances=sqrt((x-x_com).^2+(y-y_com).^2);
cone_distances=cone_distances/max(cone_distances);

mu=0;
sigma=0.3;
w=mu+sigma.*randn(size(cone_distances,1)*5,1);
w=-w(w<=-.2);
w=w/max(w);
w=sort(w,'descend');
w=w(1:10000);
[a,ai]=sort(cone_distances);
plot(w,a,'r.')
xlabel('weight')
ylabel('distance')






datarun=datarunG;
figure
centerCones=1;
b=0; % mu
c=0.1:0.5:50; % sigma parameter
radius=0.1:0.5:50; % inner cones inclusion radius
remove_surround=true;

for i=1:4
    means=zeros(100,3);
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones));    
    cellCenters=zeros(2,numel(fitted));
    
    for currentCell=find(fitted)
        current_order=sortedByCell(:,currentCell);
        
        x=datarun.cones.centers(current_order,1);
        y=datarun.cones.centers(current_order,2);
        w=allRGCsorted(:,currentCell);
        
        x_center=x(1:centerCones);
        y_center=y(1:centerCones);
        w_center=w(1:centerCones);
        x_com=sum(w_center.*x_center)/sum(w_center);
        y_com=sum(w_center.*y_center)/sum(w_center);
        cellCenters(1,currentCell)=x_com;
        cellCenters(2,currentCell)=y_com;
        
        cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2);
        
        
        check=zeros(100,sum(fitted));checkL=zeros(1,sum(fitted));checkR=check;
        old_w=w;
        if remove_surround
            old_w(w<0)=0;
        end
        for cnt=1:100
            new_wR=zeros(size(w));
            %gaussian
            new_w=exp(-(cone_distances(:,currentCell)-b).^2/(2*c(cnt)^2));
            check(cnt,currentCell)=corr(old_w,new_w);
            
            %linear
            if cnt==1
                new_wL=ones(size(cone_distances(:,currentCell)))./(cone_distances(:,currentCell)+1);
                checkL(cnt,currentCell)=corr(old_w,new_wL);
            end
            %random within certain radius
            innerCones=find(cone_distances(:,currentCell)<radius(cnt));
            new_wR(innerCones)=rand(size(innerCones));
            checkR(cnt,currentCell)=corr(old_w,new_wR);
            
        end


    end
    means(:,1)=mean(check,2);
    means(:,2)=repmat(mean(checkL),100,1);
    means(:,3)=mean(checkR,2);    
    
    subplot(2,2,i)
    hold on
    plot(means,'linewidth',2)
    legend('Gaussian','Linear','Random')
    title(titles{i})
end






currentCell=10;
current_order=sortedByCell(:,currentCell);
x=datarun.cones.centers(current_order,1);
y=datarun.cones.centers(current_order,2);
w=allRGCsorted(:,currentCell);
figure
for i=1:size(x,1)
    centerX = x(i);
    centerY = y(i);
    radius = abs(w(i))*5;
    if w(i)>0
        col=[w(i) 0 0];
    else
        col=[0 -w(i) -w(i)];
    end
    rectangle('Position',[centerX - radius/2, centerY - radius/2, radius, radius],...
        'Curvature',[1,1],...
        'FaceColor',col);
    hold on
end




plot(cellCenters(1,:),cellCenters(2,:),'*')



%% weight fraction vs max distance


datarun=datarunF;
figure
fractions=-1:0.05:1;
centerCones=1;
meanConeDistances=zeros(size(datarun.cones.weights,1),4);
stErrorData=meanConeDistances;
meanWeights=meanConeDistances;
for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones));    
    maxDist=zeros(length(fractions),size(selectRGCCones,2));
    meanDist=maxDist;
    minDist=maxDist;
    
    for currentCell=find(fitted)
        current_order=sortedByCell(:,currentCell);
        
        x=datarun.cones.centers(current_order,1);
        y=datarun.cones.centers(current_order,2);
        w=allRGCsorted(:,currentCell); % weights for current cell, sorted in descending order
        
        x_center=x(1:centerCones);
        y_center=y(1:centerCones);
        w_center=w(1:centerCones);
        x_com=sum(w_center.*x_center)/sum(w_center);
        y_com=sum(w_center.*y_center)/sum(w_center);
        
        cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights
        
        % weight fractioning
        for fcnt=1:length(fractions)-1
            a=find(w>fractions(fcnt)&w<=fractions(fcnt+1));
            if ~isempty(a)
                maxDist(fcnt,currentCell)=max(cone_distances(a,currentCell));
                meanDist(fcnt,currentCell)=mean(cone_distances(a,currentCell));
                minDist(fcnt,currentCell)=min(cone_distances(a,currentCell));
            end
        end
    end
    
    maxDist(:,~fitted)=[];
    meanDist(:,~fitted)=[];
    minDist(:,~fitted)=[];
    
    subplot(2,2,i)
    a=allRGCsorted(:,fitted);
%     a(a<0.5)=nan;    
    b=cone_distances;
%     b(isnan(a))=nan;
    plot(maxDist,fractions,'xr')
    hold on
    plot(meanDist,fractions,'xk')
    plot(minDist,fractions,'xb')
    
    xlabel('distance')
    ylabel('weight (binned)')
end



%% weight vs uneven binned distance

datarun=datarunH; % H A B C D E F work, G has too few cells
figure
sigmas=10; % how many sigmas of RF fit to take (defines the distance from the cell center to consider)
threshCS=5; % how many sigmas of noise to cut off 'noisy cones', to define 'cell center' location
col='xr';


for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones)); 
    
    rfSigmas=zeros(size(selectRGCCones,2),1);
    
    
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, datarun.cell_types{i}.cell_ids,...
    'thresh', threshCS, 'radius', [0 inf], 'polarity', 1, 'contiguity', false);
    
    for currentCell=find(fitted)
        current_order=sortedByCell(:,currentCell);
        tmp=datarun.stas.fits{selectRGCind(currentCell)}.sd;
        tmp=sqrt(prod(tmp));
        rfSigmas(currentCell)=tmp;
        
        x=datarun.cones.centers(current_order,1);
        y=datarun.cones.centers(current_order,2);
        w=allRGCsorted(:,currentCell); % weights for current cell, sorted in descending order
        
        tmp_selection=selection(current_order,currentCell);

        x_center=x(tmp_selection);
        y_center=y(tmp_selection);
        w_center=w(tmp_selection);

        x_com=sum(w_center.*x_center)/sum(w_center);
        y_com=sum(w_center.*y_center)/sum(w_center);
        
        cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2); % distances for current cell, ordered by weights

    end
    
    rfSigmas(~fitted)=[];
    cone_distances(:,~fitted)=[];
    allRGCsorted(:,~fitted)=[];
    
    a=sort(cone_distances(:));
    
    pointsPerBin=sum(fitted);
    
    meanRFSigma=mean(rfSigmas);
    
    endPoint=find(a>sigmas*meanRFSigma,1);
    meanWeight=zeros(numel(pointsPerBin:pointsPerBin:endPoint),1);
    stWeight=meanWeight;meanDist=meanWeight;stDist=meanWeight;
    maxDist=meanWeight;
    cnt=1;
    
    for bin=pointsPerBin:pointsPerBin:endPoint     
        meanWeight(cnt)=mean(allRGCsorted(cone_distances<=a(bin)));
        stWeight(cnt)=std(allRGCsorted(cone_distances<=a(bin)))/sqrt(pointsPerBin);
        meanDist(cnt)=mean(cone_distances(cone_distances<=a(bin)));
        stDist(cnt)=std(cone_distances(cone_distances<=a(bin)))/sqrt(pointsPerBin);
        maxDist(cnt)=a(bin);
        cone_distances(cone_distances<=a(bin))=1000;
        cnt=cnt+1;
    end
    meanWeight=meanWeight/max(meanWeight(:));

    subplot(2,2,i)
    hold on
    errorbar(meanDist,meanWeight,stWeight, col)
    plot(meanDist,meanWeight,'k')
    line([0 a(bin)],[0 0],'color','k')
    axis([0 a(bin) -0.2 1])
    set(gca,'ytick',[-0.2 0 0.5 1],'fontsize',16)
%     xlabel('distance')
%     ylabel('weight')
    
    [fit_res,gof]=fit(meanDist,meanWeight,'gauss2',...
        'Lower',[0.5 -0.00001 0 -1 -0.00001 0],'Upper',[5 0.00001 Inf 0 0.000001 Inf],...
        'StartPoint',[1 0 meanRFSigma -0.1 0 4*meanRFSigma]);
    
    x=meanDist;
    y=fit_res.a1*exp(-((x-fit_res.b1)/fit_res.c1).^2)+...
        fit_res.a2*exp(-((x-fit_res.b2)/fit_res.c2).^2);
    plot(x,y,'b','linewidth',3)
%     text(a(bin)/2,0.5,num2str([fit_res.a1; fit_res.b1;fit_res.c1;fit_res.a2;fit_res.b2;fit_res.c2]))
    
%     title([titles{i},', bin size (N of cells) ', int2str(pointsPerBin), ', ',int2str(cnt),' bins, sigma ',num2str(meanRFSigma,2),', Rsq=',num2str(gof.rsquare,2)])
    drawnow
    
end


% relation of noise cones and cone weight

threshCS=[2:10];
col='rgbk'

figure
for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));


    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend');
    
    cone_distances=zeros(size(selectRGCCones));
    
    minWeight=zeros(length(threshCS),size(fitted,2));
    
    for cnt=1:length(threshCS)
        [mosaic_weights, selection, extras] = select_cone_weights(datarun, datarun.cell_types{i}.cell_ids,...
            'thresh', threshCS(cnt), 'radius', [0 inf], 'polarity', 1,'contiguity', false);
        
        
        %     figure
        
        for currentCell=find(fitted)
            current_order=sortedByCell(:,currentCell);
            w=allRGCsorted(:,currentCell); % weights for current cell, sorted in descending order
            tmp_selection=selection(current_order,currentCell);
            
            x_center=x(tmp_selection);
            y_center=y(tmp_selection);
            w_center=w(tmp_selection);
            
            x_com=sum(w_center.*x_center)/sum(w_center);
            y_com=sum(w_center.*y_center)/sum(w_center);
            
            cone_distances(:,currentCell)=sqrt((x-x_com).^2+(y-y_com).^2);
            
            if sum(tmp_selection)
                minWeight(cnt,currentCell)=min(w(tmp_selection));
            else
                minWeight(cnt,currentCell)=nan;
            end
            
            
            %         plot(cone_distances(:,currentCell),w)
            %         hold on
            %         line([0 max(cone_distances(:,currentCell))],[thresh thresh],'color','k')
            %         plot(cone_distances(tmp_selection,currentCell), w(tmp_selection),'*r')
            %         a= w>thresh & ~tmp_selection;
            %         plot(cone_distances(a,currentCell), w(a),'*g')
        end
    end
    plot(threshCS,nanmean(minWeight,2),col(i))
    hold on
%     plot(minWeight,col(i));
%     hold on

end
legend(titles)


%% plotting utilities

% plot rf fits

opts=struct;
opts.clear=false;
opts.label=true;
opts.label_color='k';
opts.plot_fits=true;
opts.fit_color='k';

titles={'ON parasol','OFF parasol','ON midget','OFF midget'};

figure
for i=1:4
    subplot(2,2,i)
    plot_rf_summaries(datarun, titles{i}, opts)
    title(titles{i})
end




% plot rfs in cones
datarun=datarunH;

threshCS=5;
thresh=0.3;
% ax=[0 600 0 300];
ax=[0 200 0 200];

figure
for i=1:4
    % find indices of cells of a given type
    [~,selectRGCind,~]=intersect(datarun.cell_ids,datarun.cell_types{i}.cell_ids);
    
    selectRGCCones=datarun.cones.weights(:,selectRGCind);
    selectRGCCones=selectRGCCones./repmat(max(selectRGCCones),size(selectRGCCones,1),1);
    
    
    % all cells which have non zero cone weigths
    fitted=~isnan(sum(selectRGCCones));

    [allRGCsorted, sortedByCell]=sort(selectRGCCones,'descend'); 
    meanRGCsorted=mean(allRGCsorted(:,fitted),2);
    
    cone_distances=zeros(size(selectRGCCones)); 
    
    rfSigmas=zeros(size(selectRGCCones,2),1);
    
    w=selectRGCCones; % normalized cone weights
    x=datarun.cones.centers(:,1);
    y=datarun.cones.centers(:,2);
    
    [mosaic_weights, selection, extras] = select_cone_weights(datarun, datarun.cell_types{i}.cell_ids,...
        'thresh', threshCS, 'radius', [0 inf], 'polarity', 1,'contiguity', false);
    
    subplot(2,2,i)
    plot_rf_summaries(datarun, titles{i}, opts)
    title([titles{i},', noise sigma ',int2str(threshCS),' red - selected, green >',num2str(thresh)])
    axis(ax);
    hold on
    
    for currentCell=find(fitted);
        
        tmp=find(abs(w(:,currentCell)>thresh))';
        
        for j=tmp
            centerX = x(j);
            centerY = y(j);
            radius = abs(w(j,currentCell))*3;
            if w(j,currentCell)>thresh
                if selection(j,currentCell)
                    col='r';
                else
                    col='g';
                end
            else
                col='b';
            end
            rectangle('Position',[centerX - radius/2, centerY - radius/2, radius, radius],...
                'Curvature',[1,1],'FaceColor',col);
        end
        
        drawnow
    end

end





piece = '2012-09-13-2';
run = 'data009';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);


datarunB60_snl =datarun;
datarunB60_snl = load_cones(datarunB60_snl,'bayes-msf_60.00-snlTrue_');
datarunB60_snl = make_mosaic_struct(datarunB60_snl);
datarunB60_snl = get_sta_fits_from_vision(datarunB60_snl);  
datarunB60_snl = make_voronoi_masks(datarunB60_snl);


datarunB20_prior = datarun;
datarunB20_prior = load_cones(datarunB20_prior,'bayes-msf_20.00-prior3.5_4_');
datarunB20_prior = make_mosaic_struct(datarunB20_prior);
datarunB20_prior = get_sta_fits_from_vision(datarunB20_prior);  
datarunB20_prior = make_voronoi_masks(datarunB20_prior);



datarunB30_prior = datarun;
datarunB30_prior = load_cones(datarunB30_prior,'bayes-msf_30.00-prior3.5_4_');
datarunB30_prior = make_mosaic_struct(datarunB30_prior);
datarunB30_prior = get_sta_fits_from_vision(datarunB30_prior);  
datarunB30_prior = make_voronoi_masks(datarunB30_prior);



datarunB5_Scones = datarun;
datarunB5_Scones = load_cones(datarunB5_Scones,'bayes-msf_5.00-prior7_Scones_');
datarunB5_Scones = make_mosaic_struct(datarunB5_Scones);
datarunB5_Scones = get_sta_fits_from_vision(datarunB5_Scones);  
datarunB5_Scones = make_voronoi_masks(datarunB5_Scones);


datarun = load_sta(datarun,'load_sta','all');


a=datarunB60_snl.cones.centers;
b=datarunB20_prior.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('bayes 60','bayes 40 prior')





D=pdist2(a,b);
figure
m=min(D);
mrev=min(D')
hist(m,50)
j

y=find(m<=2); % both for bayesian matrix
y2=find(m>2); % unique bayesian cones
x2=find(mrev>2); % unique local max cones

% y1=find(m>1 & m<=3.5);

figure
hold on
plot(a(x2,1),a(x2,2),'.r','markerSize',30) % unique local max cones
plot(b(y2,1),b(y2,2),'.b','markerSize',30) % unique bayesian cones
plot(b(y,1),b(y,2),'.c','markerSize',30) % both cones (from bayesian matrix)
legend('bayes 60','bayes 10 prior','both')
plot(b(:,1),b(:,2),'+')
plot(a(:,1),a(:,2),'xr')


aw=datarunB60_snl.cones.weights;
bw=datarunB20_prior.cones.weights;



[x,y]=find(D<2);


figure
plot(b(y,1),b(y,2),'+')
hold on
plot(a(x,1),a(x,2),'xr')
title('only found by both methods')



% SBC
sbc=datarun.cell_types{5}.cell_ids;

cca=datarunB60_snl.cones.centers;
ccb=datarunB20_prior.cones.centers;

for tmp=1:length(sbc)
    
    unit=find(datarun.cell_ids==sbc(tmp),1);
    
    if nnz(aw(:,unit))
            
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[0 0 curcones(i)];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30);
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*40/(mean(curconesAbs)*30);
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[curconesA(i) 0 0];
                sizA=(curconesAbsA(i))*40/(mean(curconesAbsA)*30);
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*40/(mean(curconesAbsA)*30);
            end
            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
        end
        
        title(['SBC, i=',int2str(unit),',  blue bayes 10 prior, red bayes 10'])
    end
end



aUnique=setxor(x,1:length(a));
bUnique=setxor(y,1:length(b));


% unique-found, bayes 15
normalizedWeightsA=aw./repmat(max(aw),size(aw,1),1);
uniqueWeightsA=zeros(size(aw));
uniqueWeightsA(aUnique,:)=1;
myWeightsA=uniqueWeightsA&normalizedWeightsA>0.8;
myCellsA=find(sum(myWeightsA));

relSizeSur=50;
relSizeCentr=40;

for unit=myCellsA
    
    if nnz(aw(:,unit))&&nnz(bw(:,unit))
        
        if ~isempty(find(datarun.cell_types{1}.cell_ids==datarun.cell_ids(unit),1))
            tit='ON parasol';
        elseif ~isempty(find(datarun.cell_types{2}.cell_ids==datarun.cell_ids(unit),1))
            tit='OFF parasol';
        elseif ~isempty(find(datarun.cell_types{3}.cell_ids==datarun.cell_ids(unit),1))
            tit='ON midget';
        elseif ~isempty(find(datarun.cell_types{4}.cell_ids==datarun.cell_ids(unit),1))
            tit='OFF midget';
        elseif ~isempty(find(datarun.cell_types{5}.cell_ids==datarun.cell_ids(unit),1))
            tit='SBC';
        end
        
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        
        myWeightsAonly=myWeightsA(:,unit);
        
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[0 0 curcones(i)];
                siz=(curconesAbs(i))*relSizeCentr/(mean(curconesAbs)*30);
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*relSizeSur/(mean(curconesAbs)*30);
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[curconesA(i) 0 0];
                sizA=(curconesAbsA(i))*relSizeCentr/(mean(curconesAbsA)*30);
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*relSizeSur/(mean(curconesAbsA)*30);
            end

            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
            if myWeightsAonly(i)
                colA=[0 curconesA(i) curconesA(i)];
                plot(cca(i,1),cca(i,2),'o','markersize',sizA,'color',colA);                
            end
        end
        
        title([tit,', i=',int2str(unit),',  blue bayes 10 prior, red bayes 60 snl, cyan bayes 60 snl only'])
    end
end


% unique-found, bayes 60 snl
normalizedWeightsB=bw./repmat(max(bw),size(bw,1),1);
uniqueWeightsB=zeros(size(bw));
uniqueWeightsB(bUnique,:)=1;
myWeightsB=uniqueWeightsB&normalizedWeightsB>0.8;
myCellsB=find(sum(myWeightsB));

for unit=myCellsB
    
    if nnz(aw(:,unit))&&nnz(bw(:,unit))
        
        if ~isempty(find(datarun.cell_types{1}.cell_ids==datarun.cell_ids(unit),1))
            tit='ON parasol';
        elseif ~isempty(find(datarun.cell_types{2}.cell_ids==datarun.cell_ids(unit),1))
            tit='OFF parasol';
        elseif ~isempty(find(datarun.cell_types{3}.cell_ids==datarun.cell_ids(unit),1))
            tit='ON midget';
        elseif ~isempty(find(datarun.cell_types{4}.cell_ids==datarun.cell_ids(unit),1))
            tit='OFF midget';
        elseif ~isempty(find(datarun.cell_types{5}.cell_ids==datarun.cell_ids(unit),1))
            tit='SBC';
        end
        
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        
        myWeightsBonly=myWeightsB(:,unit);
       
        
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[0 0 curcones(i)];
                siz=(curconesAbs(i))*relSizeCentr/(mean(curconesAbs)*30);
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*relSizeSur/(mean(curconesAbs)*30);
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
            if myWeightsBonly(i)
                col=[0 curcones(i) curcones(i)];
                plot(ccb(i,1),ccb(i,2),'o','markersize',siz,'color',col);
            end
        end
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[curconesA(i) 0 0];
                sizA=(curconesAbsA(i))*relSizeCentr/(mean(curconesAbsA)*30);
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*relSizeSur/(mean(curconesAbsA)*30);
            end
            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
        end
        
        title([tit,', i=',int2str(unit),',  blue bayes 60 snl, red bayes, cyan bayes 60 snl only'])
    end
end








% SBC, only S cones run
sbc=datarun.cell_types{5}.cell_ids;
aw=datarunB5_Scones.cones.weights;
cca=datarunB5_Scones.cones.centers;

for tmp=1:length(sbc)
    
    unit=find(datarun.cell_ids==sbc(tmp),1);
    
    if nnz(aw(:,unit))
            
        figure
        sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
        colormap gray
        imagesc(sta)
        hold on
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[curconesA(i) 0 0];
                sizA=(curconesAbsA(i))*40/(mean(curconesAbsA)*30);
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*40/(mean(curconesAbsA)*30);
            end
            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
        end
        
        title(['SBC, i=',int2str(unit)])
    end
end



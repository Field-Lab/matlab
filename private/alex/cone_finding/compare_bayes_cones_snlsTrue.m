piece = '2012-09-13-2';
run = 'data009';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);

datarunB15 = datarun;
datarunB15 = load_cones(datarunB15,'bayes-msf_15');
datarunB15 = make_mosaic_struct(datarunB15);
datarunB15 = get_sta_fits_from_vision(datarunB15);  
datarunB15 = make_voronoi_masks(datarunB15);

datarunB10=datarun;
datarunB10 = load_cones(datarunB10,'bayes-msf_10');
datarunB10 = make_mosaic_struct(datarunB10);
datarunB10 = get_sta_fits_from_vision(datarunB10);  
datarunB10 = make_voronoi_masks(datarunB10);

datarunB5=datarun;
datarunB5 = load_cones(datarunB5,'bayes-msf_5');
datarunB5 = make_mosaic_struct(datarunB5);
datarunB5 = get_sta_fits_from_vision(datarunB5);  
datarunB5 = make_voronoi_masks(datarunB5);


datarunB45_snl = datarun;
datarunB45_snl = load_cones(datarunB45_snl,'bayes-msf_45.00-snlTrue_');
datarunB45_snl = make_mosaic_struct(datarunB45_snl);
datarunB45_snl = get_sta_fits_from_vision(datarunB45_snl);  
datarunB45_snl = make_voronoi_masks(datarunB45_snl);



datarunB55_snl = datarun;
datarunB55_snl = load_cones(datarunB55_snl,'bayes-msf_55.00-snlTrue_');
datarunB55_snl = make_mosaic_struct(datarunB55_snl);
datarunB55_snl = get_sta_fits_from_vision(datarunB55_snl);  
datarunB55_snl = make_voronoi_masks(datarunB55_snl);

datarunB60_snl = datarun;
datarunB60_snl = load_cones(datarunB60_snl,'bayes-msf_60.00-snlTrue_');
datarunB60_snl = make_mosaic_struct(datarunB60_snl);
datarunB60_snl = get_sta_fits_from_vision(datarunB60_snl);  
datarunB60_snl = make_voronoi_masks(datarunB60_snl);




a=datarunB15.cones.centers;
b=datarunB45_snl.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('bayes 15','bayes 45 SNL')




a=datarunB15.cones.centers;
b=datarunB55_snl.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('bayes 15','bayes 55 SNL')




a=datarunB15.cones.centers;
b=datarunB60_snl.cones.centers;
c=datarunB10.cones.centers;

figure
plot(c(:,1),c(:,2),'+g')
hold on
plot(a(:,1),a(:,2),'+r')
plot(b(:,1),b(:,2),'x')
legend('bayes 10','bayes 15','bayes 60 SNL')




%%%%%%%%%%
piece = '2012-09-13-2';
run = 'data009';

% define data path
datarun = load_data(['/Volumes/Analysis/' piece '/' run '/' run]);
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false,'load_ei',true);
datarun=load_data(datarun,opt);
datarun = load_params(datarun,struct('verbose',1));  
datarun = set_polarities(datarun);
datarun = load_sta(datarun,'load_sta',[]);

datarunB15 = datarun;
datarunB15 = load_cones(datarunB15,'bayes-msf_15');
datarunB15 = make_mosaic_struct(datarunB15);
datarunB15 = get_sta_fits_from_vision(datarunB15);  
datarunB15 = make_voronoi_masks(datarunB15);

datarunB60_snl = datarun;
datarunB60_snl = load_cones(datarunB60_snl,'bayes-msf_60.00-snlTrue_');
datarunB60_snl = make_mosaic_struct(datarunB60_snl);
datarunB60_snl = get_sta_fits_from_vision(datarunB60_snl);  
datarunB60_snl = make_voronoi_masks(datarunB60_snl);


datarun = load_sta(datarun,'load_sta','all');

% difference
a=datarunB15.cones.centers;
b=datarunB60_snl.cones.centers;

D=pdist2(a,b);
figure
m=min(D);
mrev=min(D')
hist(m,50)

figure
plot(a(:,1),a(:,2),'+')
hold on
plot(b(:,1),b(:,2),'xr')
legend('bayes 15','bayes 60 snl')



y=find(m<=1); % both for bayesian matrix
y2=find(m>1); % unique bayesian cones
x2=find(mrev>1); % unique local max cones

% y1=find(m>1 & m<=3.5);

figure
hold on
plot(a(x2,1),a(x2,2),'.r','markerSize',30) % unique local max cones
plot(b(y2,1),b(y2,2),'.b','markerSize',30) % unique bayesian cones
plot(b(y,1),b(y,2),'.c','markerSize',30) % both cones (from bayesian matrix)
legend('bayes 15','bayes 60 snl','both')
plot(b(:,1),b(:,2),'+')
plot(a(:,1),a(:,2),'xr')


aw=datarunB15.cones.weights;
bw=datarunB60_snl.cones.weights;



[x,y]=find(D<1);


figure
plot(b(y,1),b(y,2),'+')
hold on
plot(a(x,1),a(x,2),'xr')
title('only found by both methods')



% SBC
sbc=datarun.cell_types{5}.cell_ids;

cca=datarunB15.cones.centers;
ccb=datarunB60_snl.cones.centers;

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
        
        title(['SBC, i=',int2str(unit),',  blue bayes 60 snl, red bayes 15'])
    end
end



aUnique=setxor(x,1:length(a));
bUnique=setxor(y,1:length(b));


% unique-found, bayes 15
normalizedWeightsA=aw./repmat(max(aw),size(aw,1),1);
uniqueWeightsA=zeros(size(aw));
uniqueWeightsA(aUnique,:)=1;
myWeightsA=uniqueWeightsA&normalizedWeightsA>0.5;
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
        
        title([tit,', i=',int2str(unit),',  blue bayes 60 snl, red bayes 15, cyan bayes 15 only'])
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


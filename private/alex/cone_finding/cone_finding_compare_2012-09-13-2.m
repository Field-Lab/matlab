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

datarunM25 = datarun;
datarunM25 = load_cones(datarunM25,'ath25');
datarunM25 = make_mosaic_struct(datarunM25);
datarunM25 = get_sta_fits_from_vision(datarunM25);  
datarunM25 = make_voronoi_masks(datarunM25);

datarunM15 = datarun;
datarunM15 = load_cones(datarunM15,'ath15');
datarunM15 = make_mosaic_struct(datarunM15);
datarunM15 = get_sta_fits_from_vision(datarunM15);  
datarunM15 = make_voronoi_masks(datarunM15);

datarun = load_sta(datarun,'load_sta','all');

% local max difference
a=datarunM25.cones.centers;
b=datarunM15.cones.centers;

D=pdist2(a,b);
figure
m=min(D');
hist(m,50)

figure
plot(a(:,1),a(:,2),'+')
hold on
plot(b(:,1),b(:,2),'xr')
legend('local max 25','local max 15')


% bayes difference
a=datarunB5.cones.centers;
b=datarunB10.cones.centers;
c=datarunB15.cones.centers;

D=pdist2(b,c);
figure
m=min(D);
hist(m,50)

figure
plot(a(:,1),a(:,2),'+')
hold on
plot(b(:,1),b(:,2),'xr')
plot(c(:,1),c(:,2),'og')
legend('bayes 5','bayes 10','bayes 15')


% compare maximal bayes and maximal local max

a=datarunM15.cones.centers;
b=datarunB5.cones.centers;

figure
plot(a(:,1),a(:,2),'+r')
hold on
plot(b(:,1),b(:,2),'x')
legend('local max 15','bayes 5')

D=pdist2(a,b);
figure
m=min(D);
mrev=min(D');
hist(m,100)



y=find(m<=1); % both for bayesian matrix
y2=find(m>1); % unique bayesian cones
x2=find(mrev>1); % unique local max cones

% y1=find(m>1 & m<=3.5);

figure
hold on
plot(a(x2,1),a(x2,2),'.r','markerSize',30) % unique local max cones
plot(b(y2,1),b(y2,2),'.b','markerSize',30) % unique bayesian cones
plot(b(y,1),b(y,2),'.c','markerSize',30) % both cones (from bayesian matrix)
legend('localmax','bayes','both')
plot(b(:,1),b(:,2),'+')
plot(a(:,1),a(:,2),'xr')


aw=datarunM25.cones.weights;
bw=datarunB15.cones.weights;



[x,y]=find(D<1);


figure
plot(b(y,1),b(y,2),'+')
hold on
plot(a(x,1),a(x,2),'xr')
title('only found by both methods')



bwy=bw(y,:);
awx=aw(x,:);
abdiff=awx-bwy;

first_diff=nan(1,size(bwy,2));
weight_diff=nan(1,size(bwy,2));
for i=1:size(bwy,2)
    [plb,tmpb]=sort(bwy(:,i)/max(bwy(:,i)),'descend');
    [pla,tmpa]=sort(awx(:,i)/max(awx(:,i)),'descend');
    
    c=tmpb-tmpa;
    if nnz(c)
        first_diff(i)=find(c,1);
        weight_diff(i)=plb(1)-awx(tmpb(1),i)/max(awx(:,i));
%         weight_diff(i)=plb(first_diff(i))-pla(first_diff(i));
%         weight_diff(i)=sum(abs(plb(1:7)-awx(tmpb(1:7),i)/max(awx(:,i))));
%             weight_diff(i)=max(abs(plb(1:5)-awx(tmpb(1:5),i)/max(awx(:,i))));
%         figure
%         plot(plb,'-*')
%         hold on
%         plot(pla,'r-*')
% 
figure
plot(plb,'-*')
hold on
plot(awx(tmpb,i)/max(awx(:,i)),'r-*')
legend('bayes','localmax')

    end
end

sum(~isnan(weight_diff))

figure
hist(first_diff,1:max(first_diff))
figure
hist(weight_diff(abs(weight_diff)<0.1),30)

figure
hist(weight_diff(weight_diff<0.5),100)

hist(weight_diff,30)



plot(bwy(:,1))
hold on
plot(awx(:,1),'r')

plot(abdiff(:,1))






ccb=datarunB5.cones.centers(y,:);
cca=datarunM15.cones.centers(x,:);


figure
sta=squeeze(datarun.stas.stas{unit}(:,:,1,5));
colormap gray
imagesc(sta)
hold on
curcones=(bwy(:,unit)/max(bwy(:,unit))+1)/2;
curconesAbs=abs(bwy(:,unit)/max(bwy(:,unit)));
curconesA=(awx(:,unit)/max(awx(:,unit))+1)/2;
curconesAbsA=abs(awx(:,unit)/max(awx(:,unit)));
for i=1:length(ccb)
    
    if curcones(i)>0.5
        col=[curcones(i) 0 0];
        siz=(curconesAbs(i))*40;
    else
        col=[0 curcones(i) 0];
        siz=(curconesAbs(i))*50;
    end
    
    
    if curconesA(i)>0.5
        colA=[0 0 curconesA(i)];
    else
        colA=[0 curconesA(i) 0];
    end
    sizA=(curconesAbsA(i))*40;
    
    plot(ccb(i,1),ccb(i,2),'.','markersize',siz,'color',col);
    plot(cca(i,1),cca(i,2),'.','markersize',sizA,'color',colA);
end







ccb=datarunB5.cones.centers;
cca=datarunM15.cones.centers;

for unit=1:10%length(aw,2)
    if nnz(aw(:,unit))
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
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2;
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit)));
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[curcones(i) 0 0];
                siz=(curconesAbs(i))*40;
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*50;
            end
            if isempty(find(y==i, 1))
                col(3)=curcones(i);
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[0 0 curconesA(i)];
                sizA=(curconesAbsA(i))*40;
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*50;
            end
            if isempty(find(x==i, 1))
                colA(2)=curconesA(i);
            end
            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
        end
        
        title([tit,', i=',int2str(unit),',  blue bayes, red local max'])
    end
end






ccb=datarunB15.cones.centers;
cca=datarunM25.cones.centers;

for unit=1:10%length(aw,2)
    if nnz(aw(:,unit))
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
        curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2;
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit)));
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                col=[curcones(i) 0 0];
                siz=(curconesAbs(i))*40;
            else
                col=[0 curcones(i) 0];
                siz=(curconesAbs(i))*50;
            end
            if isempty(find(y==i, 1))
                col(3)=curcones(i);
            end
            plot(ccb(i,1),ccb(i,2),'+','markersize',siz,'color',col);
        end
        
        for i=1:length(cca)
            if curconesA(i)>0.5
                colA=[0 0 curconesA(i)];
                sizA=(curconesAbsA(i))*40;
            else
                colA=[0 curconesA(i) 0];
                sizA=(curconesAbsA(i))*50;
            end
            if isempty(find(x==i, 1))
                colA(2)=curconesA(i);
            end
            plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
        end
        
        title([tit,', i=',int2str(unit),',  blue bayes, red local max'])
    end
end




% SBC
sbc=datarun.cell_types{5}.cell_ids;

ccb=datarunB5.cones.centers;
cca=datarunM15.cones.centers;

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
        
        title(['SBC, i=',int2str(unit),',  blue bayes, red local max'])
    end
end



aUnique=setxor(x,1:length(a));
bUnique=setxor(y,1:length(b));


% unique-found, local max
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
        
        title([tit,', i=',int2str(unit),',  blue bayes, red local max, cyan local max only'])
    end
end


% unique-found, bayes
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
        
        title([tit,', i=',int2str(unit),',  blue bayes, red local max, cyan bayes only'])
    end
end




figure
units=[3 4 8];
tmp=-squeeze(datarun.stas.stas{units(1)}(:,:,1,5));
tmp1=-squeeze(datarun.stas.stas{units(2)}(:,:,1,5));
tmp2=-squeeze(datarun.stas.stas{units(3)}(:,:,1,5));

sta(:,:,1)=(tmp/max(abs(tmp(:)))+1)/2;
sta(:,:,2)=(tmp1/max(abs(tmp1(:)))+1)/2;
sta(:,:,3)=(tmp2/max(abs(tmp2(:)))+1)/2;

imagesc(sta)
col='rgbmyc';
mark='ooo';
relSizeSur=50;
relSizeCentr=20;

hold on
cnt=1;
plot(ccb(:,1),ccb(:,2),'x','markersize',10,'color','k');
plot(cca(:,1),cca(:,2),'+','markersize',10,'color','k');
for unit=units
    curcones=(bw(:,unit)/max(abs(bw(:,unit)))+1)/2; % normalized to mean=0.5 and max=1
        curconesAbs=abs(bw(:,unit)/max(bw(:,unit))); % normalized to max 1, mean 0, and then absolute
        curconesA=(aw(:,unit)/max(abs(aw(:,unit)))+1)/2;
        curconesAbsA=abs(aw(:,unit)/max(aw(:,unit)));
        
        for i=1:length(ccb)
            
            if curcones(i)>0.5
                ccol=col(cnt);
                siz=(curconesAbs(i))*relSizeCentr/(mean(curconesAbs)*30);
                marks=mark(cnt);
            else
                ccol=col(cnt+3);
                siz=(curconesAbs(i))*relSizeSur/(mean(curconesAbs)*30);
                marks='p';
            end
            plot(ccb(i,1),ccb(i,2),marks,'markersize',siz,'color',ccol);
        end
        cnt=cnt+1;
        
%         for i=1:length(cca)
%             if curconesA(i)>0.5
%                 colA=[curconesA(i) 0 0];
%                 sizA=(curconesAbsA(i))*relSizeCentr/(mean(curconesAbsA)*30);
%             else
%                 colA=[0 curconesA(i) 0];
%                 sizA=(curconesAbsA(i))*relSizeSur/(mean(curconesAbsA)*30);
%             end
%             plot(cca(i,1),cca(i,2),'x','markersize',sizA,'color',colA);
%         end
        

end


figure;
tmp=-squeeze(datarun.stas.stas{units(1)}(:,:,1,5));
[tmp,ic]=sort(tmp(:)/max(abs(tmp(:))),'descend');
plot(tmp.^2,'-o')

figure
plot(tmp(1:end/2).^2,'-o')
hold on
plot(tmp(end:-1:end/2).^2,'-or')
hold on
tmp1=-squeeze(datarun.stas.stas{units(2)}(:,:,1,5));
tmp1=tmp1(:)/max(abs(tmp1(:)));
plot(tmp1(ic),'-xr')
tmp2=-squeeze(datarun.stas.stas{units(3)}(:,:,1,5));
tmp2=tmp2(:)/max(abs(tmp2(:)));
plot(tmp2(ic),'-+g')



tmp=-squeeze(datarun.stas.stas{units(1)}(:,:,1,:));
tmp=reshape(tmp, 300*300,6);
a=std(tmp,0,2);
figure;
plot(a(ic),'-x')
plot(sort(a,'descend'),'-o')
m=sort(a,'descend');
m=m(1:60);


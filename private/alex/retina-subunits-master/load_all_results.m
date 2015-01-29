path2data='/Volumes/Analysis/2011-10-25-5/subunits/data001-0/anal_orig/subunit/';
a=dir([path2data,'*greedy*.mat']);

clear all_res

for i=1:length(a)
    load([path2data,a(i).name])
    all_res(i)=res;
end

clear mean_w mean_knots mean_Mspline aa bb cc
aa=zeros(6,length(all_res));
bb=zeros(length(all_res),8);
cc=zeros(length(all_res),28,6);
for i=1:length(all_res)
    
    aa(1:6,i)=all_res(i).fit_SUB.f.w;
    bb(i,1:8)=all_res(i).fit_SUB.f.knots;
    cc(i,:,:)=all_res(i).fit_SUB.f.Mspline;
    
end
mean_w=mean(aa');
mean_knots=mean(bb);
mean_Mspline=squeeze(mean(cc));

save('/Volumes/Analysis/2011-10-25-5/subunits/data001-0/orig_fit_results.mat','all_res','mean_w','mean_knots','mean_Mspline')

    
    

dataset='2011-10-25-5/data001-0';
cellTypeName='OFF midget';
workstation='bertha';


dat = loadData(workstation,dataset);
disp('now running stc')
batchAnal(workstation,dataset,cellTypeName,'stc',1,0,0.33);
disp('now running subunits')
tic
batchAnal(workstation,dataset,cellTypeName,'subunit-local',1,0,0.33)
toc
disp('now running figures')
batchAnal(workstation,dataset,cellTypeName,'subunit-local',0,1,0.33)




path2data='/Volumes/Analysis/2011-10-25-5/subunits/data001-0/anal_fixed_F/subunit/';
a=dir([path2data,'*greedy*.mat']);

clear all_res

for i=1:length(a)
    load([path2data,a(i).name])
    all_res(i)=res;
end

save('/Volumes/Analysis/2011-10-25-5/subunits/data001-0/fixed_fit_results.mat','all_res','mean_w','mean_knots','mean_Mspline')


path2data='/Volumes/Analysis/2011-10-25-5/subunits/data001-0/anal_orig/subunit/';
names=dir([path2data,'*greedy*.mat']);

clear all_res_orig

for i=1:length(names)
    load([path2data,names(i).name])
    all_res_orig(i)=res;
    tmp=regexp(names(i).name,'-');
    namings(i)=str2num(names(i).name(1:tmp(1)-1));
end


clear a b c
for i=1:length(all_res)
    
    a(i)=all_res(i).out_SUB.r2;
    b(i)=all_res_orig(i).out_SUB.r2;
    c(i)=all_res(i).out_LN.r2;
end

figure
subplot(2,1,1)
plot(a, '-*')
hold on
plot(b,'-*r')
plot(c,'-*k')
legend('fixed F','orig', 'LN')
title('r2 values, 62 OFF midgets')
subplot(2,1,2)
plot((a-b)./(a+b), '-*k')
title('r2 difference, (fixed-orig)/(fixed+orig)')
line([1 62],[0 0], 'color','k')


dat = loadData('bertha',dataset);

loadType = 1;
numCells = getCellTypeNum(dat.cellTypes,'OFF midget');
cellIds = 1:numCells;
celldat = [];
celldat.cellType = 'OFF midget';
celldat.loadType = loadType;
celldat = getDefaultOpts(celldat);
clear a b
cnt=1;
for i=1:numCells
    
    celldat.cellNum = cellIds(i);
    celldat.rgcId = cellIds(i);    
    celldat.percent = 0.33;
    
    [train test celldat] = loadCellData(dat,celldat,0);
    
    tmp=find(namings==celldat.rgcId);
    if ~isempty(tmp)
        
        diffVals = (all_res(tmp).out_SUB.Z_t - all_res(tmp).out_LN.Z_t).^2;
        inds = diffVals > prctile(diffVals,80) & diffVals < prctile(diffVals,100);
        %out_SING_SEL = evalFit(test,fit_SING,inds);
        out_SUB_SEL = evalFit(test,all_res(tmp).fit_SUB,inds);
        out_LN_SEL =  evalFit(test,all_res(tmp).fit_LN,inds);
        
        R2_SEL = [out_LN_SEL.r2; out_SUB_SEL.r2];
        
        a(1:2,cnt)=R2_SEL;
        
        
        
        diffVals = (all_res_orig(tmp).out_SUB.Z_t - all_res_orig(tmp).out_LN.Z_t).^2;
        inds = diffVals > prctile(diffVals,80) & diffVals < prctile(diffVals,100);
        %out_SING_SEL = evalFit(test,fit_SING,inds);
        out_SUB_SEL = evalFit(test,all_res_orig(tmp).fit_SUB,inds);
        out_LN_SEL =  evalFit(test,all_res_orig(tmp).fit_LN,inds);
        
        R2_SEL = [out_LN_SEL.r2; out_SUB_SEL.r2];
        
        b(1:2,cnt)=R2_SEL;
        cnt=cnt+1;
    end

end


figure
subplot(2,2,1)
plot(a', '-*')
title('subset fixed')
legend('LN', 'sub')
axis([1 62 -0.05 0.4])

subplot(2,2,2)
plot(b', '-*')
title('subset orig')
legend('LN', 'sub')
axis([1 62 -0.05 0.4])

subplot(2,2,3)
plot(diff(a), '-*k')
title('subset fixed r2 difference (LN-sub)')
axis([1 62 0 0.3])

subplot(2,2,4)
plot(diff(b), '-*k')
title('subset orig r2 difference (LN-sub)')
axis([1 62 0 0.3])

figure
plot(diff(a)-diff(b), '-*k')
title('subset r2 difference (LN-sub), FIXED-ORIG')
axis tight




clear a b c tt
cnt=1;
for i=1:length(all_res)
    
    a=all_res(i).fit_SUB.I_sc;
    b=all_res_orig(i).fit_SUB.I_sc;
    if ~prod(size(a)==size(b)) || nnz(a~=b)
        tt(cnt)=namings(i);
        cnt=cnt+1;
    end
end

path  = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/';
% load([path,'/Off_par2.mat']);
% load([path,'/Off_par.mat']);
load([path,'/On_par2.mat']);
cell_list = [121,586,650,826,841,1205,1276,1352,1426,1502,1772,1921,2101,2312,2313,2641,3152,3226,3647,3799,4021,4216,4366,4456,4711,4756,5131,5132,5134,5137,5431,5777,5778,5783,5914,5972,6257,6721,7068,7171,7291,7381,7607,7756];

cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cellsChoose([25,26,31,23])=1; 
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;

%% *
path = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/'
load([path,'/Off_type1.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cell_choose_num = [3,7,8];
cellsChoose(cell_choose_num)=1; % [25,26,31,23] or [25]
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

%% 
path = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_08_27_6/'
load([path,'/On_large_type2.mat']);
binnedSpikeResponses_coll = Y;
total_mask_log = totalMaskAccept_log;

cells = double(cells);
cellsChoose = (cells == 5931) | (cells == 5971) | (cells ==5941) | (cells == 6196) | (cells == 6109);

cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
ttf = ttf_avg;

%%
ifit=0;
for fitNnum=1:1
for Ns=10%[6,5,4,7,8,9];
%fitASM_pop = fitASM_EM_Population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask);
ifit=ifit+1
lam_start = 0.1;
gamma1=0.0000;
gamma2 = 0;
initVal=[];
%initVal.K = maskedMovdd(mask,:)*binnedSpikeResponses_coll(cellsChoose,:)';
%initVal.B = diag(ones(Ns,1));
[fitASM_pop,fval] = fitASM_EM_Population_sparse_split_admm(maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start,initVal);
B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_08_27_6/On_large1_fit_%d.mat',ifit),'Ns','fitASM_pop','fval','mask','gamma1','gamma2','lam_start','initVal','cellsChoose');
pause(0.5);
end
end
        
rho=1;
lambda=0.1;
Xref = maskedMovdd(mask,1:1000);
%[Xdecode,errMap] = decode_ASM_population3D(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:100),mask,ttf,rho,lambda,1000*maskedMovdd(mask,1:100),B_use)
[Xdecode,errMap] = decode_ASM_population(fitASM_pop,binnedSpikeResponses_coll(cellsChoose,1:1000),mask,ttf,rho,lambda,Xref,B_use)
 


[fitASM_pop,fval] = fitASM_EM_Population_sparse(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start);
[fitASM_pop,fval] = fitASM_EM_Population_sparse_split(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start);
[fitASM_pop] = fitASM_EM_Population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask);

isave=90;
for icell=cell_choose_num
    mask_Cell = sum(total_mask_log(:,icell),2)~=0; 
    for Ns = [2,3,4,6,7,8,9];
    
    fitASM_log = cell(50,1);
    fval_log = cell(50,1);
    
    for ifit=1:50
[fitGMLM,fval] = fitGMLM_EM_bias(binnedSpikeResponses_coll(icell,:)',maskedMovdd(mask_Cell,:),sum(mask_Cell),Ns,1);
    fitASM_log{ifit} =fitGMLM;
    fval_log{ifit} =fval;
    save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_Type1_fit_%d.mat',isave),'Ns','fitASM_log','fval_log','mask_Cell','icell','ttf','ifit');

    end
    isave=isave+1;
    end
end


%plotSU(fitGMLM.Linear.filter{2},mask)
plotSU(fitGMLM.Linear.filter{1},mask)

%init_value.K =  gpuArray(2*(rand(sum(mask),Ns)-0.5));
%init_value.B = gpuArray(2*(rand(Ns,sum(cellsChoose))-0.5));

%[fitASM,fval_log] = fitASM_sgd_population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),Ns,mask,gamma1,gamma2,lam_start, init_value);

%%

ifit=18;
load(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/Off_type1_fit_%d.mat',ifit));
[B_use,h]= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));
% suptitle({sprintf('nSU = %d ',Ns)});
hgexport(h,sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2008_05_13_3/results_fit_%d.eps',ifit));


path  = '/Volumes/Lab/Users/bhaishahster/GMLM_fits/pc2012_08_09_3/';
load([path,'/Off_par2.mat']);
load([path,'/Off_par.mat']);

%%
cellsChoose = zeros(size(binnedSpikeResponses_coll,1),1);
cellsChoose([1,2,3,4,5])=1;
cellsChoose=logical(cellsChoose);
mask = sum(total_mask_log(:,cellsChoose),2)~=0;
% fitASM_pop = fitASM_EM_Population(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),10,mask);


fitASM_pop = fitASM_EM_Population_sparse(1000*maskedMovdd,binnedSpikeResponses_coll(cellsChoose,:),10,mask,gamma1,gamma2);
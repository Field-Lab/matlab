% Test population 
% 2 midget cells share a sub-unit..

nCones = 14;
nSU = 11;
nCells = 2;

SU_cone_mat = zeros(nSU,nCones);
SU_cone_mat(1,[1,2])=1;
SU_cone_mat(11,[13,14])=1;
SU_cone_mat(7,[8,9])=1;
SU_cone_mat(2,3)=1;
SU_cone_mat(3,4)=1;
SU_cone_mat(4,5)=1;
SU_cone_mat(5,6)=1;
SU_cone_mat(6,7)=1;
SU_cone_mat(8,10)=1;
SU_cone_mat(9,11)=1;
SU_cone_mat(10,12)=1;
SU_cone_wt = SU_cone_mat;

cell_su_mat = zeros(nCells,nSU);
cell_su_mat(1,[1,2,3,4,5,6,7])=1;
cell_su_mat(2,[7,8,9,10,11])=1;
cell_su_wt = cell_su_mat*2;

% f =@(x)
 %f=@(x)exp(x/5);
 f = @(x) max(x/2,0);
g=@(x)x;

figure;
subplot(1,2,1);
imagesc(SU_cone_mat);
xlabel('Cones');
ylabel('SU');

subplot(1,2,2);
imagesc(cell_su_wt);
xlabel('SU');
ylabel('Cells');
%% 
T = 120*60*60;
movie_c = 2*randn(nCones,T);
su_inp = SU_cone_wt*movie_c;
cell_op = g(cell_su_wt*f(su_inp));
dt=1/120;

figure;
[x,n] = hist(su_inp(1,:),30);
plotyy(n,x,n,f(n));


resp = poissrnd(cell_op*dt);
firingRate = sum(resp,2)/(60*60);

%% fit model
icnt=15;
for Ns=[14,8,11];
    icnt=icnt+1;
mask = logical(ones(nCones,1));
lam_start = 0.1;
gamma1=0.001; % 0.0005
gamma2 = 0;
initVal=[];
%initVal.K = maskedMovdd(mask,:)*binnedSpikeResponses_coll(cellsChoose,:)';
%initVal.B = diag(ones(Ns,1));
fitASM_pop_log=cell(50,1);
fval_log = cell(50,1);

for ifit=1:50
[fitASM_pop,fval] = fitASM_EM_Population_sparse_split_admm(movie_c,resp,Ns,mask,gamma1,gamma2,lam_start,initVal);
fitASM_pop_log{ifit} = fitASM_pop;
fval_log{ifit} = fval;
%save('/Volumes/Lab/Users/bhaishahster/GMLM_fits/ts_population.mat','fitASM_pop_log','fval_log','movie_c','resp','Ns','gamma1','gamma2','mask','lam_start');
save(sprintf('/Volumes/Lab/Users/bhaishahster/GMLM_fits/ts_population%d.mat',icnt),'fitASM_pop_log','fval_log','movie_c','resp','Ns','gamma1','gamma2','mask','lam_start','f','g');

end
end

for ifit=1:50
    close all
spider((fitASM_pop_log{ifit}.K),'fit');
pause
end
%B_use= plotSU_withcells(fitASM_pop.K,mask,total_mask_log(:,cellsChoose),exp(fitASM_pop.B));


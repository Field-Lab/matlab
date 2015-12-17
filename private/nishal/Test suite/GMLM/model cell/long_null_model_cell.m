%% Explore long null stimulus on test suite .. 

%% make simple model cell, generate responses to WN and calculate STA .. from Sub_unit_test1.m

%% Long null - Can we do fast fitting using null stimulus ? 
%% Generate Original and null stimulus and responses , 30 min! 

movieLen=120*60*40;
stas{1}=zeros(size(STA,1),size(STA,2),3,size(STA,3));
stas{1}(:,:,1,:)=STA;
stas{1}(:,:,2,:)=STA;
stas{1}(:,:,3,:)=STA;
cell_params.STAlen=14;
[new_stas,mask,CellMasks]=clipSTAs(stas,cell_params)
Filtdim1=size(STA,1);Filtdim2=size(STA,2);Filtlen=30;
%null_compute_subUnit_test
 model_spatial_nulling
 
 
mov_new2 = mov_new2*1.2; % Have some gain!!
 
nTrials=1;
analyse_null_subUnit_ts

%% Initialize stuff for next tests!
% Compare actual and predicted for original and null movie.
 
% get SUs
corr_su = zeros(length(subunits),sum(mask(:)));
for isu=1:length(subunits)
xx= squeeze(subunits{isu}(:,:,6));
xx = xx(mask>0);
corr_su(isu,:)=xx;
end
comp_true_SU_str.compare=1;
comp_true_SU_str.su = gather(corr_su);
comp_true_SU_str.proj = gather(eye(size(corr_su,2)));



mov=movOrig(:,:,Filtlen:movie_new_len+Filtlen-1);
 maskedMov_orig= filterMov(mov,mask,squeeze(tf));
 maskedMov2_orig=[maskedMov_orig;ones(1,size(maskedMov_orig,2))];
nSU=4;
  initialFilters = 2*(rand(size(maskedMov_orig,1)*nSU,1)-0.5);
  
  dataStartTest = 30*120*60 + 120 +1;
  
%% initialize the filters in null space!

 % Compare actual and predicted for Null movie.
mov=movNull(:,:,Filtlen:movie_new_len+Filtlen-1);
maskedMov_null= filterMov(mov,mask,squeeze(tf));
maskedMov2_null=[maskedMov_null;ones(1,size(maskedMov_null,2))];

% make test structure for null movies
testStructNull.test=1;
testStructNull.mov =gpuArray(maskedMov_null(:,dataStartTest:end));
testStructNull.binnedResponses =gpuArray(binnedResponseNull(dataStartTest+1:end,1));

sta = maskedMov_orig*binnedResponseOrig/sum(binnedResponseOrig);
sta=sta';
[u,s,v]=svd(sta,'econ');
stainv=v*(s^-1)*u';
Prj_null = eye(size(stainv,1),size(sta,2))-(stainv*sta);
Prj_sta = (stainv*sta);

comp_true_SU_str.compare=1;
comp_true_SU_str.su = gather(corr_su);
comp_true_SU_str.proj = Prj_null;



init =  initialFilters;
for ifilt=1:nSU
    iidx=((ifilt-1)*size(maskedMov_null,1)+1:ifilt*size(maskedMov_null,1)) ;
init(iidx)= Prj_null*init(iidx);
end
interval=1;
 icnt=0;
 
for useMinNull = [0.0312,0.0625,0.125,0.25,0.5,1,2.5,5,7,10,15,20,30]
dataLenNull = useMinNull*120*60 + 120;

[fitGMLM,fval_initnull,fval_initnull_testnull,d_l_initnull_testnull] = fitGMLM_EM_power2(binnedResponseNull(1:dataLenNull,1),maskedMov_null(:,1:dataLenNull),size(maskedMov_null,1),nSU,interval,2,'initialFilters',init,'testStruct',testStructNull,'compare_true_SU',comp_true_SU_str);  
fitGMLMNull_initnull=fitGMLM;
%filterNull_initnull = showfitGMLM(fitGMLMNull_initnull,sprintf('# spks : %d Quad, null movie, random initialize in null space',sum(binnedResponseNull(1:dataLenNull,1))),mask)

icnt=icnt+1;
data_1(icnt).fitGMLM = fitGMLM;
data_1(icnt).fval = fval_initnull;
data_1(icnt).fval_test = fval_initnull_testnull;
data_1(icnt).d_true = d_l_initnull_testnull;
 data_1(icnt).useMin = useMinNull;
end
%% initialize and fit filters with lower dimensional parametrization null space!


 % Compare actual and predicted for Null movie.
mov=movNull(:,:,Filtlen:movie_new_len+Filtlen-1);
maskedMov_null= filterMov(mov,mask,squeeze(tf));
maskedMov2_null=[maskedMov_null;ones(1,size(maskedMov_null,2))];



sta = maskedMov_orig*binnedResponseOrig/sum(binnedResponseOrig);
sta=sta';
[u,s,v]=svd(sta,'econ');
stainv=v*(s^-1)*u';
Prj_null = eye(size(stainv,1),size(sta,2))-(stainv*sta);
Prj_sta = (stainv*sta);
[u,s,v] = svd(Prj_null,'econ');
P_ld = sqrt(s(1:end-1,1:end-1))*v(:,1:end-1)';
P_ld_rev = u(:,1:end-1)*sqrt(s(1:end-1,1:end-1));



comp_true_SU_str.compare=1;
comp_true_SU_str.su = gather(P_ld*corr_su')';
comp_true_SU_str.proj = Prj_null*P_ld_rev;
    
% make test structure for null movies
testStructNull.test=1;
testStructNull.mov =gpuArray(P_ld*maskedMov_null(:,dataStartTest+1:end));
testStructNull.binnedResponses =gpuArray(binnedResponseNull(dataStartTest+1:end,1));


init =  initialFilters;
init2 =  zeros((size(maskedMov_null,1)-1)*nSU,1);
for ifilt=1:nSU
    iidx=((ifilt-1)*size(maskedMov_null,1)+1:ifilt*size(maskedMov_null,1)) ;
    iidx2=((ifilt-1)*(size(maskedMov_null,1)-1)+1:ifilt*(size(maskedMov_null,1)-1)) ;
init2(iidx2)= P_ld*init(iidx);
end
maskedMov_null_ld = P_ld*maskedMov_null;

interval=1;
 icnt=0;
 data2=struct([]);
for useMinNull = [0.0312,0.0625,0.125,0.25,0.5,1,2.5,5,7,10,15,20,30]
dataLenNull = useMinNull*120*60 + 120;

[fitGMLM,fval_initnull_fitld,fval_initnull_fitld_test_null,d_l_fitld_test_null] = fitGMLM_EM_power2(binnedResponseNull(1:dataLenNull,1),maskedMov_null_ld(:,1:dataLenNull),size(maskedMov_null_ld,1),nSU,interval,2,'initialFilters',init2,'testStruct',testStructNull,'compare_true_SU',comp_true_SU_str);  
fitGMLMNull_initnull_fitld=fitGMLM;
%filterNull_initnull_fitld = showfitGMLM(fitGMLMNull_initnull_fitld,sprintf('# spks : %d Quad, null movie, random initialize and fit in null space space',sum(binnedResponseNull(1:dataLenNull,1))),mask,'preproj',P_ld_rev);
icnt=icnt+1;
data_2(icnt).fitGMLM = fitGMLM;
data_2(icnt).fval = fval_initnull_fitld;
data_2(icnt).fval_test = fval_initnull_fitld_test_null;
data_2(icnt).d_true = d_l_fitld_test_null;
 data_2(icnt).useMin = useMinNull;
end

%% Fit to null stimulus 

 % Compare actual and predicted for Null movie.
mov=movNull(:,:,Filtlen:movie_new_len+Filtlen-1);
maskedMov_null= filterMov(mov,mask,squeeze(tf));
maskedMov2_null=[maskedMov_null;ones(1,size(maskedMov_null,2))];

% make test structure for null movies
testStructNull.test=1;
testStructNull.mov =gpuArray(maskedMov_null(:,dataStartTest+1:end));
testStructNull.binnedResponses =gpuArray(binnedResponseNull(dataStartTest+1:end,1));

% compare true SU
comp_true_SU_str.compare=1;
comp_true_SU_str.su = gather(corr_su);
comp_true_SU_str.proj =Prj_null;

% Fit model 

  interval=1;
 icnt=0;
 data_3=struct([]);
for useMinNull = [0.0312,0.0625,0.125,0.25,0.5,1,2.5,5,7,10,15,20,30];
dataLenNull = useMinNull*120*60 + 120;

  [fitGMLM,fval_null,fval_null_test_null,d_l_null] = fitGMLM_EM_power2(binnedResponseNull(1:dataLenNull,1),maskedMov_null(:,1:dataLenNull),size(maskedMov_null,1),nSU,interval,2,'initialFilters',initialFilters,'testStruct',testStructNull,'compare_true_SU',comp_true_SU_str);  
  fitGMLMNull=fitGMLM;
  %filterNull = showfitGMLM(fitGMLMNull,sprintf('# spks : %d Quad, null movie, random initialize',sum(binnedResponseNull(1:dataLenNull,1))),mask)

  icnt=icnt+1;
data_3(icnt).fitGMLM = fitGMLM;
data_3(icnt).fval = fval_null;
data_3(icnt).fval_test = fval_null_test_null;
data_3(icnt).d_true =d_l_null;
 data_3(icnt).useMin = useMinNull;
end

%% Original movie


 % Compare actual and predicted for Null movie.
mov=movNull(:,:,Filtlen:movie_new_len+Filtlen-1);
maskedMov_null= filterMov(mov,mask,squeeze(tf));
maskedMov2_null=[maskedMov_null;ones(1,size(maskedMov_null,2))];

% make test structure for null movies
testStructNull.test=1;
testStructNull.mov =gpuArray(maskedMov_null(:,dataStartTest+1:end));
testStructNull.binnedResponses =gpuArray(binnedResponseNull(dataStartTest+1:end,1));

comp_true_SU_str.compare=1;
comp_true_SU_str.su = gather(corr_su);
comp_true_SU_str.proj = Prj_null;

mov=movOrig(:,:,Filtlen:movie_new_len+Filtlen-1);
 maskedMov_orig= filterMov(mov,mask,squeeze(tf));
 maskedMov2_orig=[maskedMov_orig;ones(1,size(maskedMov_orig,2))];

 % Fit model
   interval=1;
    data_4=struct([]);
    icnt=0;
 for useMinOrig = [0.0312,0.0625,0.125,0.25,0.5,1,2.5,5,7,10,15,20,30];
dataLenOrig = useMinOrig*120*60 + 120;

  [fitGMLM,fval_log_orig,fval_log_orig_on_null,d_l_orig] = fitGMLM_EM_power2(binnedResponseOrig(1:dataLenOrig,1),maskedMov_orig(:,1:dataLenOrig),size(maskedMov_orig,1),nSU,interval,2,'initialFilters',initialFilters,'testStruct',testStructNull,'compare_true_SU',comp_true_SU_str);  
  fitGMLMOrig=fitGMLM;
 % filterOrig = showfitGMLM(fitGMLMOrig,sprintf('# spks : %d Quad, WN movie',sum(binnedResponseOrig(1:dataLenOrig,1))),mask);
 %  filterOrigproj_Null = showfitGMLM(fitGMLMOrig,sprintf('# spks : %d Quad, WN movie, projected into null space',sum(binnedResponseOrig(1:dataLenOrig,1))),mask,'proj',gather(P_null));
   
 icnt=icnt+1;
data_4(icnt).fitGMLM = fitGMLM;
data_4(icnt).fval = fval_log_orig;
data_4(icnt).fval_test = fval_log_orig_on_null;
data_4(icnt).d_true =d_l_orig;
 data_4(icnt).useMin = useMinOrig;
end

   %% plots
figure;
plot(fval_null_test_null);
hold on;
plot(fval_initnull_testnull);
hold on;
plot(fval_initnull_fitld_test_null);
hold on;
plot(fval_log_orig_on_null);
legend('fval null','fval initnull','fval initnull fit ld','fval orig on null');




figure;
  plot(d_l_null)
hold on;
plot(d_l_initnull_testnull)
hold on;
plot(d_l_fitld_test_null)
hold on;
plot(d_l_orig);
legend('fval null','fval initnull','fval initnull fit ld','fval orig on null');

%% 

[data_len1,fval_test1,d1] = collect_ends(data_1);
[data_len2,fval_test2,d2] = collect_ends(data_2);
[data_len3,fval_test3,d3] = collect_ends(data_3);
[data_len4,fval_test4,d4] = collect_ends(data_4);

figure;
plot(data_len1,fval_test1);
hold on;
plot(data_len2,fval_test2);
hold on;
plot(data_len3,fval_test3);
hold on;
plot(data_len4,fval_test4);
legend('init null','init null fit null','random init but fit using null data','orig fit');



figure;
plot(data_len1,d1);
hold on;
plot(data_len2,d2);
hold on;
plot(data_len3,d3);
hold on;
plot(data_len4,d4);
legend('init null','init null fit null','random init but fit using null data','orig fit');


%% combined - mix original and null !

fracNull_list = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
mov_len_list = [0.5,1,2,3,5,7,10,20,30];
filt_dist = zeros(length(mov_len_list),length(fracNull_list),20);
fitLog = cell(length(mov_len_list),length(fracNull_list),20);
for ifracNull = 9:length(fracNull_list)
for imov_len = 1:length(mov_len_list)

   mov_len =  mov_len_list(imov_len)*60*120+1;
fracNull =fracNull_list(ifracNull);  
fracOrig = 1-fracNull;
[fracNull,imov_len]
mov1=movOrig(:,:,Filtlen:mov_len*fracOrig+Filtlen-1);
mov2=movNull(:,:,Filtlen:mov_len*fracNull+Filtlen-1); 
mov = cat(3,mov1,mov2);
binnedResponses = [binnedResponseOrig(1:mov_len*fracOrig);binnedResponseNull(1:mov_len*fracNull)];
maskedMov= filterMov(mov,mask,squeeze(tf));
  
  interval=1;
 
  nSU=4;
  for niter = 1:20
  [fitGMLM,output] = fitGMLM_EM_accelerate_power2(binnedResponses,maskedMov,7,nSU,interval,2);  
  fitGMLM.fval = output;
  fitGMLMcombined=fitGMLM;
  filterCombined = showfitGMLM(fitGMLMcombined,sprintf('# spks : %d Quad, null movie',sum(binnedResponses)),mask);
fitLog{imov_len,ifracNull,niter} = fitGMLMcombined;
 filt_dist(imov_len,ifracNull,niter)= find_dist_filters(filter_true,filterCombined)
 close all
  end
end
end
  
   filterCombined = showfitGMLM(fitLog{end,end},sprintf('# spks : %d Quad, null movie',sum(binnedResponses)),mask);
  find_dist_filters(filter_true,filterOrig)
  
  figure;
  surf(fracNull_list',mov_len_list',log(mean(filt_dist,3)));
  %contour(fracNull_list',mov_len_list',log(filt_dist),20);
  xlabel('fraction of Null');
  ylabel('Minutes of recording ');
  
  
  % prediction quality
  
movieLen=120*5*60;
null_compute_subUnit_test
 
mov_new2 = mov_new2*3; % Have some gain!!
 
nTrials=1;
analyse_null_subUnit_ts

  
  
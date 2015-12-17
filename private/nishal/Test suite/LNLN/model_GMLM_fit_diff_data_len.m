

function [fitGMLM_log,mse_data] = model_GMLM_fit_diff_data_len(movie,response,model,nSU,initialData,mask2,P_mat)


 maskedMov= filterMov(movie,mask2,squeeze(model.ttf));
 maskedMov2=[maskedMov;ones(1,size(maskedMov,2))];

 binnedResponses = response';
 filteredStimDim =size(maskedMov,1);
 
 %  EM like Max Expected Likelihood .. 
 interval=1;
 
 fitGMLM_log = cell(50,1);
 fval_log = zeros(50,1);

 %nSU = 6;
 
initialFilters =initialData.initialFilters;
initalbias = initialData.initalbias;

ifit = 0;
frac_list = 1%0.1:0.1:1;
  allidx = randperm(length(binnedResponses));
for frac=frac_list
    ifit= ifit+1;
trainData=allidx(1:floor(length(binnedResponses)*frac));
maskedMov_input = maskedMov(:,trainData);
binnedResponses_input = binnedResponses(trainData);
[fitGMLM,f_val] = fitGMLM_EM_bias_initPoint(binnedResponses_input,maskedMov_input,filteredStimDim,nSU,interval,initialFilters,initalbias); 

fitGMLM_log{ifit} = fitGMLM;
fval_log(ifit) = f_val;
end


%% show SUs - original
szz =[1,1;
    1,2;
    2,2;
    2,2;
    2,3;
    2,3;
    3,3;
    3,3;
    3,3;
    3,4;
    3,4;
    3,4;
    4,4;
    4,4;
    4,4];

for ifit=1:length(frac_list)
 fitGMLM = fitGMLM_log{ifit};   
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

h=figure;
for ifilt=1:nSU
subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial = reshape_vector(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)),masked_frame,indexedframe);

strongestFrame = -1*u_spatial/max(abs(u_spatial(:)))+0.5;
szstr = size(strongestFrame,1);
ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);
mmask = sum(model.totalConeMap3D,3)==0;
xx = repelem(mmask,1,1,3).*ssf; 
aa = model.totalConeMap3D;  aa(repelem(sum(aa,3),1,1,3)<0.5)=aa(repelem(sum(aa,3)<0.5,1,1,3))+0.5;aa =(aa-0.5)/max(aa(:)-0.5) + 0.5;
mxt = ((aa(repelem(sum(aa,3)~=1.5,1,1,3))));mx= max(mxt(:));
xx(repelem(sum(aa,3)~=1.5,1,1,3))= mxt;
xx = xx(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

% mag = repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf + model.totalConeMap3D*2;
% mag = mag(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

imagesc(xx);
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);

title(sprintf('SU # %d',ifilt));
%caxis([-max(mag(:)),max(mag(:))]);
%caxis([0,max(aa(:))])
end
suptitle(sprintf('Fraction of data %0.02d',frac_list(ifit)));
end


%% %% show SUs - in null space

P_null = P_mat.P_null;
P_sta = P_mat.P_sta;

for ifit=1:length(frac_list)
 fitGMLM = fitGMLM_log{ifit};   
sta_dim1 = size(mask2,1);
sta_dim2 = size(mask2,2);
indexedframe = reshape(1:sta_dim1*sta_dim2,[sta_dim1,sta_dim2]);
masked_frame = indexedframe(logical(mask2));

h=figure;
for ifilt=1:nSU
subplot(szz(nSU,1),szz(nSU,2),ifilt)
u_spatial = reshape_vector(gather(P_null*(fitGMLM.Linear.filter{ifilt}(1:length(masked_frame)))),masked_frame,indexedframe);

strongestFrame = -1*u_spatial/max(abs(u_spatial(:)))+0.5;
szstr = size(strongestFrame,1);
ssf = repelem(strongestFrame,model.gridSzX/szstr,model.gridSzX/szstr,3);
mmask = sum(model.totalConeMap3D,3)==0;
xx = repelem(mmask,1,1,3).*ssf; 
aa = model.totalConeMap3D;  aa(repelem(sum(aa,3),1,1,3)<0.5)=aa(repelem(sum(aa,3)<0.5,1,1,3))+0.5;aa =(aa-0.5)/max(aa(:)-0.5) + 0.5;
mxt = ((aa(repelem(sum(aa,3)~=1.5,1,1,3))));mx= max(mxt(:));
xx(repelem(sum(aa,3)~=1.5,1,1,3))= mxt;
xx = xx(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

% mag = repelem(sum(model.totalConeMap3D,3)==0,1,1,3).*ssf + model.totalConeMap3D*2;
% mag = mag(1:ceil(max(model.conesX))+40,1:ceil(max(model.conesY))+40,:);

imagesc(xx);
axis image
set(gca,'xTick',[]);
set(gca,'yTick',[]);

title(sprintf('SU # %d',ifilt));
%caxis([-max(mag(:)),max(mag(:))]);
%caxis([0,max(aa(:))])
end
suptitle(sprintf('Null projected, Fraction of data %0.02d',frac_list(ifit)));
end

%% make convergence plot .. 

P_null = P_mat.P_null;
P_sta = P_mat.P_sta;

filters_end = fitGMLM_log{length(frac_list)}.Linear.filter;

sqe_log=[]; sqe_null_log=[];sqe_sta_log=[];proj_sta=zeros(length(frac_list),1);proj_null=zeros(length(frac_list),1);

for ifit =1:length(frac_list)
    fitGMLM = fitGMLM_log{ifit};   

    sqe = 0;sqe_null=0;sqe_sta=0;
    for ifilt=1:nSU
        sqe = sqe + norm(fitGMLM.Linear.filter{ifilt} - filters_end{ifilt})^2;
        sqe_null = sqe_null + (norm(P_null*(fitGMLM.Linear.filter{ifilt} - filters_end{ifilt}))^2) / (norm(filters_end{ifilt}))^2;
        sqe_sta = sqe_sta +  (norm(P_sta*(fitGMLM.Linear.filter{ifilt} - filters_end{ifilt}))^2) / (norm(filters_end{ifilt}))^2; 
        proj_sta(ifit) = proj_sta(ifit) +   (norm(P_sta*(fitGMLM.Linear.filter{ifilt}))^2)/norm(fitGMLM.Linear.filter{ifilt})^2;
        proj_null(ifit) = proj_null(ifit) + (norm(P_null*(fitGMLM.Linear.filter{ifilt}))^2)/norm(fitGMLM.Linear.filter{ifilt})^2;
    end

    sqe_log(ifit) = sqe;
    sqe_null_log(ifit)=sqe_null;
    sqe_sta_log(ifit)=sqe_sta;
end

figure;
plot(frac_list,(sqe_log));
hold on;
plot(frac_list,(sqe_null_log));
hold on;
plot(frac_list,(sqe_sta_log));
xlabel('Fraction of data used');
ylabel('squared error between filters');
legend('SQE','SQE in null direction','SQE in STA direction');

mse_data.frac_list = frac_list;
mse_data.sqe_log = sqe_log;
mse_data.sqe_null_log = sqe_null_log;
mse_data.sqe_sta_log = sqe_sta_log;
mse_data.proj_sta = proj_sta;
mse_data.proj_null = proj_null;
end

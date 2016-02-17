% compare python LNLN fits with matlab ASM fits

%% load data
off_pars = load('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_analysis/Off_parasol.mat')
cellID=3348;
mask = off_pars.totalMaskAccept_log(:,off_pars.cells==3348);

mm = off_pars.maskedMovdd(logical(mask),:);
%% fit exponential and quadratic ASM
fitLNLN_exp=cell(5,1);
for nSU=1:5
    nSU
[fitLNLN_exp{nSU},f_val] = fitGMLM_MEL_EM_bias(off_pars.Y(off_pars.cells==3348,:)',mm,sum(mask(:)),nSU,1)

K = zeros(sum(mask)+1,nSU);
for isu=1:nSU
K(1:sum(mask(:)),isu) = fitLNLN_exp{nSU}.Linear.filter{isu} ; 
K(sum(mask(:))+1,isu) = fitLNLN_exp{nSU}.Linear.bias{isu};
end
fitLNLN_exp{nSU}.Linear.K = K;

end

fitLNLN_quad=cell(5,1);
for nSU=1:5
    nSU
[fitLNLN_quad{nSU},f_val] = fitGMLM_MEL_EM_power2(off_pars.Y(off_pars.cells==3348,:)',mm,sum(mask(:)),nSU,1,2)

K = zeros(sum(mask),nSU);
for isu=1:nSU
K(:,isu) = fitLNLN_quad{nSU}.Linear.filter{isu} ;   
end
fitLNLN_quad{nSU}.Linear.K = K;

end

%% see sub-unit filters

% LNLN model
for nSU=3:5
fitLNLN = loadLNLN_matlab(sprintf('fit_softplus_Cid=[3348]_Nsub=%d_.mat',nSU));
W = fitLNLN.W;
h=plotSU(W,logical(mask))
suptitle('softplus LNLN')
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/filters_softplus_nSU_%d.fig',nSU));

% LNLN model
fitLNLN = loadLNLN_matlab(sprintf('fit_sigmoid_Cid=[3348]_Nsub=%d_.mat',nSU));
W = fitLNLN.W;
h=plotSU(W,logical(mask))
suptitle('sigmoid LNLN')
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/filters_sigmoid_nSU_%d.fig',nSU));


K = fitLNLN_quad{nSU}.Linear.K;
h=plotSU(K,logical(mask))
suptitle('Quadratic ASM')
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/filters_quadratic_nSU_%d.fig',nSU));


K = fitLNLN_exp{nSU}.Linear.K(1:end-1,:);
h=plotSU(K,logical(mask))
suptitle('Exponential ASM')
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/filters_exp_nSU_%d.fig',nSU));

end
%% See bipolar NL
figure;
% LNLN model
fitLNLN = loadLNLN_matlab('fit_softplus_Cid=[3348]_Nsub=1_.mat')
W = fitLNLN.W;
nl_fcn = @(x) fitLNLN.ganglionNL_fcn(fitLNLN.bipolarNL_fcn(x)+repelem(fitLNLN.sig,size(x,1),1))/120;
%nl_fcn = @(x) (fitLNLN.bipolarNL_fcn(x));
plot_bipolarNL_normalized(W,mm,nl_fcn,'kb',false,[]);

nSU=1;
K = fitLNLN_quad{nSU}.Linear.K;
nl_fcn = @(x) ((x.*(x>0)).^2 + fitLNLN_quad{nSU}.mu) /120;
plot_bipolarNL_normalized(K,mm,nl_fcn,'rb',false,[])

nSU=1;
K = fitLNLN_exp{nSU}.Linear.K;
nl_fcn = @(x) exp(x)/120;
plot_bipolarNL_normalized(K,[mm;ones(1,size(mm,2))],nl_fcn,'gb',false,[])

%% See lam function impose by the three models

nSU=5;

fitLNLN = loadLNLN_matlab('fit_softplus_Cid=[3348]_Nsub=5_.mat');

K = zeros(sum(mask),nSU);
for isu=1:nSU
K(:,isu) = fitLNLN_quad{nSU}.Linear.filter{isu} ;   
end

for idim=1:nSU
    idim
    for jdim=idim+1:nSU
        jdim
dim1=K(:,idim);dim2=K(:,jdim); filters = [dim1,dim2];
var_mm = zeros(size(mm,1),1);
for i=1:size(mm,1)
var_mm(i) = var(mm(i,:));
end
var_mm = mean(var_mm);
mean_mm = 0;

sz = size(mm,1);
sig11 = diag(var_mm*ones(sz,1));
sig22 = filters'*diag(var_mm*ones(sz,1))*filters;
sig12 = diag(var_mm*ones(sz,1))*filters;
sig21 = sig12';
post_var = sig11 - sig12*(sig22\sig21);
[E,V] = eig(post_var);post_var_sqrt = E*sqrt(abs(V));
n_samples =100000;

X_log=[];Y_log=[];lam_log=[];lam_quad_log=[];lam_exp_log=[];
ix=0;
for proj_dim1 = [-5:0.5:5]*sqrt(sig22(1,1))
    ix=ix+1;iy=0;
    proj_dim1
    for proj_dim2 = [-5:0.5:5]*sqrt(sig22(2,2))
        iy=iy+1;
        
        x = [proj_dim1;proj_dim2];
        post_mean = sig12*(sig22\x);
        
        samples = randn(sz,n_samples);
        samples_transformed = post_var_sqrt*samples + repelem(post_mean,1,n_samples);
        lam_LNLN = (fitLNLN.lam_fcn(samples_transformed'));
        X_log(ix,iy) = proj_dim1;Y_log(ix,iy)=proj_dim2;
        lam_log(ix,iy) =mean(lam_LNLN);
        
        % quadratic
        inpt_quad = samples_transformed'*fitLNLN_quad{nSU}.Linear.K;
        lam_quad = (sum((inpt_quad.*(inpt_quad>0)).^2,2)+fitLNLN_quad{nSU}.mu)/120;
        lam_quad_log(ix,iy) = mean(lam_quad);
        
        % exponential
        inpt_exp = [samples_transformed;ones(1,n_samples)]'*fitLNLN_exp{nSU}.Linear.K;
        lam_exp = sum(exp(inpt_exp),2)/120;
        lam_exp_log(ix,iy)=mean(lam_exp);
        
        
    end
end
h = figure;
contourf(X_log,Y_log,log(real(lam_log)));
colorbar;
hold on;
caxis([-6,2]);
plotcov2([0;0],sig22,'conf',0.9)
title('LNLN model softplus');
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/stim_rate_map_sigmoid_nSU_%d_idim_%d_jdim_%d.fig',nSU,idim,jdim));


h=figure;
contourf(X_log,Y_log,log(real(lam_quad_log)));
colorbar
hold on;
plotcov2([0;0],sig22)
title('LNLN model quad');
caxis([-6,2]);
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/stim_rate_map_quad_nSU_%d_idim_%d_jdim_%d.fig',nSU,idim,jdim));


h= figure;
contourf(X_log,Y_log,log(real(lam_exp_log)));
colorbar
hold on;
plotcov2([0;0],sig22)
title('LNLN model exp');
caxis([-6,2]);
savefig(h,sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/stim_rate_map_exp_nSU_%d_idim_%d_jdim_%d.fig',nSU,idim,jdim));
close all

    end
end


%% See Null response predictions

%%
% Condition strings
% Condition strings
nConditions=8;
condDuration=1200/4;
cond_str=cell(nConditions,1);
cond_str{1}='Original';
cond_str{2}='Null ';

interestingConditions=[1,2,3,4,5,6,7];

dataRuns_OFF_additivity = [3,4,6,7,9,11,13];
dataRuns_ON_additivity = [3,5,6,8,10,12,13];
movies_OFF_addivitiy =[1,2,5,6,10,14,13];
movies_ON_additivity = [1,4,5,8,12,16,13];
location = '/Volumes/Analysis/2015-10-29-2/d00_36-norefit';

 
movies = movies_OFF_addivitiy;
dataRuns = dataRuns_OFF_additivity;

% load movies
interval=4;
condMov=cell(nConditions,1);
rawMovFrames=1200/(4);
icnt=0;
% make pixel histogram
for imov=movies
    [stim,height,width,header_size] = get_raw_movie(sprintf('/Volumes/Data/2015-10-29-2/Visual/%d.rawMovie',imov),rawMovFrames,1);
    subtract_movies{3}=mean(stim,1);
    subtract_movies{3}=mean(stim,1)*0+127.5;
    movie=stim-repmat(subtract_movies{3},[rawMovFrames,1,1]);
    movie=movie/255;
    
    icnt=icnt+1;
    qq=permute(movie,[2,3,1]);
    ifcnt = 0;
    condMov{icnt}=zeros(size(qq,1),size(qq,2),size(qq,3)*interval);
    for iframe=1:size(qq,3)
        for irepeat=1:interval
            ifcnt=ifcnt+1;
            condMov{icnt}(:,:,ifcnt)=qq(:,:,iframe)+0.5; % cond mov is between 0 and 1 now!
        end
        
    end
    
end

condDuration=10;
nConditions=1;
    for idata=1:length(dataRuns)
    Null_datafile = sprintf('%s/data0%02d',location,dataRuns(idata));
    neuronPath = [Null_datafile,sprintf('/data0%02d.neurons',dataRuns(idata))];
    [spkColl,spkCondColl{idata},h]=plot_raster_script_pc2015_09_23_0_light(cellID,nConditions,condDuration,cond_str,neuronPath);
    end
    
for nSU = 1:5
for icond=3:4;
movd = (condMov{icond}-0.5)*2;  % *2 if it is a low contrast movie!
 % movie is twice big for LNLN python fits compared to ASM fits!! for LNLN
 % fits, multiply the movie by 2!
maskedMov_cond= filterMov(movd,mask,off_pars.ttf_avg);
       
% predict responses
nTrials=30;
% LNLN

fitLNLN = loadLNLN_matlab(sprintf('fit_softplus_Cid=[3348]_Nsub=%d_.mat',nSU));

predictLNLN = poissrnd(repelem(fitLNLN.lam_fcn(maskedMov_cond'),1,nTrials)');

% Quad 
predictQuad = predictGMLM_gamma2_lr(fitLNLN_quad{nSU},maskedMov_cond,nTrials,2,1)';

% exp 
predictExp = predictGMLM_bias_lr(fitLNLN_exp{nSU},maskedMov_cond,nTrials,1)';

[x_lnln,y_lnln] = plotSpikeRaster(boolean(predictLNLN),'PlotType','vertline')
[x_e,y_e] = plotSpikeRaster(boolean(predictExp),'PlotType','vertline')
[x_q,y_q] = plotSpikeRaster(boolean(predictQuad),'PlotType','vertline')

close all
h= figure();
plot(spkCondColl{icond}.xPoints/20000,spkCondColl{icond}.yPoints,'k');
hold on;
plot(x_lnln/120,y_lnln-30,'r');
hold on;
plot(x_q/120,y_q-60,'k');
hold on;
plot(x_e/120,y_e-90,'r');
title(sprintf('%d SU, %d cond',nSU,icond));
xx = sprintf('/Volumes/Lab/Users/bhaishahster/pc2015_10_29_2_analysis_fits/data002_LNLN_compare_movie_has higher contrast/predictions_softplus_quad_exp_nSU_%d_icond_%d.fig',nSU,icond)
(ylim([-90,30]));set(gca,'yTick',[]);
savefig(h,xx);

end
end
%% AKHeitman 2014-04-27
% This will be admittedly long and ugly.  But more self contained.
% GLMPars = GLMParams  or GLMPars = GLMParams(GLMType.specialchange_name)
% Sensitive to naming changes in GLMParams.
% Only saves stuff (and calls directories) if we are in a troubleshooting mode
% Heavily GLMType dependent computations will be carried out here
% Outsourcable computations will be made into their own functions
% troubleshoot optional
% need troublshoot.doit (true or false), 
% troubleshoot.plotdir,
% troubleshoot.name



% CALLS which use GLMType:
%  prep_paramindGP
%  prep_stimcelldependentGPXV


function glm_dataexport(GLMType,fitspikes,fitmovie,testspikes_raster,testmovie,inputstats,glm_cellinfo, BD)


%% NON STIMULUS STUFF, can leave alone, will be added in linearly 
% setup
fittedGLM.cell_savename = glm_cellinfo.cell_savename;
fittedGLM.d_save        = glm_cellinfo.d_save;
fittedGLM.cellinfo      = glm_cellinfo;
GLMPars           = GLMParams;
if isfield(GLMType, 'specialchange') && GLMType.specialchange
    GLMPars = GLMParams(GLMType.specialchange_name);
end
fittedGLM.GLMPars = GLMPars;
fittedGLM.GLMType = GLMType;
if isfield(GLMType, 'debug') && GLMType.debug
    GLMPars.optimization.tolfun = 1;
end

% Timing
frames = size(fitmovie,3);
bins   = frames * GLMPars.bins_per_frame;
t_bin  = glm_cellinfo.computedtstim / GLMPars.bins_per_frame; % USE THIS tstim!! %
fittedGLM.t_bin = t_bin;
fittedGLM.bins_per_frame = GLMPars.bins_per_frame;

% % Make bases for the other inputs
% bin_size      = t_bin;
% if GLMType.PostSpikeFilter
%     basis_params  = GLMPars.spikefilters.ps;
%     ps_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
% end
% if GLMType.CouplingFilters
%     basis_params  = GLMPars.spikefilters.cp;
%     cp_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
% end
% if isfield(GLMType,'Contrast') && GLMType.Contrast
%     basis_params  = GLMPars.spikefilters.C;
%     C_basis      = prep_spikefilterbasisGP(basis_params,bin_size);
% end
% if isfield(GLMType, 'Saccades')
%    basis_params = GLMPars.saccadefilter;
%    sa_basis = prep_spikefilterbasisGP(basis_params, bin_size);
% end
% clear bin_size basis_params

% Convolve with Bases to get covariates
home_sptimes = fitspikes.home';
home_spbins  = ceil(home_sptimes / t_bin);
home_spbins = home_spbins(home_spbins < bins);
train_spikes = zeros(bins,1);
train_spikes(home_spbins) = 1;


test_sptimes = testspikes_raster.home{10}; % just for now, I will pick a random raster to test
test_spbins  = ceil(test_sptimes / t_bin);
test_spbins = test_spbins(test_spbins < size(testmovie,3)*GLMPars.bins_per_frame);
test_spikes = zeros(size(testmovie,3)*GLMPars.bins_per_frame, 1);
test_spikes(test_spbins) = 1;
% if GLMType.PostSpikeFilter
%     basis         = ps_basis';
%     PS_bin        = prep_convolvespikes_basis(home_spbins,basis,bins);
% end
center_coord       = glm_cellinfo.slave_centercoord;
% if isfield(GLMType,'Contrast') && GLMType.Contrast
%     stimsize.width  = size(fitmovie,1);
%     stimsize.height = size(fitmovie,2);
%     ROIcoord        = ROI_coord(GLMPars.spikefilters.C.range, center_coord, stimsize);
%     contrast = imresize(squeeze(mean(mean(double(fitmovie(ROIcoord.xvals,ROIcoord.yvals, :))-0.5))), [bins 1],'nearest');
%     C_bin = zeros(GLMPars.spikefilters.C.filternumber,bins);
%     for i = 1:GLMPars.spikefilters.C.filternumber
%        tmp = conv(contrast, C_basis(:,1), 'full');
%        C_bin(i,:) = tmp(1:bins); 
%     end
% end
% % NBCoupling 05-28-14
% if GLMType.CouplingFilters;
%     basis = cp_basis';
%     for j_pair=1:GLMPars.spikefilters.cp.n_couplings
%         %spikes of neighbor neurons NB
%         neighbor_sptimes = neighborspikes.home{j_pair}';
%         neighbor_spbins  = ceil(neighbor_sptimes / t_bin);
%         neighbor_spbins = neighbor_spbins(find(neighbor_spbins < bins) );
%         CP_bin{j_pair}=prep_convolvespikes_basis(neighbor_spbins,basis,bins);
%     end
% end
% if isfield(GLMType, 'Saccades')
%    basis = sa_basis';
%    spbins = 1:120:bins;
%    SA_bin = prep_convolvespikes_basis(spbins, basis, bins);
% end
% if GLMType.TonicDrive
%     MU_bin = ones(1,bins);
% end
% 
% % Make Covariate Vector of NON stimulus stuff
% [paramind] =  prep_paramindGP(GLMType, GLMPars);
% p_init     = .01* ones(paramind.paramcount,1);
% rawfit.init = p_init;
% glm_covariate_vec = NaN(paramind.paramcount , bins );  % make sure it crasheds if not filled out properly
% if isfield(paramind, 'MU')
%     glm_covariate_vec( paramind.MU , : ) = MU_bin;
% end
% if isfield(paramind, 'PS')
%     glm_covariate_vec( paramind.PS , : ) = PS_bin;
% end
% % NBCoupling 05-28-14
% if isfield(paramind, 'CP')
%     for j_pair=1:GLMPars.spikefilters.cp.n_couplings
%         glm_covariate_vec( paramind.CP{j_pair} , : ) = CP_bin{j_pair};
%     end
% end
% if isfield(paramind, 'SA')
%     glm_covariate_vec(paramind.SA, :) = SA_bin;
% end
% if isfield(paramind, 'C')
%     glm_covariate_vec(paramind.C, :) = C_bin;
% end

GLMType.stimfilter_mode = 'rk1';
% STIM
% WN_STA             = double(glm_cellinfo.WN_STA);
[~,train_movie]    = prep_stimcelldependentGPXV(GLMType, GLMPars, fitmovie, inputstats, center_coord);%, WN_STA);
train_movie = train_movie';
[~,test_movie]    = prep_stimcelldependentGPXV(GLMType, GLMPars, testmovie, inputstats, center_coord);%, WN_STA);
test_movie = test_movie';

 glm_cellinfo.slave_centercoord = center_coord;
 center(1) = center_coord.x_coord;
 center(2) = center_coord.y_coord;

 WN_STA = glm_cellinfo.WN_STA;
 stimsize.width  = size(fitmovie,1);
 stimsize.height = size(fitmovie,2);
 ROIcoord        = ROI_coord(GLMPars.stimfilter.ROI_length, center_coord, stimsize);
 clear stimsize
 STA = WN_STA(ROIcoord.xvals,ROIcoord.yvals,:);
 klen = size(STA,1);
 duration = size(STA, 3);
 STA = reshape(STA, [klen^2,duration])  - mean(STA(:));
 [U,S,V]  = svd (STA);
 imagesc(reshape(U(:,1), [klen, klen]))
 title('Initial Space Filter')
 axis image
 clear X_frame_temp
 plot(V(:,1))
 tempFilt = V(:,1);
 clear STA U S V
 
fit_batch = length(test_movie)/2;
test_batch = length(test_movie)/2;


disp('Saving files')



save([BD.GLM_develop_output_raw '/DataExport/' glm_cellinfo.cell_savename '.mat'], 'test_spikes', 'train_spikes', 'test_movie', 'train_movie','center','tempFilt','fit_batch', 'test_batch', '-v7.3')


end

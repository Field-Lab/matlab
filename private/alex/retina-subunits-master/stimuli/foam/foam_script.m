function foam_script(dur,sd,tau,contrast,threshold,plotting)

% foam_script(dur,sd,tau,contrast,threshold,plotting)
%
% generates a Gaussian-filtered stimulus and writes it out to
% a binary movie file. options for varied spatial and temporal
% smoothing, along with variable contrast, or thresholding (binary
% and ternary)
%
% inputs:
% dur -- duration in seconds (assuming a refresh interval of 8 frames
%        on a 120 Hz monitor)
% sd -- standard deviation of the gaussian (for spatial smoothing)
% tau -- time constant of exponential (for temporal smoothing)
% contrast -- rms contrast, set to 1 for thresholded (binary or ternary)
% threshold -- type of threhsold, 'b' - binary, 't' - ternary,
%              only used if contrast = 1
% plotting -- show stimuli during generation?


%% settings for debugging
 dur = 3;
 sd = 3;
 tau = 2;
 plotting = 1;
 contrast = 1;
 threshold = 'b';
%%
p.seed = 100;
p.refresh = 120; % hard coded
p.interval = 8; % hard coded
p.dur_sec = dur;
p.dims = [300 300 (p.refresh/p.interval)*p.dur_sec];
p.gauss_spatial_sd = sd;
p.tau_temp_interval = tau;
p.contrast_sd = contrast;
p.thresholding = threshold;
testing = 0;
%%
fprintf('(foam_script) making stimuli\n');
randn('seed',p.seed);

if p.gauss_spatial_sd
  filt = fspecial('gaussian',roundOdd(8*p.gauss_spatial_sd),p.gauss_spatial_sd);
  filt = filt/norm(filt(:));
end

if p.tau_temp_interval
  dur = p.dims(3)+10;
else
  dur = p.dims(3);
end
stim_vec = randn(p.dims(1),p.dims(2),dur);

fprintf('(foam_make_frames) making frames\n');
T=text_waitbar;

for i=1:dur
  T=text_waitbar(T,i/dur - 0.01);

  if testing
    tmp1 = 1;
    tmp2 = 1;
    tmp3 = 0;
    tmp4 = 0;
    im1 = [tmp1 tmp3; tmp4 tmp2];
    if mod(i,2) == 1
      im1 = double(~im1);
    end
    im1 = imresize(im1,p.dims(1:2),'nearest');
    
  else
    if p.gauss_spatial_sd
      stim_vec(:,:,i) = cconv2(stim_vec(:,:,i),filt);
    end
    if p.contrast_sd == 1 && p.tau_temp_interval
      stim_vec(:,:,i) = stim_vec(:,:,i)*0.1 + 0.5;
    elseif p.contrast_sd == 1 && p.tau_temp_interval == 0
      if strcmp(p.thresholding,'b')
        tmp = zeros(size(stim_vec(:,:,i)));
        tmp(stim_vec(:,:,i)<=0) = 0;
        tmp(stim_vec(:,:,i)>0) = 1;
        stim_vec(:,:,i) = tmp;
      end
      if strcmp(p.thresholding,'t')
        tmp = zeros(size(stim_vec(:,:,i)));
        tmp(stim_vec(:,:,i)<norminv(1/3)) = 0;
        tmp(stim_vec(:,:,i)>=norminv(1/3) & stim_vec(:,:,i)<norminv(2/3)) = 0.5;
        tmp(stim_vec(:,:,i)>=norminv(2/3)) = 1;
        stim_vec(:,:,i) = tmp;
      end
    else
      stim_vec(:,:,i) = stim_vec(:,:,i)*p.contrast_sd + 0.5;
      clip_frac(i) = sum(vector(stim_vec(:,:,i)<0 | stim_vec(:,:,i)>1))/length(vector(stim_vec(:,:,i)));
      stim_vec(:,:,i) = max(stim_vec(:,:,i),0);
      stim_vec(:,:,i) = min(stim_vec(:,:,i),1);
    end
    
  end
  
  %stim_vec(:,i) = im1(:);
  
end
fprintf('\n');
%%
%clip_count = sum(im1(:)<0|im1(:)>1);
%if clip_count > 10
%    warning('clipping %g pixels\n',clip_count);
%end
%im1 = clip(im1,0,1);

%%

if p.tau_temp_interval
  
  T = (p.refresh/p.interval);
  tau = (1000/T)*p.tau_temp_interval;
  fprintf('(foam_make_frames) doing temporal filtering\n');
  T_=text_waitbar;
  for i=1:size(stim_vec,1);
    T_=text_waitbar(T_,i/size(stim_vec,1) - 0.01);
    for j=1:size(stim_vec,2)
      stim_vec(i,j,:) = expsmooth(squeeze(stim_vec(i,j,:)),T,tau);
    end
  end
  fprintf('\n');
  %%
  fprintf('(foam_make_frames) getting new contrast\n');
  
  if p.contrast_sd == 1
    %%
    if strcmp(p.thresholding,'b');
      bnd1 = prctile(vector(stim_vec(:,:,11:min(size(stim_vec,3),1000))),50);
    end
    if strcmp(p.thresholding,'t');
      bnd1 = prctile(vector(stim_vec(:,:,11:min(size(stim_vec,3),1000))),(1/3)*100);
      bnd2 = prctile(vector(stim_vec(:,:,11:min(size(stim_vec,3),1000))),(2/3)*100);
    end
    fprintf('(foam_make_frames) fixing contrast and doing clipping\n');
    T_ = text_waitbar;
    for i=1:dur
      T_=text_waitbar(T_,i/dur - 0.01);
      tmp = stim_vec(:,:,i);
      if strcmp(p.thresholding,'b')
        tmp(stim_vec(:,:,i)<=bnd1) = 0;
        tmp(stim_vec(:,:,i)>bnd1) = 1;
      end
      if strcmp(p.thresholding,'t');
        tmp(stim_vec(:,:,i)<bnd1) = 0;
        tmp(stim_vec(:,:,i)>=bnd1 & stim_vec(:,:,i)<bnd2) = 0.5;
        tmp(stim_vec(:,:,i)>=bnd2) = 1;
      end
      stim_vec(:,:,i) = tmp;
    end
    fprintf('\n');
  else
    con_factor = p.contrast_sd/mystd2((stim_vec(:,:,11:min(size(stim_vec,3),200))));
    %%
    fprintf('(foam_make_frames) fixing contrast and doing clipping\n');
    T_ = text_waitbar;
    for i=1:dur
      T_=text_waitbar(T_,i/dur - 0.01);
      stim_vec(:,:,i) = (stim_vec(:,:,i) - 0.5)*con_factor + 0.5;
      stim_vec(:,:,i) = max(stim_vec(:,:,i),0);
      stim_vec(:,:,i) = min(stim_vec(:,:,i),1);
    end
    fprintf('\n');
  end
end


%%
vidObj = VideoWriter('test.avi');
vidObj.FrameRate = 15;
open(vidObj);
clf
if plotting;
  %set(gcf,'Position',[873 386 560 420]);
  colormap(gray)
  for i=11:50
    imshow((stim_vec(:,:,i)),'DisplayRange',[0 1]);
    axis image off;
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
end
close(vidObj);
%%
fprintf('(foam_script) saving file \n');
if testing
  saveName = sprintf('foam_testing_checkerboard_simple.rawMovie');
else
  if p.contrast_sd == 1
    saveName = sprintf('crest_s%g_t%g_%s_d%g.rawMovie',p.gauss_spatial_sd,p.tau_temp_interval,p.thresholding,p.dur_sec);
  else
    saveName = sprintf('foam_s%g_t%g_c%02g_d%g.rawMovie',p.gauss_spatial_sd,p.tau_temp_interval,p.contrast_sd*100,p.dur_sec);
  end
end
if p.contrast_sd == 1
  savePath = sprintf('/share/users-archive/freeman/crest_%s/',p.thresholding);
else
  savePath = sprintf('/share/users-archive/freeman/foam_c%02g/',p.contrast_sd*100);
end
%%
if p.tau_temp_interval
  foam_export_movie_obvius(fullfile(savePath,saveName),stim_vec(:,:,11:end),p);
else
  foam_export_movie_obvius(fullfile(savePath,saveName),stim_vec,p);
end


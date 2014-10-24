function foam_mk_movie(

p.seed = 100;
p.refresh = 120;
p.target_interval = 4;
p.dur_sec = 0.1;
p.dims = [300 300 (p.refresh/p.target_interval)*60*p.dur_sec];
p.gauss_spatial_sd = 1;
p.tau_temp_interval = 0;
p.contrast_sd = 0.1;
testing = 0;
%%
randn('seed',p.seed);
if p.gauss_spatial_sd
    filt = fspecial('gaussian',roundOdd(8*p.gauss_spatial_sd),p.gauss_spatial_sd);
    filt = filt/norm(filt(:));
end
stim_vec = zeros(prod(p.dims(1:2)),p.dims(3));
%imagesc(filt);
%%

if p.tau_temp_interval
    dur = p.dims(3)+10;
else
    dur = p.dims(3);
end
T=text_waitbar;

for i=1:dur
    T=text_waitbar(T,i/dur - 0.01);
    im1 = randn(p.dims(1),p.dims(2));
    
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
            im1 = cconv2(im1,filt);
        end
        im1 = im1*p.contrast_sd + 0.5;
        clip_count = sum(im1(:)<0|im1(:)>1);
        if clip_count > 10
            warning('clipping %g pixels\n',clip_count);
        end
        im1 = clip(im1,0,1);
    end
    
    stim_vec(:,i) = im1(:);
    
end
%%
if p.tau_temp_interval
    stim_vec_smooth = zeros(size(stim_vec));
    T = (p.refresh/p.target_interval);
    tau = (1000/T)*p.tau_temp_interval;
    for i=1:size(stim_vec,1);
        tmp = expsmooth(stim_vec(i,:)',T,tau);
        stim_vec_smooth(i,:) = tmp';
    end
    %%
    stim_vec_smooth = (stim_vec_smooth - 0.5) * ...
        (p.contrast_sd/std(stim_vec_smooth(:))) + 0.5;
    %%
    clip_count = sum(stim_vec_smooth(:)<0|stim_vec_smooth(:)>1);
    if clip_count
        warning('clipping %g pixels\n',clip_count);
    end
    stim_vec_smooth = stim_vec_smooth(:,11:end);
    stim_vec_smooth = clip(stim_vec_smooth,0,1);
    
    stim = reshape(stim_vec_smooth,p.dims);
else
    stim = reshape(stim_vec,p.dims);
end
%%


%%
clear F
colormap(gray)
for i=1:10%p.dims(3);
    imagesc((stim(:,:,i)),[0 1]);
    axis image off;
    pause(1/(p.refresh/p.target_interval));
end

%%
if testing
    saveName = sprintf('saved/foam_testing_checkerboard_simple.rawMovie');
else
    saveName = sprintf('saved/foam_testing_sd_%g_tau_%g.rawMovie',p.gauss_spatial_sd,p.tau_temp_interval);
end
foam_export_movie_obvius(saveName,stim,p);



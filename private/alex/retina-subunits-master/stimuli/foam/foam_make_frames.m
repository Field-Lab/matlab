function stim_vec = foam_make_frames(p,testing)

if p.gauss_spatial_sd
    filt = fspecial('gaussian',roundOdd(8*p.gauss_spatial_sd),p.gauss_spatial_sd);
    filt = filt/norm(filt(:));
end
stim_vec = zeros(prod(p.dims(1:2)),p.dims(3));

if p.tau_temp_interval
    dur = p.dims(3)+10;
else
    dur = p.dims(3);
end
fprintf('(foam_make_frames) making frames\n');
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
        %if clip_count > 10
        %    warning('clipping %g pixels\n',clip_count);
        %end
        im1 = clip(im1,0,1);
    end
    
    stim_vec(:,i) = im1(:);
    
end
fprintf('\n');
%%

if p.tau_temp_interval
    
    T = (p.refresh/p.target_interval);
    tau = (1000/T)*p.tau_temp_interval;
    fprintf('(foam_make_frames) doing temporal filtering\n');
    T_=text_waitbar;
    for i=1:size(stim_vec,1);
        T_=text_waitbar(T_,i/size(stim_vec,1) - 0.01);
        tmp = expsmooth(stim_vec(i,:)',T,tau);
        stim_vec(i,:) = tmp';
    end
    fprintf('\n');
    %%
    stim_vec = stim_vec(:,11:end);
    stim_vec = (stim_vec - 0.5) * ...
        (p.contrast_sd/std(stim_vec(:))) + 0.5;
    %%
    clip_count = sum(stim_vec(:)<0|stim_vec(:)>1);
    if clip_count
        warning('clipping %g pixels\n',clip_count);
    end
    stim_vec = clip(stim_vec,0,1);
    
    stim_vec = reshape(stim_vec,p.dims);
else
    stim_vec = reshape(stim_vec,p.dims);
end
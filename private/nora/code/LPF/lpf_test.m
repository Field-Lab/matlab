
frames = 100;
mine = lpf_make_mat('frames', frames);

%% Compare mean and SD
clear a b m b2 m2
for i = 1:frames
   a = mine(:,:,i);
   b(i) = std(a(:));
   m(i) = mean(a(:));
end

for i = 1:frames
   a = double(X(i,:,:));
   b2(i) = std(a(:));
   m2(i) = mean(a(:));
end


plot(m)
hold on
plot(m2)
legend('mine', 'orig')
title('mean')

figure;
plot(b)
hold on
plot(b2)
legend('mine', 'orig')
title('std')

%% Distributions
a = hist(double(X(:)), 0:255);
a2 = hist(double(mine(:)), 0:255);
plot(a2/frames)
hold on
plot(a/frames)
legend('mine', 'orig')
title('Distribution of Pixel Values')
 
%% Compare correlations

% Time
tlag = 100;
t = zeros(2*tlag+1, 1);
t2 = zeros(2*tlag+1, 1);
for i = 1:320
    for j = 1:160
        t = t + xcorr(squeeze(mine(i,j,:)), squeeze(mine(i,j,:)), tlag, 'unbiased');
        t2 = t2 + xcorr(squeeze(X(:,i,j)), squeeze(X(:,i,j)), tlag, 'unbiased');
    end
end
plot(t/max(t))
hold on
plot(t2/max(t2))
title('Temporal Correlation')
legend('mine', 'orig')

% Space
s = zeros(2*tlag+1, 1);
s2 = zeros(2*tlag+1, 1);
for i = 1:160
    for j = 1:frames
        s = s + xcorr(squeeze(mine(:,i,j)), squeeze(mine(:,i,j)), tlag, 'unbiased');
        s2 = s2 + xcorr(squeeze(X(j,:,i)), squeeze(X(j,:,i)), tlag, 'unbiased')';
    end
end
figure;
plot(s/max(s))
hold on
plot(s2/max(s2))
title('Spatial Correlation')
legend('mine', 'orig')



%% Watch movies
idx = 1:10;
for i = 1:frames
   subplot(2,1,1)
   imagesc(mine(:,:,i)')
   colormap gray
   caxis([0 255])
   axis image
   subplot(2,1,2)
   imagesc(squeeze(X(i,:,:))')
   colormap gray
   caxis([0 255])
   axis image
   pause(0.1)
end

%% are the movies off by some offset?

me = squeeze(mine(1,1,:));
orig = double(squeeze(X(:,1,1)));
xcor = zeros(901, 1);
for i = 0:(frames-100)
    idx = i+(1:100);
    xcor(i+1) = corr(orig, me(idx));
end
plot(xcor)

%%
frame_10 = mine(1,:,5);
plot(frame_10, 'k')
hold on
for i = 6
    fr = squeeze(X(i,1,:));
    a = plot(fr);
    pause()
    delete(a)
end

%%
line = 1;
movie_alex = gauss_Test('frames', 2);
scalar_previous = exp(-1/6);
scalar_current = sqrt(1 - scalar_previous^2);
movie_alex = movie_alex - 127.5;
alex = scalar_previous*movie_alex(:,line,1)+scalar_current*movie_alex(:,line,2);
alex = alex + 127.5;
% alex = conv(alex, gausswin(10), 'same');
plot(alex)
hold on
plot(X(1,:,line));
plot(mine(:,line))
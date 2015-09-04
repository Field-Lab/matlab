datarun = load_data(['/Volumes/Analysis/2010-09-24-1/data034/data034']);
datarun = load_params(datarun,'verbose',1);
datarun = load_neurons(datarun);

[inputs_all, refresh, duration] = get_wn_movie_ath_rgb(datarun, 'RGB-gaussian-8-2-0.48-11111.xml');
inputs = inputs_all(:);

writerObj = VideoWriter('/Users/alexth/Desktop/video.avi');
writerObj.FrameRate = 15;
writerObj.Quality=100;
open(writerObj);

figure
set(gcf, 'position', [-796   643   427   427] )
subplot(1,1,1)
set(gca, 'position', [0 0 1 1])
colormap('gray')
for j=1:60
    a = inputs(1:80*80);
    a = reshape(a, 80, 80);
    inputs(1:80*80) = [];

    tmp = repmat(a, 1,1,3);

    imagesc(tmp)

    imagesc(a)
    
    set(gca, 'DataAspectRatio', [1 1 1])
    axis off
    movieFrames(j)=getframe;
    writeVideo(writerObj,movieFrames(j));
end
close(writerObj);






path2data = '/Volumes/Analysis/2015-03-09-2/d05-27-norefit/';
for_use = ['data005'; 'data008'; 'data014'; 'data018'; 'data022'; 'data026'];
datarun = cell(6,1);
for i=6
    datarun{i} = load_data(fullfile(path2data, for_use(i,:), for_use(i,:)));
    datarun{i} = load_params(datarun{i},'verbose',1);
    datarun{i} = load_sta(datarun{i},'load_sta',[]);
    datarun{i} = get_sta_fits_from_vision(datarun{i});
    datarun{i} = load_sta(datarun{i},'load_sta', 'all');
end

sta = squeeze(datarun{6}.stas.stas{18});
tc = datarun{6}.vision.timecourses(18).g(11:end);
sta = sta(:,:,11:end);
interp_sta = zeros(80,80, 40);
for i=1:80
    for j=1:80
        interp_sta(i,j,:) = interp1(1:2:2*20, squeeze(sta(i,j,:))',1:2*20);
    end
end
interp_sta = interp_sta(:,:,3:38);
interp_sta = interp_sta/(max(abs(interp_sta(:)))*2);
interp_sta = interp_sta + 0.5;

my_sta = sta.*reshape(repmat(tc', size(sta,1)^2,1), size(sta,1), size(sta,2),20 );
my_sta = sum(my_sta,3);
my_sta = my_sta/(max(my_sta(:))*2)+0.5;
my_sta = 1-repmat(my_sta, 1,1,3);
    
figure
set(gcf, 'position', [-796   643   427   427] )
subplot(1,1,1)
set(gca, 'position', [0 0 1 1])
writerObj = VideoWriter('/Users/alexth/Desktop/video_sta.avi');
writerObj.FrameRate = 10;
writerObj.Quality=100;
open(writerObj);
% 
% imagesc(my_sta)
% set(gca, 'DataAspectRatio', [1 1 1])
% axis off
% movieFrames(1)=getframe;
% writeVideo(writerObj,movieFrames(1));
for j=1:35
    tmp = interp_sta(:,:,j);
    tmp = repmat(tmp, 1,1,3);

    imagesc(tmp)
    
    set(gca, 'DataAspectRatio', [1 1 1])
    axis off
    movieFrames(j)=getframe;
    writeVideo(writerObj,movieFrames(j));
end
imagesc(my_sta)
set(gca, 'DataAspectRatio', [1 1 1])
axis off
movieFrames(37)=getframe;
writeVideo(writerObj,movieFrames(37));
close(writerObj);

mvpath='/Volumes/Data/Stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/';
my_movie=zeros(320,320,3600);
cnt=1;
for i=1:30
    load([mvpath, 'movie_chunk_', int2str(i), '.mat']);
    movie=movie(81:240,:,:);
    movie=imresize(movie,2,'nearest');
    my_movie(:,:,cnt:cnt+119)=movie;
    cnt=cnt+120;
end
my_movie=my_movie/255-0.5;
clear mvpath movie




data='008';

for i=1
    switch data
        case '004'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            ndf=4;
            nmdata='009';
        case '008'
            wn_movie_name = 'BW-10-8-0.48-11111-32x32.xml';
            ndf=4;
            nmdata='009';
        case '011'
            wn_movie_name = 'BW-16-8-0.48-11111-20x20.xml';
            ndf=3;
            nmdata='012';
        case '014'
            wn_movie_name = 'BW-8-8-0.48-11111-40x40.xml';
            ndf=3;
            nmdata='012';
        case '015'
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            ndf=2;
            nmdata='016';
        case '018'  
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            ndf=2;
            nmdata='016';
        case '019' 
            wn_movie_name = 'BW-16-6-0.48-11111-20x20.xml';
            ndf=1;
            nmdata='020';
        case '022' 
            wn_movie_name = 'BW-8-6-0.48-11111-40x40.xml';
            ndf=1;
            nmdata='020';
        case '023' 
            wn_movie_name = 'BW-16-4-0.48-11111-20x20.xml';
            ndf=0;
            nmdata='024';
        case '027' 
            wn_movie_name = 'BW-4-2-0.48-11111-80x80.xml';
            ndf=0;
            nmdata='024';
    end  
    starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
    starun = load_params(starun,'verbose',1);
    starun = set_polarities(starun);
    starun = load_neurons(starun);
    starun = load_sta(starun,'load_sta','all','keep_java_sta',true);

    nmrun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',nmdata,'-from-d05-d27/data',nmdata,'-from-d05-d27']);
    nmrun = load_params(nmrun,'verbose',1);
    nmrun = load_neurons(nmrun);
end
clear data nmdata

addpath(genpath('/Users/alexth/test4/matlab/private/nora/'))

data='008';
starun = load_data(['/Volumes/Analysis/2015-03-09-2/d05-27-norefit/data',data,'-from-d05-d27/data',data,'-from-d05-d27']);
starun = load_params(starun,'verbose',1);
starun = set_polarities(starun);
starun = load_neurons(starun);
starun = load_sta(starun,'load_sta','all','keep_java_sta',true);
wn_movie_name = 'BW-10-8-0.48-11111-32x32';
cells = 800;
glm = glm_fit_from_WN(cells, starun, wn_movie_name);
plotfilters(glm)


cells = [800];
datarunID = find(starun.cell_ids == cells);


[inputs, refresh, duration] = get_wn_movie_ath(starun, wn_movie_name); 
stim = zeros(32*32,size(inputs,2)*8);
for i=1:8
    stim(:,i:8:end) = inputs;
end
stim = reshape(stim, 32,32,size(stim,2));

spikes = starun.spikes{datarunID};

center = starun.stas.fits{datarunID}.mean;


STA_Test(spikes, stim, center);




if wn_movie_name(end)=='l'
    wn_movie_name = wn_movie_name(1:end-4);
end
glm = glm_fit_from_WN(cells, starun, wn_movie_name);
plotfilters(glm)


stix_size=10;
stim=zeros(320/stix_size,320/stix_size, 3600);
cnt_i=1;
for i=1:stix_size:320
    cnt_j=1;
    for j=1:stix_size:320
        tmp=reshape(my_movie(i:i+stix_size-1,j:j+stix_size-1, :), stix_size*stix_size, 3600);
        stim(cnt_i, cnt_j,:)=mean(tmp);
        cnt_j=cnt_j+1;
    end
    cnt_i=cnt_i+1;
end

stim1 = imresize(my_movie,1/stix_size, 'method','box');
a = stim1 - stim;
mean(a(:))

stim1 = permute(stim1,[2,1,3]);

glm_prediction = glm_predict(glm, stim1);
nnz(glm_prediction.rasters.glm_sim)


xval = struct;
xval.rasters = '';
plotraster(xval, glm)


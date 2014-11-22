function mov = generate_movie_ts(mov_params)

if(mov_params.type=='bw')

movie_spec=mov_params.movie_spec;%'/Volumes/Analysis/movie-xml/RGB-8-1-0.48-11111.xml' 
mdf_file=movie_spec;%'/Volumes/Analysis/deprecated/movie-xml2/RGB-8-1-0.48-11111.xml';
mov_time = mov_params.movie_len;
refresh=mov_params.refresh; % in ms
triggers=[0:100*refresh/1000:mov_time];

frames =ceil(mov_time*1000/refresh); % 10 minutes
[mov,height,width,duration,refresh] = get_movie(mdf_file, triggers,frames);
mov=(mov-0.5);
mov2=zeros(size(mov,2),size(mov,1),size(mov,3),size(mov,4));

for itime=1:frames
for icol=1:3
    mov2(:,:,icol,itime) = mov(:,:,icol,itime)';
end
end

mov=mov2;
end





end
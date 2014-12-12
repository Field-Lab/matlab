movie_time=120*10;
movie_log=zeros(320,160,1,movie_time);
icnt=0;
for imov=1:10
    imov
    icnt=icnt+1;
movv=load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_%d',imov));
movv=movv.movie;
movie_log(:,:,1,(icnt-1)*120+1:icnt*120)=movv;
end

movie_log=uint8(movie_log);
writerObj = VideoWriter('~/Dropbox/Lab/Transfer/naturual_movie.avi');
writerObj.FrameRate=120;
open(writerObj);
writeVideo(writerObj,movie_log);
close(writerObj);

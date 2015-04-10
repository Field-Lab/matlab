function [mov_buffered,movie_params]=generate_movie(movie_params)

display('Generating movies');

if(strcmp(movie_params.mov_type,'nsem'))

no_mov=movie_params.no_images_per_movie;%10;
start_image=movie_params.start_image;

movie_time=120*no_mov;%size(movie,3);


movie_log=zeros(320,160,movie_time);
icnt=0;
for imov=start_image:no_mov+start_image-1
    imov
    icnt=icnt+1;
movv=load(sprintf('/Volumes/Data/stimuli/movies/eye-movement/current_movies/NSbrownian_6000/matfiles/movie_chunk_%d',imov));
movv=movv.movie;
movie_log(:,:,(icnt-1)*120+1:icnt*120)=movv;
end
movie2=movie_log;
clear movie_log
var64=movie_params.var64;
var32=movie_params.var32;
mov=zeros(var64,var32,movie_time);
for itime=1:movie_time
mov(:,:,itime)=(imresize(movie2(:,:,itime),[var64,var32],'bilinear','Antialiasing',true)-movie_params.mean); % Doubt!!
%mov(:,:,itime)=movie2(var64,var32,itime);
display('Resizing a snippet');
end

% % Show movie
% figure
% for itime=1:movie_time
%     subplot(2,1,1);
%     imagesc(movie2(:,:,itime)');
%     caxis([0,255]);
%     axis image
%     %title(sprintf('Time: %d',itime));
%     colorbar
%     
%     subplot(2,1,2);
%    imagesc(mov(:,:,itime))
%    caxis([-128,128]);
%    colormap gray
%    axis image
%    %title(sprintf('Time: %d',itime));
%    colorbar
%    
%    pause(1/120);
% end


%mov=mov*127.5/max(abs(mov(:))); %contrast correction

% Add 
% Add one second of dark screen in front and back of movie .. 
% This will take care of initial and final conditions and make proper
% buffers
if(isfield(movie_params,'interval'))
blankFrames = double(uint8(0/movie_params.interval));

mov_buffered=zeros(var64,var32,movie_time+2*blankFrames)-0.25*(255);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
else
    
blankFrames = 0;
mov_buffered=zeros(var64,var32,movie_time+2*blankFrames)-0.25*(255);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
end
end

% BW noise
if(strcmp(movie_params.mov_type,'bw'))
movie_time=movie_params.movie_time;
var64=movie_params.var64;
var32=movie_params.var32;
mov=double(rand(var64,var32,movie_time)>0.5);
mov=mov-0.5;
mov=mov*movie_params.deviation/max(abs(mov(:)));
% White noise, but it's +1 or -1


% Add 
% Add one second of dark screen in front and back of movie .. 
% This will take care of initial and final conditions and make proper
% buffers
if(isfield(movie_params,'interval'))
blankFrames = double(uint8(0/movie_params.interval));

mov_buffered=zeros(var64,var32,movie_time+2*blankFrames);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
else
    
blankFrames = 0;
mov_buffered=zeros(var64,var32,movie_time+2*blankFrames);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
end
end

% BW noise used in before
if(strcmp(movie_params.mov_type,'bw-precomputed'))
movie_time=movie_params.movie_time;
var64=movie_params.var64;
var32=movie_params.var32;

if(~isfield(movie_params,'mdf_file'))
mdf_file='/Volumes/Analysis/stimuli/white-noise-xml/BW-10-1-0.48-11111.xml';
else
    mdf_file=movie_params.mdf_file;
end
mov=zeros(var64,var32,movie_time);
triggers=[0:100/120:movie_time*10/120]; % Or maybe, movie_time/120 ?? Doubt!
[mvi] = load_movie(mdf_file,triggers);
[~,height,width,duration,refresh] = get_movie_ath(mdf_file,...
    triggers, 1,2);

    for itime=1:movie_time
        F=mvi.getFrame(itime-1).getBuffer;
        mov(:,:,itime)=double(reshape(F(1:3:end),width,height)); % Flip in each direction ? Doubt!

    end
mov=mov-0.5;
mov=mov*movie_params.deviation/max(abs(mov(:)));


% Add 
% Add one second of dark screen in front and back of movie .. 
% This will take care of initial and final conditions and make proper
% buffers
if(isfield(movie_params,'interval'))
blankFrames = double(uint8(0/movie_params.interval));

mov_buffered=zeros(var64,var32,movie_time+2*blankFrames);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
else
    
blankFrames = 0;
mov_buffered=zeros(var64,var32,movie_time+2*blankFrames);
mov_buffered(:,:,blankFrames+1:end-blankFrames)=mov;
movie_time=movie_time+2*blankFrames;
movie_params.movie_time=movie_time;
end
end



end
function [stim,width,height,frames_generated]=load_movie_rawMovie(movie_file, frames)

%movie_file='/jacob/snle/lab/Development/Stimulus/movies/eye-movement/eye-120-10-s-5_5-30.rawMovie';


disp('Getting movie');tic

%get header
    tab=hex2dec('09');
    new1=hex2dec('0d');
    new2=hex2dec('0a');

    fclose all;
    fid=fopen(movie_file,'r');
    file=fscanf(fid, '%c', 2000);


    header_length=strfind(file,[new1 new2 new1 new2])+3;

    t=strfind(file,'width')+5+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    width=str2num(file(t:t+tt-1));
    
    t=strfind(file,'height')+6+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    height=str2num(file(t:t+tt-1));
    
    t=strfind(file,'frames-generated')+16+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    frames_generated=str2num(file(t:t+tt-1));
    
   % if 

%get frames

    %input file
        fid=fopen(movie_file,'r');
    %scipp header
        fread(fid, header_length);
    

    stim = zeros(samples, prod(xy/zoom));
    for i=1:samples
        t=fread(fid,width*height*3,'ubit8');
        tt=reshape(t,3,width,height);   
        ttt=squeeze(sum(tt([2 3],:,:)))'/2;
        m=ttt(index);

        if ~n==1
            m = reshape(m,[n xy(2)/n xy(1)]);
            m = permute(m,[1 3 2]);
            m = reshape(m,[n*n xy(1)/n xy(2)/n]);
            m = squeeze(sum(m,1))';
        end
        stim(i,:)=m(:);
    end
    
    fclose(fid);
    
disp(sprintf('\tdone (%.2f seconds)\n',toc));


    
    
    
    
function [stim,width,height,frames_generated,header_length]=load_rawmovie(movie_file, frames, varargin)
%load rawmovie file
%
%   stim(height,width,RGB,frames)
%
%   movie_file: full path
%   [frames]: desired number of frames from start, default all, 0 -> only header  
%
%greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('tab', '09');
    p.addParamValue('new1', '0d');
    p.addParamValue('new2', '0a');
    p.parse(varargin{:});
    params = p.Results;
    

%disp('Getting movie');tic

%get header
    tab=hex2dec(params.tab);
    new1=hex2dec(params.new1);
    new2=hex2dec(params.new2);

    fclose all;
    fid=fopen(movie_file,'r');
    file=fscanf(fid, '%c', 8000);
    fclose(fid);

    header_length=strfind(file,[new1 new2 new1 new2])+3;
    if isempty(new1)
        header_length=header_length-2;
    end

    t=strfind(file,'width')+5+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    width=str2num(file(t:t+tt-1));
    
    t=strfind(file,'height')+6+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    height=str2num(file(t:t+tt-1));
    
    t=strfind(file,'frames-generated')+16+1;
    tt=strfind(file(t:t+10),[new1 new2]);
    frames_generated=str2num(file(t:t+tt-1));
    
if ~exist('frames','var');
    frames=frames_generated;
end 


%get frames
if frames>0
    %input file
        fid=fopen(movie_file,'r');
    %scipp header
        fread(fid, header_length);
    

    stim = single(zeros(height,width,3,frames));
    nr=[];
    for i=1:frames
        t=fread(fid,width*height*3,'ubit8');
        if length(t)==width*height*3
            f=reshape(t,3,width,height);
            stim(:,:,:,i)=permute(f,[3 2 1]);
        else
            nr=[nr i];
        end
    end 
    if ~isempty(nr)
        stim=stim(:,:,:,1:nr(1)-1);
        disp(sprintf('\trawmovie has only %d frames\n',nr(1)-1));
    end
       
    fclose(fid);
    
else
    stim=[];
end  

%disp(sprintf('\tdone (%.2f seconds)\n',toc));


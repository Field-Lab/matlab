function [stim,height,width,header_size] = get_raw_movie(movie_file, frames, color, index, zoom)
%greschner

% load movie
fid=fopen(movie_file,'r');

t=fscanf(fid, '%s', 1);
if ~isequal(t,'header-size')
    error('no header-size')
else
    header_size=str2num(fscanf(fid, '%s', 1));
end

height=[];width=[];
while isempty(height)|isempty(width)
    t=fscanf(fid,'%s',1);
    switch t
        case 'height'
         	height=str2num(fscanf(fid,'%s',1));
        case 'width'
            width=str2num(fscanf(fid,'%s',1));
        otherwise
            fscanf(fid,'%s',1);
    end
end

stim=[];
if exist('frames','var') & frames
    T=text_waitbar('get_raw_movie');
    
    color=color/sum(color);
    
    %input file
        fid=fopen(movie_file,'r');
    %scipp header 
        fread(fid, header_size);

    if ~exist('index','var')
        stim = zeros(frames,width,height);
        for i=1:frames
            t=fread(fid,width*height*3,'ubit8');
            tt=reshape(t,3,width,height); 
            %stim(i,:,:)=((tt(1,:,:).*color(1)+tt(2,:,:).*color(2)+tt(3,:,:).*color(3))/sum(color));
            stim(i,:,:)=tt(1,:,:);
        end    
    else    
        stim = zeros(frames, length(index));
        for i=1:frames
            t=fread(fid,width*height*3,'ubit8');
            tt=reshape(t,3,width,height); 
            %ttt=squeeze(sum(tt([2 3],:,:)))'/2;
            m=squeeze(tt(1,:,:).*color(1)+tt(2,:,:).*color(2)+tt(3,:,:).*color(3))';
            mm=m;
            if zoom~=1
                m = reshape(m,[zoom height/zoom width]);
                m = permute(m,[1 3 2]);
                m = reshape(m,[zoom*zoom width/zoom height/zoom]);
                m = squeeze(sum(m,1))';
            end
            stim(i,:)=m(index);
            T=text_waitbar(T,i/frames);
        end
    end
    
    fclose(fid);
end


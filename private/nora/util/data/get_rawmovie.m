function X=get_rawmovie(moviefile, frames, start_frame)
% Warning: can take awhile to load
% Roughly 30  seconds to load up in Bertha for a 30 second movie
% Roughly 90 seconds on Alligator for a 30 second movie

% start frame is indexed from 0, so if you set start_frame = 0, will start
% at the beginning of the movie

% sort through the header to find the movie size
fid = fopen(moviefile,'r');
t = fscanf(fid,'%s',1);
if ~strcmpi(t,'header-size')
    error('no header-size')
else
    header_size = str2double(fscanf(fid, '%s', 1));
end
height = [];
width = [];
while ( isempty(height) || isempty(width) )
    t = fscanf(fid,'%s',1);
    switch t
        case {'height', 'HEIGHT'}
            height = str2double(fscanf(fid,'%s',1));
        case {'width', 'WIDTH'}
            width = str2double(fscanf(fid,'%s',1));
        otherwise
            fscanf(fid,'%s',1);
    end
end

% now it's time to actually read the movie
fid = fopen(moviefile,'r');
fread(fid, header_size); % skip header

if ~exist('X','var')
    X = zeros(frames,width,height,'uint8');
end

% read to the point we actually want 
fread(fid,start_frame*width*height*3,'ubit8');

% Loading up the Raw Movie
for i = 1:frames
    t = fread(fid,width*height*3,'ubit8');  % I think the 3 relates to RGB guns
    if isempty(t)
        error(['Raw movie is only ' num2str(i) ' frames.'])
    end
    tt = reshape(t,3,width,height);
    X(i,:,:) = tt(1,:,:);
end

fclose(fid);

end

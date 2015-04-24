% Check the raw movies!

% Enter the raw movie file, the number of frames to load, and the frame to
% start at here.
moviefile='NSbrownian_6000_A_mean025.rawMovie';
frames=150;
start_frame=1;
% Warning: can take awhile to load
% Roughly 30  seconds to load up in Bertha for a 30 second movie
% Roughly 90 seconds on Alligator for a 30 second movie

% sort through the header to find the movie size
fid = fopen(moviefile,'r');
t = fscanf(fid,'%s',1);
if ~isequal(t,'header-size')
    error('no header-size')
else
    header_size = str2double(fscanf(fid, '%s', 1));
end
height = [];
width = [];
while ( isempty(height) || isempty(width) )
    t = fscanf(fid,'%s',1);
    switch t
        case 'height'
            height = str2double(fscanf(fid,'%s',1));
        case 'width'
            width = str2double(fscanf(fid,'%s',1));
        otherwise
            fscanf(fid,'%s',1);
    end
end

% now it's time to actually read the movie
fid = fopen(moviefile,'r');
fread(fid, header_size); % skip header
framenums = (1:frames)+start_frame-1;

if ~exist('X','var')
    X = zeros(frames,width,height,'uint8');
end

% Loading up the Raw Movie
for i = 1:frames
    f = i+start_frame-1;
    t = fread(fid,width*height*3,'ubit8');  % I think the 3 relates to RGB guns
    tt = reshape(t,3,width,height);
    X(i,:,:) = tt(1,:,:);
end

for i = 1:frames
    imagesc(squeeze(X(i,:,:))')
    colormap gray
    axis image
    pause(0.005)
end

fclose(fid);

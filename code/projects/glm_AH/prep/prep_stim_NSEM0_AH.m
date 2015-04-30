%%% NEW X_RESCALED SHOULD WORK BRING THE NSEM INTO A  -.5 TO .5 SCALE !!!

%% JUST USING EIZABURO'S AND MARTIN'S CODE  DIRECTLY

%the same movie will be used for 2012-04-13 or 14.
%the file is called eye-120-10-3-0_3600.rawMovie
% saving a single long NSEM movie with dozens of MAT files
   % so that MATLAB could handle the movie more efficiently.
% 2012-01-27-4, edoi.

clear
glmhome = '/snle/home/snl-e/glm';

d_mov = '/Volumes/Rat/Data/edoi/2012-01-27-4/stimulus/';
cd(d_mov)
f_mov = 'eye-3600.rawMovie';
% this is 30 sec, 120 frames/sec movie (hence named 3600).
% 30 and 60 sec/block, for repeat and non-repeat trials, resp.
% here, save movie for every 30 sec.
d_sav = sprintf('%s/output_AH_STIMplusminus/2012-09-27-3/NSEM/Movie',glmhome)
if ~exist(d_sav,'dir')
   unix(sprintf('mkdir %s',d_sav));
end

fid = fopen(f_mov,'r'); % load movie
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

%%

fid = fopen(f_mov,'r'); % load movie
fread(fid, header_size); % skip header

nfb = 120*30; % number of frames per block
Xr = nan(3600,80,40); % reduced by 4x4 pixels -> 1 pixel
for k = 1:2*3*600/30
   F = (k-1)*nfb+1:k*nfb;
   %X = get_raw_movie2(f_mov,FR,[1,0,0]); % note: X is uint8
   
   %----------------
   if ~exist('index','var')
      if ~exist('X','var') % initialize
         X = zeros(nfb,width,height,'uint8');
      end
      for i = 1:length(F)
         f = F(i);
         t = fread(fid,width*height*3,'ubit8');
         tt = reshape(t,3,width,height);
         X(i,:,:) = tt(1,:,:);
      end
   else
      1;
   end
   %----------------
   
   for f = 1:nfb
      tmp = double(squeeze(X(f,:,:)));
      for i = 1:80
         for j = 1:40
            tmp2 = tmp(4*(i-1)+1:4*i,4*(j-1)+1:4*j);
            Xr(f,i,j) = sum(tmp2(:))/16;
         end
      end
      if rem(f,360) == 0
         fprintf('now%d done\r',f)
      end
   end
   X_rescaled = ( (Xr-1) /255 ) - .5; 
   eval(sprintf('save %s/eyemov%d X Xr X_rescaled',d_sav,k))
   fprintf('\rNow %d/%d done\n',k,2*3*600/30)
end

fclose(fid);

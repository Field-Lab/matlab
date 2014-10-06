function mm = makemovie(framevec,title,fps,reps,movieformat,itwm_logo,blockpat)

% function mm = makemovie(framevec,title[,framespersec[,reps[,movieformat]]])
%
% makemovie takes the files __MovieD/<framevec(:)>.ppm (usually) created
% with makeframe and writes them either into a file <title>.mpg, or
% into a file <title>.gif. The first is an MPEG animation, the second a
% gif89a animation. 
%
% Selectable movieformats are 'gif', 'mpg' or 'all', movies are written
% accordingly. On exit, the __MovieD and all files in it are removed.
%
% For the gif animation the framespersec argument specifies the
% frames per second (default 10) and the reps argument the
% repetitions (default 1). Gif animations can be viewed with
% Mediaplayer, NETSCAPE, animate, etc.
%	
% If called without arguments, makemovie returns the directory that the
% movie package is installed in.
%
% To create multiple movies in parallel, create a subdirectory
% for each movie in which you execute the makeframe and makemovie
% commands.
%
% See also makeframe, moviedemo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by ALW 11/24/95                              %
% updated to Matlab 5 by Don Bovee 4/97                %
%                                                      %
%      modified by ALW, 05/15/98 for Win95 & NT        %
%      modified by ALW, 03/22/00 for Linux             %
%      modified by ALW, 02/19/02 for multiple movies   %
%      modified by ALW, 10/09/02 for itwm logo options %
%                                and correct gzip use  %
% last modified by ALW, 11/07/02 use user dependent    %
%                                temporary directory   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Site dependent settings -------------------------------------------------
%

% Linux/Unix shared directory
packagedir = '/homes/wiegmann/pc/matlab/Movie/';

%
% No more editing needed beyond this point ---------------------------------
%

% ---------------------------------------------------
% however, you may want to include your architecture
% and or system paths here.
% --------------------------------------------------
comp = computer;
if    sum(comp(1:3) ~= 'ALP') & ...
      sum(comp(1:3) ~= 'SGI') & ...
      sum(comp(1:3) ~= 'HP7') & ...
      sum(comp(1:3) ~= 'SUN') & ...
      sum(comp(1:3) ~= 'SOL') & ...
      sum(comp(1:3) ~= 'LNX') & ...
      sum(comp(1:3) ~= 'PCW')
  disp('makemovie currently set up for SUN SPARC, SUN SOLARIS,')
  disp('Intel PCs (Win 95, NT and Linux), DEC ALPHA, SGI and HP only')
  mm = -1 ; error = -2;
  return
  % windows differs from UNIX: \/, command names
elseif sum(comp(1:3) == 'PCW'),
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%   Windows case  %%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Special network path for ALW at ITWM, 10/9/02
  [s,w] = unix('set COMPUTERNAME');
  if length(w) == 20
    if w(14:19)== 'HYDRA2'
      packagedir = '\\fs1pre\wiegmann\matlab\Movie\';
    end    
  else
    % Standard Windows path as described in documentation
    packagedir = 'c:\Matlab\Movie\';
  end
  Mdir1     = getenv('TEMP');
  if ~isempty(Mdir1)
    if exist(Mdir1,'dir')==7
      Mdir1 = [Mdir1,'\'];
    else
      Mdir1 = [];
    end
  end
  encoder    = [packagedir,computer,'\mpeg_encode ']; % -float_dct
  %merge      = [packagedir,computer,'\gifmerge -255,255,255 '];
  merge      = [packagedir,computer,'\gifmerge '];
  ppmtogif   = [packagedir,computer,'\ppmtogif '];
  giftopnm   = [packagedir,computer,'\giftopnm '];
  gunzip     = [packagedir,computer,'\gunzip -c '];
  rm         = ['del  '];  % let string end in space!
  mv         = ['move '];  % let string end in space!

	% Enable gunzip as an option in the future!
	% ALW, 8/22/03
  gunzip = 'type ';
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%   Unix/Linux case  %%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Mdir1       = '/var/tmp';

  if exist(Mdir1,'dir')==7
		[s,user] = unix('whoami');
		% Attention: last character in user is carriage return!
    Mdir1    = [Mdir1 '/' user(1:end-1)];
		if exist(Mdir1,'dir')~=7
			[s,w] = unix(['mkdir ' Mdir1]);
			if s~=0
			  Mdir1 = pwd;
			end
		end
		Mdir1 = [Mdir1 '/'];
  else
    Mdir1 = [];
  end

  encoder    = [packagedir,computer,'/mpeg_encode ']; % -float_dct
  [s,w] = unix(encoder); % New Linux library needed
  if s == 128
    if sum(comp(1:3) == 'LNX')
      encoder    = [packagedir,computer,...
                    '/mr_mpeg_encode '];
    else
      disp(M2);
    end
  end
  %merge      = [packagedir,computer,'/gifmerge -255,255,255 '];
	merge      = [packagedir,computer,'/gifmerge '];
  ppmtogif   = [packagedir,computer,'/ppmtogif '];
  giftopnm   = [packagedir,computer,'/giftopnm '];
  if exist('/usr/bin/gunzip')
    gunzip   = ['/usr/bin/gunzip -c ']; % let string end in space!
  elseif exist('/bin/gunzip')
    gunzip   = ['/bin/gunzip -c ']; % let string end in space!
  elseif exist('/usr/ucb/gunzip')
    gunzip   = ['/usr/ucb/gunzip -c ']; % let string end in space!
  else disp('Having trouble finding gunzip, proceed with fingers crossed')
    gunzip   = ['gunzip -c '];      % let string end in space!
  end
  if exist('/usr/bin/rm')
    rm         = ['/usr/bin/rm ']; % let string end in space!
  elseif exist('/bin/rm')
    rm         = ['/bin/rm ']; % let string end in space!
  elseif exist('/usr/ucb/rm')
    rm         = ['/usr/ucb/rm ']; % let string end in space!
  else disp('Having trouble finding rm, proceed with fingers crossed')
    rm         = ['rm '];      % let string end in space!
  end
  if exist('/usr/bin/mv')
    mv         = ['/usr/bin/mv ']; % let string end in space!
  elseif exist('/bin/mv')
    mv         = ['/bin/mv ']; % let string end in space!
  elseif exist('/usr/ucb/mv')
    mv         = ['/usr/ucb/mv ']; % let string end in space!
  else disp('Having trouble finding mv, proceed with fingers crossed')
    mv         = ['mv '];      % let string end in space!
  end

	% Enable gunzip as an option in the future!
	% ALW, 8/22/03
  gunzip = 'cat ';

end %if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REALLY NO EDITING needed beyond this point ---------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 0, mm = packagedir; return; end;
M1 = ['makemovie: error --- could not open movie directory.   '
      '           please run makeframe first                  '
      '           returning to calling routine, no movie made.'];
M2 = ['Libraries missing in makemovie for mpeg_encoder.'
      'Try to use old mpeg_encoder in LNX86 directory! '];

if exist('fps'),
  if fps == 'all' | fps == 'gif' | fps == 'mpg',
    tmp = fps;
    if exist('reps')
      fps = reps;
      if exist('movieformat')
        reps = movieformat;
      end
    else
      clear fps
    end
    movieformat = tmp;
  end
end

if exist('fps'),
  delay = round(100/fps);
  if delay ~= 100/fps,
    disp(['warning: requested frames per second ',...
    num2str(fps),' changed to ',num2str(100/delay)])
  end
else delay = 10;
end;
if ~exist('reps'),
  reps = 1;
else tempreps = reps;
  reps = round(reps);
  if tempreps ~= reps
    disp(['warning: requested repetitions ',...
    num2str(tempreps),' changed to ',num2str(reps)])
  end;
end

merge = [merge, ' -',num2str(delay),' -l',num2str(reps)];

% NEED SAME Mdir as in makeframe.m
Mdir   = [Mdir1,'__',title(1:min(6,length(title)))];
Odir   = pwd;
if sum(comp(1:3) == 'PCW'),
  Odir = [Odir,'\'];
else
  Odir = [Odir,'/'];
end
mm     = 0;
errors = 0;

framevec = framevec(:);
images   = size(framevec,1);

patterncounter = images;
pattern        = '';
if exist('blockpat')~=1
        blockpat       = 'IBBPBBPBBPBBPBBPBB';
end
lbp            = length(blockpat);
restpat        = blockpat(1:lbp-1);
while patterncounter  > lbp,
  pattern = [pattern,blockpat];
  patterncounter= patterncounter -lbp;
end
pattern = [pattern,restpat(1:patterncounter-1),'I'];

titleg  = [title,'.gif'];
title   = [title,'.mpg'];

if ~exist('movieformat'),
  movieformat = 'all';
else
  if movieformat ~= 'gif' & movieformat ~= 'mpg' & movieformat ~= 'all',
    error(['unknown format: ',movieformat])
  end
end

if movieformat == 'all' | movieformat == 'mpg',
  if size(dir(title),1)
    disp(['File ',title,' exists, '...
    'is renamed to ',title,'.orig'])
    if size(dir([title,'.orig']),1)
      unix([rm,' ', title,'.orig']);
    end
    unix([mv,' ',Odir,title,' ',Odir,title,'.orig']);
  end
end

if movieformat == 'all' | movieformat == 'gif',
  frames = '';
  for frame = framevec'
    framestr = num2str(frame);
    frames = [frames,' ',framestr];
  end % if
  if size(dir(titleg),1)
    disp(['File ',titleg,' exists, '...
    'is renamed to ',titleg,'.orig'])
    if size(dir([titleg,'.orig']),1)
      unix([rm,' ', titleg,'.orig']);
    end
    unix([mv,' ',Odir,titleg,' ',Odir,titleg,'.orig']);
  end
end

eval(['cd ',Mdir],['errors = 1;']);
if errors == 1,
  disp(M1);
  return;
end;

if movieformat == 'all' | movieformat == 'gif'
  [s,w] = unix([merge,' ',frames,' > ',titleg]);
  unix([mv,' ',titleg,' ',Odir]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The parameters for the MPEG movie are
% written (and could be changed) here!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if movieformat == 'all' |movieformat == 'mpg'
  load width
  load height
  fp = fopen('movie.par','w');
  fprintf(fp,'# parameter file for mpeg movie\n\n');
  fprintf(fp,'PATTERN              %s\n',pattern);
  fprintf(fp,'OUTPUT               %s\n\n',title);
  fprintf(fp,'YUV_SIZE             %dx%d\n\n',width,height);
  fprintf(fp,'BASE_FILE_FORMAT     PNM\n');
	
	inputconv = [gunzip ' *']
	[s,user] = unix('whoami');
	if length(user)==9
		if user(1:8)=='wiegmann';
			inputconv = [gunzip ' * | pnmcomp -align right -valign top -alpha ',...
									 packagedir,'AWiegmann.pgm ',packagedir,'AWiegmann.ppm ']
		end
	end
	
  if exist('itwm_logo') == 1
    switch itwm_logo
     case +1
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGO.pnm 1 1 \n'],inputconv);
     case +2
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGOs.pnm 1 1 \n'],inputconv);
     case +3
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGOss.pnm 1 1 \n'],inputconv);
     case -1
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGO.pnm -89 -90 \n'],inputconv);
     case -2
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGOs.pnm -44 -45 \n'],inputconv);
     case -3
      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
                  packagedir,'ITWMLOGOss.pnm -22 -23 \n'],inputconv);
     otherwise
			fprintf(fp,'INPUT_CONVERT        %s *\n',inputconv);
%      fprintf(fp,['INPUT_CONVERT        %s | pnmpaste -replace ',...
%                  packagedir,'ITWMLOGO.pnm -89 -90 \n'],inputconv);
    end
  else
    fprintf(fp,'INPUT_CONVERT        %s *\n',inputconv);
  end
  fprintf(fp,'GOP_SIZE             30\n');
  fprintf(fp,'SLICES_PER_FRAME     1\n\n');
  fprintf(fp,'INPUT_DIR            .\n\n');
  fprintf(fp,'INPUT\n');
  for i = 1:images
		% Enable gunzip as an option in the future!
		% ALW, 8/22/03
		% fprintf(fp,'%s.ppz\n',num2str(framevec(i)));
		fprintf(fp,'%s.ppm\n',num2str(framevec(i)));
  end
  fprintf(fp,'END_INPUT            \n\n');
  fprintf(fp,'PIXEL                HALF\n');
  fprintf(fp,'RANGE                15\n\n'); % was: 10
  fprintf(fp,'PSEARCH_ALG          LOGARITHMIC\n');   % was: EXHAUSTIVE
  fprintf(fp,'BSEARCH_ALG          CROSS2\n\n'); % was: EXHAUSTIVE
  fprintf(fp,'IQSCALE              1\n');   % was:  8
  fprintf(fp,'PQSCALE               1\n');  % was: 10
  fprintf(fp,'BQSCALE               1\n\n');% was: 25
  fprintf(fp,'REFERENCE_FRAME      DECODED\n');
  fprintf(fp,'FORCE_ENCODE_LAST_FRAME\n');
  fclose(fp);
	
  [s,w] = unix([encoder,' movie.par > mpg.out']);
	disp(w)
  unix([mv,' ',title,' ',Odir]);
end

% clean up __MovieD
if 0 % disable by setting to 0 if desired for debugging puposes
  files = dir;
  for i = 1:size(files,1)
    if ~ files(i).isdir,
      unix([rm ,files(i).name]);
    end
  end;
  chdir('..');
  unix(['rmdir ',Mdir]);
end
chdir(Odir)

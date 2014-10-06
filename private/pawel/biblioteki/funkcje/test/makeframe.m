function mf = makeframe(frame,title,geom,pmethod,q)

% function mf = makeframe(frame,title[,geom,[pmethod,[q]]])
%
% Takes the figure in the currently active window, and writes
% it into a file __title(1:min(6,length(title)))/<frame> inside a
% temporary directory determined from getenv('TEMP') under windows
% and  in /var/tmp under UNIX/Linux if they exist, otherwise in the
% current directory.
%
% The size of the images defaults to 240x192 pixels, unless
% width and height are specified in the optional vector geom.
% For best results, the user should make sure that all images
% are on the same scale, though technically it is not necessary.
% For best compression, image sizes should be multiples of 8
% (pixels).
%
% The image generation method is declared by the optional argument
% pmethod. By default, a method is used that creates inferior
% quality images, but works even if the user obscurs the figure window.
% For matlab versions 5.3 and higher pmethod can be set to 1 and the
% getframe command is used to generate exactly the image as seen on
% screen. This method has difficulties when a user switches screens
% away from the matlab figure. The latter feature is unavailable on 
% SUN and SGI machines, if you wish to use it please contact 
%
% andreas@familie-wiegmann.de
%	
% If the optional argument q is given then the gif animation is
% quantified to use at most q colors. This is neccessary because
% the gif file format supports only 256 colors.
%
% See also makemovie.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by ALW 11/24/95                              %
% updated to Matlab 5 by Don Bovee 4/97                %
%                                                      %
%      modified by ALW, 05/15/98 for Win95 & NT        %
%      modified by ALW, 03/22/00 for Linux             %
%      modified by ALW, 02/19/02 for multiple movies   %
%      modified by ALW, 10/09/02 for zbuffer rendering %
%                                under windows         %
%      modified by ALW, 10/26/02 for tifftopnm use     %
%      modified by ALW, 11/07/02 use user dependent    %
%                                temporary directory   %
% last modified by ALW, 06/08/03 default to old print  %
%                                mode                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Site dependent settings -------------------------------------------------
%

% Linux/Unix shared directory
packagedir = '/homes/wiegmann/pc/matlab/Movie/';

%
% No editing needed beyond this point ---------------------------------
%
% ---------------------------------------------------
% however, you may want to include your architecture
% and or system paths here.
% ---------------------------------------------------
comp = computer;
if    sum(comp(1:3) ~= 'ALP') & ...
      sum(comp(1:3) ~= 'SGI') & ...
      sum(comp(1:3) ~= 'HP7') & ...
      sum(comp(1:3) ~= 'SUN') & ...
      sum(comp(1:3) ~= 'SOL') & ...
      sum(comp(1:3) ~= 'PCW') & ...
      sum(comp(1:3) ~= 'LNX')
  disp('makeframe currently set up for SUN SPARC, SUN SOLARIS,')
  disp('Intel PCs (Win 95, NT, Linux), DEC ALPHA, SGI and HP only')
  mm = -1 ; error = -2;
  return
elseif sum(comp(1:3) == 'PCW')
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%% PCWin special case %%%%%%%%%%%%%%%
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

  Mdir1       = getenv('TEMP');
  if ~isempty(Mdir1)
    if exist(Mdir1,'dir')==7
      Mdir1 = [Mdir1,'\'];
    else
      Mdir1 = [];
    end
  end
  pnmcut     = [packagedir,computer,'\pnmcut'];
  ppmtogif   = [packagedir,computer,'\ppmtogif'];
  ppmquant   = [packagedir,computer,'\ppmquant'];
  gzip       = [packagedir,computer,'\gzip'];
  tifftopnm  = [packagedir,computer,'\tifftopnm'];
  rm         = ['del '];
else
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%   Unix/Linux case  %%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('q') ...
        & sum(comp(1:3) ~= 'SOL') ...
        & sum(comp(1:3) ~= 'SUN') ...
        & sum(comp(1:3) ~= 'ALP') ...
        & sum(comp(1:3) ~= 'HP7') ...
        & sum(comp(1:3) ~= 'LNX')
    disp (['Quantification not yet implemented for ',comp]);
    q = 0;
  end
	 Mdir1       = '/var/tmp';

  if exist(Mdir1,'dir')==7
		[s,user] = unix('whoami');
		% Attention: last character in user is carriage return!
    Mdir1    = [Mdir1 '/' user(1:end-1)];
		if exist(Mdir1,'dir')~=7
			[s,w] = unix(['mkdir ' Mdir1]);
			if s ~= 0
			  Mdir1 = pwd;
			end
		end
		Mdir1 = [Mdir1 '/'];
  else
    Mdir1 = [];
  end

  pnmcut     = [packagedir,computer,'/pnmcut'];
  if exist('/usr/bin/pnmcut')
    pnmcut   = ['/usr/bin/pnmcut ']; % let string end in space!
  end
	ppmtogif   = [packagedir,computer,'/ppmtogif'];
  if exist('/usr/bin/ppmtogif')
    ppmtogif = ['/usr/bin/ppmtogif ']; % let string end in space!
  end
  ppmquant   = [packagedir,computer,'/ppmquant'];
  if exist('/usr/bin/ppmquant')
    ppmquant = ['/usr/bin/ppmquant ']; % let string end in space!
  end
  if exist('/usr/bin/tifftopnm')
    tifftopnm  = ['/usr/bin/tifftopnm ']; % let string end in space!
  elseif exist('/bin/tifftopnm')
    tifftopnm  = ['/bin/tifftopnm ']; % let string end in space!
  elseif exist('/usr/ucb/tifftopnm')
    tifftopnm  = ['/usr/ucb/tifftopnm ']; % let string end in space!
  else
    disp('Cannot find tifftopnm, using older print -zbuffer version')
  end

  if exist('/usr/bin/gzip')
    gzip   = ['/usr/bin/gzip ']; % let string end in space!
  elseif exist('/bin/gzip')
    gzip   = ['/bin/gzip ']; % let string end in space!
  elseif exist('/usr/ucb/gzip')
    gzip   = ['/usr/ucb/gzip ']; % let string end in space!
  else disp('Having trouble finding gzip, proceed with fingers crossed')
    gzip   = ['gzip '];      % let string end in space!
  end

  if exist('/usr/bin/mv')
    mv         = ['/usr/bin/mv ']; % let string end in space!
  elseif exist('/usr/ucb/mv')
    mv         = ['/usr/ucb/mv ']; % let string end in space!
  elseif exist('/bin/mv')
    mv         = ['/bin/mv ']; % let string end in space!
  else disp('Having trouble finding mv, proceed with fingers crossed')
    mv         = ['mv '];      % let string end in space!
  end
  if exist('/usr/bin/rm')
    rm         = ['/usr/bin/rm ']; % let string end in space!
  elseif exist('/usr/ucb/rm')
    rm         = ['/usr/ucb/rm ']; % let string end in space!
  elseif exist('/bin/rm')
    rm         = ['/bin/rm ']; % let string end in space!
  else disp('Having trouble finding rm, proceed with fingers crossed')
    rm         = ['rm '];      % let string end in space!
  end
end %if

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REALLY NO EDITING needed beyond this point ---------------------------------
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEED SAME Mdir as in makemovie.m

Mdir       = ['__',title(1:min(6,length(title)))];
firstframe = (mkdir(Mdir1,Mdir) == 1);
Mdir       = [Mdir1,Mdir];
Odir       = pwd;
eval(['cd ',Mdir]);

framest = num2str(frame);

% default, possible to input these
width = 240;
height= 192;

if exist('geom'),
  if (size(geom) ~= 2)
    disp(['makeframe: incorrect geometry input size ',...
      num2str(size(geom))]);
    error = -4;
    chdir('..');
    unix([rm,' -R ',Mdir]);
    chdir(Odir);
    return
  else
    if geom(1) <= 0
      disp(['makeframe: incorrect geometry input (x non-positive): ',...
        num2str(geom(1))]);
      error = -8;
      chdir('..');
      unix([rm,' -R ',Mdir]);
      chdir(Odir);
      return
    end
    if geom(2) <= 0
      disp(['makeframe: incorrect geometry input (y non-positive): ',...
        num2str(geom(2))]);
      error = -8;
      chdir('..');
      unix([rm,' -R ',Mdir]);
      chdir(Odir);
      return
    end
  end
  width  = geom(1);
  height = geom(2);
end;
if ~exist('pmethod')
  pmethod = 0;
end
if ~exist('q')
  q = 0;
end

OLDVERSION = 5.2;
VERSION    = version;
VERSION    = str2num(VERSION(1:3));

if pmethod
	if OLDVERSION >= VERSION
		OLDMATLAB = 1;
	else    
		OLDMATLAB = 0;
	end
	if (all(comp(1:3) == 'SUN') | ...
	    all(comp(1:3) == 'SOL'))
		OLDMATLAB = 1;
	end	
else
	OLDMATLAB = 1;
end

if (exist('tifftopnm') ~= 1) | OLDMATLAB
  oldPaperUnits    = get(gcf,'PaperUnits');
  oldPaperPosition = get(gcf,'PaperPosition');

  set(gcf,'PaperUnits','points')
  PaperSize        = get(gcf,'PaperSize'); %size needed in points
  set(gcf,'PaperPosition',[0 PaperSize(2)-height width height])

  %flush the event queue and draw current figure, then print
  drawnow
  print -zbuffer -dppmraw image.ppm

  %restore settings for current figure
  set(gcf,'PaperUnits',   oldPaperUnits)
  set(gcf,'PaperPosition',oldPaperPosition)
else
  oldPosition = get(gcf,'Position');
  set(gcf,'Position',[oldPosition(1:2) width height]);
  figure(gcf); % ensure getframe gets whole figure window
  gf  = getframe(gcf);
  gf3 = gf.cdata;
  imwrite(gf3,'image.tif','tif');
  clear gf gf3 cmap
  [s w] = unix([tifftopnm,' image.tif > image.ppm']);
  unix([rm,' image.tif']);
  set(gcf,'Position',oldPosition);
  drawnow
end

geom = [' 0 0 ',int2str(width),' ',int2str(height),' '];

if exist('firstframe')
  save width width -ascii
  save height height -ascii
  [s,w] = unix(['echo "',num2str(width),'x',...
    num2str(height),'" > firstframesize']);
  [s,w] = unix(['echo "',num2str(width),' ',...
    num2str(height),'" > framesize']);
else
  [s,w] = unix(['echo "',num2str(width),' ',...
    num2str(height),'" >> framesize']);
end

outf = [framest];
if q
  [s,w] = unix([...
    ppmquant,' ',num2str(q),' image.ppm | ',  ...
    pnmcut,geom,' | ',...
    ppmtogif,' > ',outf]);
else
  [s,w] = unix([...
    pnmcut,geom,' image.ppm | ',...
    ppmtogif,' > ',outf]);
end
% Enable gzip as an option in the future!
% ALW, 8/22/03
%[s,w] = unix([pnmcut,geom,' < image.ppm | ',...
%              gzip,' -f - >', outf,'.ppz']);
[s,w] = unix([pnmcut,geom,' < image.ppm >', outf,'.ppm']);
[s,w] = unix([rm,' image.ppm ']);
chdir(Odir);
mf = 0;

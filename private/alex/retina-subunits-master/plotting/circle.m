function circle(xPos,yPos,sz,FaceColor,EdgeColor,LineWidth);

%-------------------------------------------
%
% h = circle(xPos,yPos,sz,FaceColor,EdgeColor,LineWidth)
%
% wrapping for matlab's "rectangle" for plotting
% a circle centered on (xPos,yPos)
%
% required inputs:
% xPos -- location of center in x
% yPos -- location of center in y
% sz -- size of circle
%
% optional inputs:
% FaceColor -- fill color of circle
% EdgeColor -- color of edge line
% LineWidth -- width of edge line
%
% freeman, 7-12-2012
%-------------------------------------------

if nargin < 3
  error('(circle) must provide position and size\r');
end

if ~exist('FaceColor','var') || isempty(FaceColor)
  FaceColor = 'w';
end

if ~exist('EdgeColor','var') || isempty(FaceColor)
  EdgeColor = 'k';
end

if ~exist('LineWidth','var') || isempty(FaceColor)
  LineWidth = 5;
end

hDot = rectangle('Position',[xPos-sz/2 yPos-sz/2 sz sz],'Curvature',[1 1]);

set(hDot,'FaceColor',FaceColor);
set(hDot,'LineWidth',LineWidth);
set(hDot,'EdgeColor',EdgeColor);
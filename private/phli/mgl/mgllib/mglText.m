% mglText.m
%
%        $Id: mglText.m 138 2007-01-17 15:49:49Z justin $
%      usage: mglText('string')
%         by: justin gardner
%       date: 05/10/06
%  copyright: (c) 2006 Justin Gardner, Jonas Larsson (GPL see mgl/COPYING)
%    purpose: Returns image of string.
%
%       e.g.: 
%
%mglOpen
%mglVisualAngleCoordinates(57,[16 12]);
%mglTextSet('Helvetica',32,[1 1 1],0,0,0);
%thisText = mglText('Hello')
%mglBltTexture(thisText,[0 0],'left','top');
%mglFlush;


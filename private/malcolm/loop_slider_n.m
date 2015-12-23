function b1 = loop_slider_n(k,kmin,kmax,n)
%LOOP_SLIDER_N is a quick and easy way to make a simple slider at the bottom of your plot
%   LOOP_SLIDER_N differs from loop_slider by allowing you to have more
%   than one slider on a plot.  loop_slider(k,kmin,kmax) is the same as
%   loop_slider_n(k,kmin,kmax,1)
%       k -> starting value for slider
%       kmin -> minimum value for slider
%       kmax -> maximum value for slider
%       n -> the number of the loop slider, 1 for the first loop slider, 2
%            for the second, etc.
%
%   For Example:
%      b1 = loop_slider(1,kmin,kmax,1);
%      k=get(b1,'Value')
%      b2 = loop_slider(1,jmin,jmax,2);
%      j=get(b2,'Value')
%
% Marvin Thielk 2013
% mthielk@salk.edu

b1= uicontrol(gcf,...
    'Style','slider',...
    'Min' ,kmin,'Max',kmax, ...
    'Position',[200*(n-1),0,200,15], ...
    'Value', k,...
    'SliderStep',[1/(kmax-kmin) 1/(kmax-kmin)],...
    'CallBack', 'uiresume;');
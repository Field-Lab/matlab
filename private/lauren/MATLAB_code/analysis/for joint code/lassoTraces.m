function [selectedTraces]=lassoTraces(traces)

% traces must be a 2D array of nTraces x nSamples
% lassoTraces -  enables the selection/encircling parts of traces in a plot by hand 
%          using the mouse
% 
% Input:    traces              - a set of traces in an array (traces x samples).
% Output:   selectedTraces      - a set of selected traces in 1 column vector
% 
% Note:   After the plot of traces is given, selection by mouse is started immediately.
%         Encircling is done by pressing subsequently the LEFT button mouse at the requested positions 
%         in a scatter plot.
%         Closing the loop is done by a RIGHT button press.
%         
% original lasso.m written by T.Rutten V2.0/9/2003
% modified by Lauren Hruby, SNL-E (lhruby@salk.edu)
% 2008-10-24

nTraces = size(traces, 1);
nSamples = size(traces, 2);


traceColor = hsv(size(traces,1));

hold on
for i = 1:nTraces
    current = plot(traces(i,:));
    set(findobj(current,'Type','line'), 'Color', traceColor(i,:))
end
hold off

las_x=[];
las_y=[];

c=1;

key=0;

% disp('press a KEY to start selection by mouse, LEFT mouse button for selection, RIGHT button closes loop')
% while key==0
% key=waitforbuttonpress;
% pause(0.2)
% end

while c==1 

[a,b,c]=ginput(1);
las_x=[las_x;a];las_y=[las_y;b];
line(las_x,las_y)
end;

las_x(length(las_x)+1)=las_x(1); %"closes loop" by making adding first coordinates to end
las_y(length(las_y)+1)=las_y(1);

line(las_x,las_y) % redraws closed loop
pause(.2)

selectedTracesIndeces = zeros(1, nTraces);
for i = 1:nTraces
    in=inpolygon((1:nSamples),traces(i,:),las_x,las_y); %matrix of zeros and ones specifying which points are in polygon
    selectedTracesIndeces(i) = max(in) > 0;
end

ev_in=find(selectedTracesIndeces>0);

% selx=x(ev_in);
% sely=y(ev_in);
% 
% figure,plot(x,y,'b.',selx,sely,'g.');
% legend(num2str([length(x)-length(selx);length(selx)]));

selectedTraces=ev_in;
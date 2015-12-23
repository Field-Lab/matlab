function [selectedPoints]=lassoPoints(points)

% traces must be a matrix of nPoints x 2
% lassoPoints -  enables the selection/encircling points in a plot by hand 
%          using the mouse
% 
% Input:    points              - a set of x,y coordinates traces in an array (points x 2).
% Output:   selectedTraces      - a set of selected traces in 1 column vector
% 
% Note:   After the plot of points is given, selection by mouse is started immediately.
%         Encircling is done by pressing subsequently the LEFT button mouse at the requested positions 
%         in a scatter plot.
%         Closing the loop is done by a RIGHT button press.
%         
% original lasso.m written by T.Rutten V2.0/9/2003
% modified by Lauren Hruby, SNL-E (lhruby@salk.edu)
% 2008-10-24


plot(points(:,1), points(:,2), 'k.')

las_x=[];
las_y=[];

c=1;

key=0;

while c==1 

[a,b,c]=ginput(1);
las_x=[las_x;a];las_y=[las_y;b];
line(las_x,las_y)
end;

las_x(length(las_x)+1)=las_x(1); %"closes loop" by making adding first coordinates to end
las_y(length(las_y)+1)=las_y(1);

line(las_x,las_y) % redraws closed loop
pause(.2)


in = inpolygon(points(:,1), points(:,2), las_x, las_y); %binary matrix specifying which points are in polygon


selectedPoints = find(in);

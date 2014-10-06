function dl = dotted_line(points,num_dots)
% dotted_line     interpolate values between two given end points
%                   also useful for generating color gradients
%
% usage:  dl = dotted_line(points,num_dots)
%
% arguments:     points - 2xN matrix, each line a point in N-dim space
%              num_dots - number of points to interpolate
%
% outputs:     dl - (num_dots x N) matrix of interpolated points
%                       dl(1,:) = points(1,:)
%                       dl(2,:) = points(2,:)
%
% 2010-01  gauthier
% 2010-01  phli - switched to use LINSPACE; not sure how efficiency
% compares, should be decent and handled the num_dots == 1 situation more
% smoothly.



% get dimensionality
N = size(points,2);

% initialize
dl = zeros(num_dots,N);

% for each dimension...
for nn=1:N
    % get difference
    d = points(2,nn)-points(1,nn);
    
    % if there is any difference
    if d ~= 0
        % interpolate points
%        dl(:,nn) = points(1,nn):(d/(num_dots-1)):points(2,nn);
        dl(:,nn) = linspace(points(1,nn), points(end, nn), num_dots);
    else
        % otherwise, just replicate
        dl(:,nn) = repmat(points(1,nn),[num_dots 1]);
    end
end
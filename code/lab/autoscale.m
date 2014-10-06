function [xrange yrange]=autoscale(xdata,ydata,varargin)
% sets XLim and YLim in current axes
%
% greschner

% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('outlier', .0);%ignor outlier
    p.addParamValue('border', .15);%border region
    p.parse(varargin{:});
    params = p.Results;
    

if ~isempty(xdata);    
    data=sort(xdata);
    temp=ceil(length(data)*[params.outlier 1-params.outlier]+eps);
    temp=data(temp);

    xrange=[temp(1)-(temp(2)-temp(1))*params.border temp(2)+(temp(2)-temp(1))*params.border];
    set(gca,'XLim',xrange);
else
    xrange=[];
end
 

if exist('ydata','var') && ~isempty(ydata);
    data=sort(ydata);
    temp=ceil(length(data)*[params.outlier 1-params.outlier]+eps);
    temp=data(temp);

    yrange=[temp(1)-(temp(2)-temp(1))*params.border temp(2)+(temp(2)-temp(1))*params.border];
    set(gca,'YLim',yrange);
else
    yrange=[];
end
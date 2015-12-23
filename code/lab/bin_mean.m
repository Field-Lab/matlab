function [r_mean r_bin_center r_std]=bin_mean(x,y,varargin)
%bins y based on x, returns mean and center of bins
%vector only
%
%             params    
%               fixed_nr_in_bin default=0  assums gaussian distribution
%               robust_mean default=1
%               robust_std default=1
%               edges default=[]
%               nr_bins default=30
%               ignor_percent_outlier default=0
%
% greschner


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('fixed_nr_in_bin', 0);% 
    p.addParamValue('nr_bins', 25);%
    p.addParamValue('gaussian', 0);%
    p.addParamValue('robust_mean', 1);% 
    p.addParamValue('robust_std', 1);% 
    p.addParamValue('edges', []);%
    p.addParamValue('ignor_percent_outlier', 0);%
  
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

    

    
%get edges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if ~params.fixed_nr_in_bin
    
    if isempty(params.edges)
        t=sort(x);
        tt=length(t)*params.ignor_percent_outlier;
        params.edges=[t(1+tt):(t(end-tt)-t(1+tt))/params.nr_bins:t(end-tt)];
    end        
    r_bin_center=(params.edges(1:end-1)+params.edges(2:end))/2;
    
else
  
    if params.gaussian
        params.nr_bins=params.nr_bins+mod(params.nr_bins,2);
        params.edges=norminv([0:1/params.nr_bins:1],0,std(x));
        r_bin_center=norminv([0:1/params.nr_bins:1-1/params.nr_bins]+1/params.nr_bins/2,0,std(x));
    else
        t=sort(x);
        tt=length(x)/params.nr_bins;
        
        e=[1 : (length(x)/2)/(params.nr_bins/2) : length(x)/2];
        e=[e length(x)/2 : (length(x)-(length(x)/2))/(params.nr_bins/2) : length(x)];
        e=round(unique(e));
        
        params.edges=t(round(e));
        r_bin_center=(params.edges(1:end-1)+params.edges(2:end))/2;
    end
    
end         
    


r_mean=zeros(length(params.edges)-1,1);
r_std=zeros(length(params.edges)-1,1);

for i=1:length(params.edges)-1
    t=find(x>=params.edges(i) & x<params.edges(i+1));

    if params.robust_mean
        r_mean(i)=robust_mean(y(t));
    else
        r_mean(i)=mean(y(t));        
    end

    if params.robust_std
        r_std(i)=robust_std(y(t));
    else
        r_std(i)=std(y(t));        
    end
end
    
    
       
 
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
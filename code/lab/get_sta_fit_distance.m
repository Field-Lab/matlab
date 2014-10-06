function [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun, cell_specification_1, cell_specification_2, cell_specification_norm)
%returns matrix of cell_specification_1 x cell_specification_2 with the
%distance
%
%  [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun,
%  cell_specification_1, cell_specification_2)
%  [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun,
%  cell_specification_1)
%
%   dist_stixel in stixel
%   dist_sd     in mean sd
%   normdist    normalized by sd along line between centers 
%
% see sta_fit_distance
%
% greschner

% get cell numbers
index_1=get_cell_indices(datarun,cell_specification_1);
if ~exist('cell_specification_2','var');
    index_2=index_1;
else
    index_2 = get_cell_indices(datarun,cell_specification_2);
end



dist_stixel=zeros(length(index_1),length(index_2));
dist_sd=zeros(length(index_1),length(index_2));
normdist=zeros(length(index_1),length(index_2));
  
if isequal(index_1, index_2)
    for i=1:length(index_1)
        for ii=i+1:length(index_2)
            [dist_stixel(i,ii) normdist(i,ii)] = sta_fit_distance(datarun.(datarun.default_sta_fits).sta_fits{index_1(i)}, datarun.(datarun.default_sta_fits).sta_fits{index_1(ii)});
            dist_stixel(ii,i)=dist_stixel(i,ii);
            normdist(ii,i)=normdist(i,ii);
        end
    end
else
    for i=1:length(index_1)
        for ii=1:length(index_2)
            sta_fit_a=datarun.(datarun.default_sta_fits).sta_fits{index_1(i)};
            sta_fit_b=datarun.(datarun.default_sta_fits).sta_fits{index_2(ii)};
            if isempty(sta_fit_a)|isempty(sta_fit_b)
                warning(sprintf('sta_fit empty: %d %d',datarun.cell_ids(index_1(i)),datarun.cell_ids(index_2(ii)) ))
                dist=NaN; 
                normdist=NaN;
            else            
                [dist_stixel(i,ii) normdist(i,ii)] = sta_fit_distance(sta_fit_a,sta_fit_b);
            end
        end
    end
end


%mean sd
if ~exist('cell_specification_norm','var');
    if length(index_1)<10 || length(index_2)<10
        warning(sprintf('get_sta_fit_distance: mean sd computed on small numbers: %d,%d', length(index_1),length(index_2)))
    end
    
    temp=zeros(size(index_1));
    for i=1:length(index_1)
        temp(i)=sqrt(prod(datarun.(datarun.default_sta_fits).sta_fits{index_1(i)}.sd));
    end
    r1=robust_mean(temp);
    
    temp=zeros(size(index_2));
    for i=1:length(index_2)
        temp(i)=sqrt(prod(datarun.(datarun.default_sta_fits).sta_fits{index_2(i)}.sd));
    end
    r2=robust_mean(temp);
    
    r=mean([r1 r2]);
else
    ind=get_cell_indices(datarun,cell_specification_norm);
    temp=zeros(size(ind));
    for i=1:length(ind)
        temp(i)=sqrt(prod(datarun.(datarun.default_sta_fits).sta_fits{ind(i)}.sd));
    end
    r=robust_mean(temp);
end

   
dist_sd=dist_stixel/r;

















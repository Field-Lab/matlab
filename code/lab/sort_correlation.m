function pairs_sorted=sort_correlation(datarun,cell_specification_1,cell_specification_2,varargin)
% returns list of pairs in order of correlation
%
% synchrony_index=correlation(datarun, cell_specification_1, [cell_specification_2], [params]) 
%
% defaults.bin = 10; see correlation


% SET UP OPTIONAL ARGUMENTS
    p = inputParser;
    p.addParamValue('synchrony_index_field', 'synchrony_index');% 
    p.addParamValue('method', 'synchrony');%
  
    % parse inputs
    p.parse(varargin{:});
    params = p.Results;

        
% get cell numbers
index_1=get_cell_indices(datarun,cell_specification_1);
if ~exist('cell_specification_2','var');
    index_2=index_1;
    cell_specification_2=cell_specification_1;
else
    index_2=get_cell_indices(datarun,cell_specification_2);
end


if isequal(params.method,'synchrony')    
    %get synchrony_index        
    if isfield(datarun, params.synchrony_index_field)
        synchrony_index=full(datarun.(params.synchrony_index_field));
    else
        error(sprintf('sort_correlation: field %s does not exit',params.synchrony_index_field));
    end

    if isequal(index_1,index_2)
        syn=ones(length(index_1))*NaN;
        for i=1:length(index_1)
            for ii=i+1:length(index_1)   
                syn(i,ii)=synchrony_index(index_1(i),index_2(ii));
            end
        end
    else
        syn=synchrony_index(index_1,index_2);
    end 
    syn=full(syn);


    [junk t]=sort(syn(:));
    tt=find(~isnan(junk));
    tt=flipud(tt);
    [t1,t2]=ind2sub(size(syn),t(tt));

    pairs_sorted=[datarun.cell_ids(index_1(t1))' datarun.cell_ids(index_2(t2))' junk(tt)];

end


if isequal(params.method,'distance')    
    
    [dist_stixel dist_sd normdist] = get_sta_fit_distance(datarun, cell_specification_1, cell_specification_2);

    if isequal(index_1,index_2)
        dist=ones(length(index_1))*NaN;
        for i=1:length(index_1)
            for ii=i+1:length(index_1)   
                dist(i,ii)=normdist(i,ii);
            end
        end
    else
        dist=normdist;
    end 


    [junk t]=sort(dist(:));
    tt=find(~isnan(junk));
    [t1,t2]=ind2sub(size(dist),t(tt));

    pairs_sorted=[datarun.cell_ids(index_1(t1))' datarun.cell_ids(index_2(t2))' junk(tt)];

end





















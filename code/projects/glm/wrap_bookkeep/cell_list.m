% AKHeitman 2014-03-29
% Trying to bring better organization to our cell selection process
% Avoid annoying loop structures over and over again


function [exp_nm,cid_cell,expname,badcells]  = cell_list( expnumber, cellselectiontype)


%%% Identify exp_nm  %%%
if expnumber == 1, exp_nm ='2012-08-09-3'; expname = 'expA'; end
if expnumber == 2, exp_nm ='2012-09-27-3'; expname = 'expB'; end
if expnumber == 3, exp_nm ='2013-08-19-6'; expname = 'expC'; end 
if expnumber == 4, exp_nm ='2013-10-10-0'; expname = 'expD'; end 


%%% Find appropriate CID %%%
if strcmp(cellselectiontype, 'shortlist')    
    if expnumber == 1, cid_cell = {5161 841 5086 1471 1426 1772 2101 1276 3676}; end % 3676 and 1772 doesnt work for power raise on Bertha
    if expnumber == 2, cid_cell = {1 31 301 1201 1726 91 1909 2360 6858};  end
    if expnumber == 3, cid_cell = {2824 3167 3996 5660 6799 737 1328 1341 2959 5447}; end %  5447}; end 
    if expnumber == 4, cid_cell =  {32 768 2778 4354 5866 7036 346 1233 3137  5042 5418};end %cid_cell = { 3137  5042 5418}
end


if strcmp(cellselectiontype, 'debug')    
    if expnumber == 1, cid_cell = {1471 841}; end %{5086 841}
    if expnumber == 2, cid_cell = {1 6858};  end
    if expnumber == 3, cid_cell = {2824  5447}; end 
    if expnumber == 4, cid_cell = { 5866 5418};end 
end

if strcmp(cellselectiontype, 'all')
    if expnumber == 1, cid_cell = [];badcells =  [7756]; end % 3676 and 1772 doesnt work for power raise on Bertha
    if expnumber == 2, cid_cell = [];badcells =  [];  end
    if expnumber == 3, cid_cell = [];badcells =  []; end 
    if expnumber == 4, cid_cell = [];badcells =  [];end %cid_cell = { 3137  5042 5418}
end


end
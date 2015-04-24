function datarun = load_obvius_sta_fits(datarun, varargin)
% LOAD_OBVIUS_STA_FITS     Load obvius sta fits from the server
%
%
%
% usage:  datarun = load_obvius_sta_fits(datarun, params)
%
% arguments:  datarun - datarun struct with field specifying the path to the fits:
%                           datarun.names.obvius_fit_path
%
% outputs:    datarun - datarun struct with the following fields added, as possible
%
% 	piece information:
%     	datarun.obvius.stas_fits
%
% optional fields in params, their default values, and what they specify:
%
% sort_ids       false      when determining cell ID corresponding to an
%                           obvius STA fit id, sort cells recorded by a particular electrode
%                           according to id (this corrects for bug that
%                           occurs when cell ids in neurons or params file
%                           are not already sorted)
%
% 2009-04  greschner
% 2010-04  gauthier
% 2012-03  jepson - added sort_ids parameter
%



% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('sort_ids', false);

% resolve user input and default values
p.parse(varargin{:});

% get params struct
params = p.Results;



disp('load obvius STA fits')

if ~isempty(datarun.names.obvius_fit_path)
    sta_fit_full_path = datarun.names.obvius_fit_path
end


% access directory %%%%%%%%%%%%%%%%%%%%%
if ~isdir(sta_fit_full_path)
    error(['load_sta_:' sta_fit_full_path 'is no directory']);
else
    temp=dir([sta_fit_full_path '/stas*']);
    sta_fit_full_path = [sta_fit_full_path '/' temp.name '/'];
    if ~isdir(sta_fit_full_path)
        error(['load_sta_:' sta_fit_full_path 'is no directory']);
    end
end
files=dir([sta_fit_full_path 'sta-fit-*']);

% loop through files %%%%%%%%%%%%%%%%%%%%%
sta_fits = cell(length(datarun.cell_ids), 1);
for j=1:length(files)


    % open file %%%%%%%%%%%%%%%%%%%%%
    fid=fopen([sta_fit_full_path files(j).name]);
    if (fid == -1)
        error(['load_index: can not open ' full_path files(j).name]);
    end

    file=textscan(fid, '%s');
    fclose(fid);
    tex=file{1};

    i=0;
    while i<length(tex)
        i=i+1;
        if ~isempty(strfind(tex{i},':')); %scan until found new field name
            t=strrep(tex{i},':','');
            t=strrep(t,'(','');
            t=strrep(t,'-','_');
            t=lower(t);
            sta(j).(t)=[]; % make index entry new struct field name
        else
            tex{i}=strrep(tex{i},'(','');
            tex{i}=strrep(tex{i},')','');
            if ~isempty(str2num(tex{i}))
                sta(j).(t)=[sta(j).(t) str2num(tex{i})];
            else
                sta(j).(t)=[sta(j).(t) tex{i}];
            end
        end
    end
    sta(j).obvius_id=files(j).name;

    t = str2num(files(j).name(10:12));
    t = find(datarun.channels == t);
        
    if params.sort_ids
        [~,sortInd] = sort(datarun.cell_ids(t));
        t = datarun.cell_ids(t(sortInd));
    else 
        t = datarun.cell_ids(t);
    end

    if     files(j).name(13)=='A'; id=t(1);
    elseif files(j).name(13)=='B'; id=t(2);
    elseif files(j).name(13)=='C'; id=t(3);
    elseif files(j).name(13)=='D'; id=t(4);
    elseif files(j).name(13)=='E'; id=t(5);
    elseif files(j).name(13)=='F'; id=t(6);
    elseif files(j).name(13)=='G'; id=t(7);
    elseif files(j).name(13)=='H'; id=t(8);
    elseif files(j).name(13)=='I'; id=t(9);
    elseif files(j).name(13)=='J'; id=t(10);
    elseif files(j).name(13)=='K'; id=t(11);
    elseif files(j).name(13)=='L'; id=t(12);
    elseif files(j).name(13)=='M'; id=t(13);
    elseif files(j).name(13)=='N'; id=t(14);
    elseif files(j).name(13)=='O'; id=t(15);
    elseif files(j).name(13)=='P'; id=t(16);
    elseif files(j).name(13)=='Q'; id=t(17);
    elseif files(j).name(13)=='R'; id=t(18);
    elseif files(j).name(13)=='S'; id=t(19);
    elseif files(j).name(13)=='T'; id=t(20);
    elseif files(j).name(13)=='U'; id=t(21);
    elseif files(j).name(13)=='V'; id=t(22);
    elseif files(j).name(13)=='W'; id=t(23);
    elseif files(j).name(13)=='X'; id=t(24);
    elseif files(j).name(13)=='Y'; id=t(25);
    elseif files(j).name(13)=='Z'; id=t(26);
    else error(['Unrecognized letter at end of ' files(j).name]);
    end
        
    sta(j).cell_id = id;
    temp_index = get_cell_indices(datarun, id);
    sta_fits{temp_index} = sta(j);

end


datarun.obvius.sta_fits=sta_fits;













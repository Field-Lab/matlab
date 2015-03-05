function headerField=read_header_field_heka(pathway, fileName, field)
% read_header_field_heka reads heka file and extracts a field from the
% header (doesn't read data)
% INPUT
%   pathway: pathway to HEKA file
%   fileName: file name
%   field: string, name of the field you want to get. E.g. 'Stimulus Protocol'.
% OUTPUT
%   headerField: string with everything contained in this field.
% SPECIAL CASES:
%   if field is 'Stimulus Protocol': returns a matrix with n rows
%   and 5 columns (Time(ms),bmpID,posX(um),posY (um),Voltage(mV)).
%   Time in the current string corresponds to the time when this bmpID
%   finished.

fid=fopen(fullfile(pathway,fileName),'r');
if fid==-1
    fprintf('\nread_header_heka: CANNOT OPEN FILE\n%s\n',fullfile(pathway,fileName))
end
h=fread(fid,1,'int',0,'b');
header=fread(fid,h,'char=>char',0,'b');
fclose (fid);
header=header';

[~, headerField]=regexp(header,['\r',field,': ([^\r]*)\r'],'match','tokens'); % returns all characters of the stim protocol field (cell array of characters)
headerField=char(headerField{1,1}(1)); % converts to string

if strcmp(field,'Stimulus Protocol')
    [~, headerField]=regexp(headerField,'\n(.*)','match','tokens'); % throws away headline
    headerField=char(headerField{1,1}(1)); % converts to string again
    headerField=str2num(headerField); % converts to numbers (doubles)
end
function fitToJSON_sub(celldat,fit,fid,locs_s,subNum)

% create json object for a subunit

fprintf(fid,'"sub_id":%g,',subNum-1);
fprintf(fid,'"x":%g,',locs_s(subNum,1));
fprintf(fid,'"y":%g',locs_s(subNum,2));
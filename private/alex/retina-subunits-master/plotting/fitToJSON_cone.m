function fitToJSON_cone(celldat,fit,fid,cellNum,coneNum,opts)

% produce a json object for a cone 

fprintf(fid,'"cell_id":%g,',cellNum-1);
fprintf(fid,'"cone_id":%g,',celldat.coneIds(coneNum));
switch opts.fitType
  case {'greedy','greedy-local'}
    fprintf(fid,'"sub_id":%g,',find(fit.I_sc(:,coneNum))-1);
    fprintf(fid,'"sub_count":%g,',sum(fit.I_sc(find(fit.I_sc(:,coneNum)),:)));
    fprintf(fid,'"cone_weight":%g,',sum(fit.A_sc(:,coneNum)));
  case {'LN'}
    fprintf(fid,'"sub_id":%g,',0);
    fprintf(fid,'"sub_count":%g,',1);
    fprintf(fid,'"cone_weight":%g,',abs(fit.d_s(coneNum)));
end
fprintf(fid,'"x":%g,',celldat.locs_c(coneNum,1));
fprintf(fid,'"y":%g',celldat.locs_c(coneNum,2));


%subIds = find(sum(fit.I_sc,2)>1)-1;
%if sum((find(fit.I_sc(:,coneNum))-1)==subIds)
%  subCount = find((find(fit.I_sc(:,coneNum))-1)==subIds);
%else
%  subCount = 0;
%end
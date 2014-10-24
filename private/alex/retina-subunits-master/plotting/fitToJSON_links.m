function fitToJSON_links(celldat,fit,fid,cellNum,subNum)

fprintf(fid,'"source_cell":%g,',cellNum-1);
fprintf(fid,'"source_sub":%g,',subNum-1);
fprintf(fid,'"target_cell":%g,',cellNum-1);
fprintf(fid,'"value":%g',fit.B_s(subNum));
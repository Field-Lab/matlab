function fitToJSON_celltype(dat,fid,opts)

%%
fprintf(fid,'"cones":[');
%%
allinds = 1:dat.nCones;
for cellNum = 1:opts.nCells
  [celldat fit] = fitToJSON_load_cell(dat,opts,cellNum);
  allinds(celldat.keepInds) = 0;
  for coneNum=1:celldat.nConesFit
    fprintf(fid,'{');
    fitToJSON_cone(celldat,fit,fid,cellNum,coneNum,opts);
    if coneNum == celldat.nConesFit && cellNum == opts.nCells
      fprintf(fid,'}');
    else
      fprintf(fid,'},');
    end
  end
  
end
%%
fprintf(fid,'],');
%%
fprintf(fid,'"remCones":[');
%%
remcones = 1:dat.nCones;
for i=1:length(remcones)
  coneNum = remcones(i);
  fprintf(fid,'{');
  fprintf(fid,'"x":%g,',dat.locations(coneNum,1));
  fprintf(fid,'"y":%g',dat.locations(coneNum,2));
  if i == length(remcones)
    fprintf(fid,'}');
  else
    fprintf(fid,'},');
  end
end
%%
fprintf(fid,'],');
%%
fprintf(fid,'"cellCenters":[');
%%
for cellNum = 1:opts.nCells
  [celldat, fit] = fitToJSON_load_cell(dat,opts,cellNum);
  
  impath = sprintf('/Users/freemanj11/Dropbox/hacking/subunits/public/pngs/%s/',dat.dataSet);
  filepath = fullfile(impath,sprintf('%g.png',celldat.rgcId));
  if ~isdir(impath)
    mkdir(impath);
  end
  if ~isfile(filepath) % save out the image of the rf
    clf;
    image(norm_image(mean(celldat.rf,3)));
    hold on
    plot(celldat.locs_c(:,1)-celldat.xRange(1)+1,celldat.locs_c(:,2)-celldat.yRange(1)+1,'k.','MarkerSize',20);
    axis image off;
    g = getframe();
    %imwrite(matrix_scaled_up(norm_image(mean(celldat.rf,3)),20),filepath);
    imwrite(frame2im(g),filepath);
  end
  fprintf(fid,'{');
  fitToJSON_cell(celldat,fid,cellNum,fit);
  locs_s = normMat(fit.I_sc,'L1',2)*celldat.locs_c;
  fprintf(fid,',"subCenters":[');
  for subNum=1:fit.S
    fprintf(fid,'{');
    fitToJSON_sub(celldat,fit,fid,locs_s,subNum)
    if subNum == fit.S
      fprintf(fid,'}');
    else
      fprintf(fid,'},');
    end
  end
  fprintf(fid,']');
  if cellNum == opts.nCells
    fprintf(fid,'}');
  else
    fprintf(fid,'},');
  end
end
%%
fprintf(fid,'],');
%%
fprintf(fid,'"subCenters":[');
%%
for cellNum = 1:opts.nCells
  [celldat fit] = fitToJSON_load_cell(dat,opts,cellNum);
  locs_s = normMat(fit.I_sc,'L1',2)*celldat.locs_c;
  for subNum=1:fit.S
    fprintf(fid,'{');
    fitToJSON_sub(celldat,fit,fid,locs_s,subNum)
    if subNum == fit.S && cellNum == opts.nCells
      fprintf(fid,'}');
    else
      fprintf(fid,'},');
    end
  end
end
%%
fprintf(fid,'],');
%%
fprintf(fid,'"links":[');
%%
for cellNum = 1:opts.nCells
  [celldat fit] = fitToJSON_load_cell(dat,opts,cellNum);
  for subNum=1:fit.S
    fprintf(fid,'{');
    fitToJSON_links(celldat,fit,fid,cellNum,subNum);
    if subNum==fit.S && cellNum == opts.nCells
      fprintf(fid,'}');
    else
      fprintf(fid,'},');
    end
  end
end
%%
fprintf(fid,']');
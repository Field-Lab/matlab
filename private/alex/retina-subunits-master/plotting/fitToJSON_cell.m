function fitToJSON_cell(celldat,fid,cellNum,fit)

% produce a json object for a cell

center = getCentroid(celldat.locs_c);

fprintf(fid,'"cell_id":%g,',cellNum-1);
fprintf(fid,'"rgc_id":%g,',celldat.rgcId);
fprintf(fid,'"rf_x_px":%g,',size(celldat.rf,2));
fprintf(fid,'"rf_y_px":%g,',size(celldat.rf,1));
fprintf(fid,'"r2_LN":%g,',fit.r2_LN);

if isfield(fit,'r2_SUB')
  fprintf(fid,'"r2_SUB":%g,',fit.r2_SUB);
  fprintf(fid,'"r2_LN_SEL":%g,',fit.r2_LN_SEL);
  fprintf(fid,'"r2_SUB_SEL":%g,',fit.r2_SUB_SEL);
end

if isfield(fit,'f');
  f_in = fit.f.knots;
  f_out = evalNonLin(f_in,fit.f);
  fprintf(fid,'"f":[');
  for i=1:length(f_in)
    fprintf(fid,'{');
    fprintf(fid,'"x":%g,',f_in(i));
    fprintf(fid,'"y":%g',f_out(i));
    if i == length(f_in)
      fprintf(fid,'}');
    else
      fprintf(fid,'},');
    end
  end
  fprintf(fid,'],');
end

fprintf(fid,'"x":%g,',center(1));
fprintf(fid,'"y":%g',center(2));


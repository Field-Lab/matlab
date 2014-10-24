function fitToJson(filename)

fid = fopen(filename,'w');
%%
fprintf(fid,'{');
%%
fprintf(fid,'"retina":[');

%%
list = masterList();
dataSets = list.names;
dataSetsCells = [ones(1,length(list.midget.off)) ones(1,length(list.midget.on))*2];

for i=1:length(dataSets)
  checkData('laptop',dataSets{i});
end

%%
opts = [];

% get greedy subunit fits, off midget
for i=1:length(dataSets)
  opts_tmp.dataSet = dataSets{i};
  opts_tmp.percent = 0.33;
  opts_tmp.r2thresh = 0.05;
  if dataSetsCells(i) >= 1
    opts_tmp.cellType = 'off midget';
    opts_tmp.fitType = 'greedy';
    opts = cat(1,opts,opts_tmp);
    opts_tmp.fitType = 'LN';
    opts = cat(1,opts,opts_tmp);
  end
  if dataSetsCells(i) >= 2
    opts_tmp.cellType = 'on midget';
    opts_tmp.fitType = 'greedy';
    opts = cat(1,opts,opts_tmp);
    opts_tmp.fitType = 'LN';
    opts = cat(1,opts,opts_tmp);
  end
end


%%

dat = [];

for i=1:length(opts)
  
  if isempty(dat) || ~strcmp(dat.dataSet,opts(i).dataSet)
    clear dat
    fprintf('(fitToJson) loading data set %s\n',opts(i).dataSet);
    dat = loadData('laptop',opts(i).dataSet);
  end
  
  fprintf('(fitToJson) jsoning %s, %s, %s\n',...
    opts(i).dataSet,opts(i).cellType,opts(i).fitType);

  
  opts(i).cellInds = checkFits(dat,opts(i));
  opts(i).nCells = length(opts(i).cellInds);
  
  fprintf(fid,'{');
  fprintf(fid,'"cell_type":"%s",',replaceSpaces(opts(i).cellType));
  fprintf(fid,'"fit_type":"%s",',opts(i).fitType);
  fprintf(fid,'"data_set":"%s",',opts(i).dataSet);
  fprintf(fid,'"percent":"%g",',opts(i).percent);
  fprintf(fid,'"num_cells":"%g",',opts(i).nCells);
  fitToJSON_celltype(dat,fid,opts(i));
  
  if i == length(opts)
    fprintf(fid,'}');
  else
    fprintf(fid,'},');
  end
  
end

%%
fprintf(fid,']');
%%
fprintf(fid,'}');

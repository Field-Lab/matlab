% input: celldat
% output: celldat
% 
% add a field with the size of each gaussian fit to the rf

% input is the size of each gaussian
% error is the squared difference between real rf and gaussian bumps

function celldat = fitGaussCones(celldat)

%%
options = optimset('Display','off','Algorithm','active-set');
x0 = ones(1,celldat.nConesFit);
%%
szEst = fmincon(@(x) fitGaussCones_errFun(celldat,x),x0,...
  [],[],[],[],zeros(1,celldat.nConesFit)+0.01,ones(1,celldat.nConesFit)*3,...
  [],options);

celldat.rfFit.sz = szEst;
celldat.rfFit.amp = abs(celldat.sta(celldat.coneIds));
celldat = fitGaussCones_fitFun(celldat);

function err = fitGaussCones_errFun(celldat,p)
celldat.rfFit.sz = p(1:celldat.nConesFit);
celldat.rfFit.amp = abs(celldat.sta(celldat.coneIds));
celldat = fitGaussCones_fitFun(celldat);
err = sum(vector((celldat.rfFit.fit/max(celldat.rfFit.fit(:))-celldat.rf/max(celldat.rf(:))).^2));


function celldat = fitGaussCones_fitFun(celldat)

dims = size(celldat.rf);
[xramp,yramp] = meshgrid([1:dims(1)],[1:dims(2)]);

fit_rf = zeros(size(celldat.rf));
for icone=1:celldat.nConesFit
  center = [celldat.locs_c(icone,1)-min(celldat.xRange)+1 ...
    celldat.locs_c(icone,2)-min(celldat.yRange)+1]; 
  xramp_tmp = xramp - center(2);
  yramp_tmp = yramp - center(1);
  ampl = 1/(2*pi*celldat.rfFit.sz(icone));
  e = (xramp_tmp.^2 + yramp_tmp.^2)/(-2 * celldat.rfFit.sz(icone));
  res = ampl .* exp(e);
  gauss = res';
  gauss = gauss/max(gauss(:));
  fit_rf = fit_rf + gauss*celldat.rfFit.amp(icone);
end
celldat.rfFit.fit = fit_rf;
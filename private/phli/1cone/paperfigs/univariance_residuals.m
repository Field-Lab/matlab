% 2013-05-04 PHLI
% First run cr_compound_plot and capture all outputs, then run below
% to collect statistics for GDoc:
% https://docs.google.com/spreadsheet/ccc?key=0Am2KjZWyH8AndDZ4M2hGbERyUy01dkRTaDJKd0lBbEE#gid=0

id = gcf;
cones = [1 3 4]; 
nozero = [1:2 4:7]; 
pols = 4:7;

rgc = find(cellfun(@(C) isequal(C, id), rasterrun.rgcs));
stdres =    std(reshape(residuals{rgc}(cones, nozero), [], 1));
stdrespol = std(reshape(residuals{rgc}(cones, pols),   [], 1));
maxresp = max(reshape(crs(:,:,rgc), [], 1));

disp(sprintf('%d\t%f\t%f\t%f', id, stdres, stdrespol, maxresp));


%% Results brought back from GDoc
% Previous values before Neuron revision
% offM = [0.0300 0.0213 0.0330 0.0208 0.0244 0.0314 0.0158 0.0238 0.0234 0.0501 0.0232 0.0178 0.0113 0.0106 0.0164 0.0205 0.0162 0.0158 0.0186 0.0273 0.0092 0.0287];
% offP = [0.0131 0.0194 0.0288 0.0210 0.0231 0.0201 0.0244 0.0293 0.0278 0.0219 0.0245 0.0192 0.0301];
% onM  = [0.0699 0.0650 0.0525 0.0291 0.0321 0.0150 0.0163 0.0238 0.0095 0.0090 0.0161 0.0318 0.0001];

% Recalculated and checked
offM = [0.0205    0.0213    0.0330    0.0208    0.0165    0.0244    0.0158    0.0238    0.0234    0.0501    0.0232    0.0178    0.0113    0.0106    0.0164    0.0205    0.0162    0.0158    0.0186    0.0273    0.0092    0.0287];
offP = [0.0131    0.0194    0.0288    0.0210    0.0231    0.0201    0.0244    0.0293    0.0278    0.0219    0.0245    0.0192    0.0301];
onM  = [0.0097    0.0650    0.0557    0.0291    0.0277    0.0122    0.0179    0.0238    0.0095    0.0090    0.0161    0.0318];

cla
hold on;
plot(onM, zeros(size(onM)), 'ko')
plot(offM, 0.0015.*ones(size(offM)), 'ks')
plot(offP, 0.003.*ones(size(offP)), 'k^')
axis equal
set(gca, 'XLim', [-0.004 0.08], 'XTick', 0:0.02:0.085, 'YLim', [-0.0015 0.0045], 'YTick', [], 'YColor', [1 1 1])

% Figure cells
plot([0.0097 0.0106 0.0158 0.0244], [0 0.0015 0.0015 0.003], 'rx');
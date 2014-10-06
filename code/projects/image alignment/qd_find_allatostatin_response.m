

% load datarun
if ~exist('datarun','var')
    switch 4
        case 1
            drun_spec = '2009-06-03-1/data004/data004';
            act = 2; % allatostatin-sensitive cell type
            timepoints = [880 2640 3340];
        case 2
            drun_spec = '2009-07-23-0/data000/data000';
            act = 2; % allatostatin-sensitive cell type
            timepoints = [3600 5630 8400];
        case 3
            drun_spec = '2009-07-23-1/data002/data002';
            act = 2; % allatostatin-sensitive cell type
            timepoints = [4520 6400 9000];
        case 4
            drun_spec = '/Volumes/Creampuff/Analysis/Gauthier/2009-07-27-0/data000/data000';
            act = 1; % allatostatin-sensitive cell type
            timepoints = [3880 7340 10090 12410];
    end
    datarun = load_data(drun_spec);
    datarun = load_neurons(datarun);
    datarun = load_params(datarun,'order_cell_types',0,'verbose',1);
end


% parameters
edges = [0:50:datarun.duration];


% get expected spike rate vector
expected = histc(vertcat(datarun.spikes{get_cell_indices(datarun,datarun.cell_types{act}.cell_ids)}),edges);
% plot it
figure(1);clf;plot(edges,expected);hold on
for tt=1:length(timepoints)
    plot([1 1]*timepoints(tt),get(gca,'ylim'),'r');end
title('mean spike rate of cells with apparent allatostatin response')

% compute projections
cell_indices = get_cell_indices(datarun,'all');

% initialize
proj = zeros(length(cell_indices),1);
srv = cell(length(cell_indices),1);

% for each cell
for cc=1:length(cell_indices)
    % get spike rate vector
    srv{cc} = histc(datarun.spikes{cell_indices(cc)},edges);
    % project
    proj(cc) = sum(srv{cc}.*expected)/norm(srv{cc});
end

% plot projection of each cell
figure(2);clf;hist(proj,40)

% plot spike rate of cells, with bars showing when allatostatin was added
% first the cells in the allatostatin class, then others, sorted by similarity to the allatostatin class spike rate
figure(3);clf;hold off;set(3,'color','w')
% get allatostatin list
cell_indices_allatostatin = get_cell_indices(datarun,{act});
% get sorted list
[junk,srtd]=sort(proj);
% concatenate lists
ci_to_plot = [cell_indices_allatostatin fliplr(srtd')];

% plot timecourses
ha=loop_slider_list(1,1,length(ci_to_plot));

while 1
    pp = ci_to_plot(round(get(ha,'Value')));
    cla;
    plot(edges,srv{pp});hold on
    for tt=1:length(timepoints)
        plot([1 1]*timepoints(tt),get(gca,'ylim'),'r');end
    title(sprintf('cell id %d',datarun.cell_ids(cell_indices(pp))))
    uiwait
end

%    list=[]; list2=[]; 
%    [ha list list2]=loop_slider_list(k,kmin,kmax,list,list2);
%    while k
%       k=round(get(ha,'Value'));
%       plot(mat(:,k));
%       uiwait;
%    end


% for pp=[cell_indices_allatostatin fliplr(srtd')]
%     clf;
%     plot(edges,srv{pp});hold on
%     for tt=1:length(timepoints)
%         plot([1 1]*timepoints(tt),get(gca,'ylim'),'r');end
%     title(sprintf('cell id %d',datarun.cell_ids(cell_indices(pp))))
%     pause
% end
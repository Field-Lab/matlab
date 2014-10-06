% load axons paths from an excel spreadsheet
% NOTE: NOW OBSELETE, COPIED TO "qd_rgc_alignment.m"


% read spreadsheet
sp_path = '/snle/lab/Experiments/Array/Analysis/2007-09-18-4/images/2007-09-18-4-axons.xls';
[junk, junk, axons_xls] = xlsread(sp_path);

% load coordinates of the path for each axon
num_axons = size(axons_xls,2)/2;
axons = cell(num_axons,1);
for aa=1:num_axons
    start_col = aa*2-1;
    % determine if contains numeric data
    if isnumeric(axons_xls{2,start_col})
        % load entries
        axons{aa}(:,1) = cell2mat({axons_xls{2:end,start_col}});
        axons{aa}(:,2) = cell2mat({axons_xls{2:end,start_col+1}});
        % remove the trailing NaN entries, which correspond to empty fields
        axons{aa} = axons{aa}(~isnan(axons{aa}(:,1)),:);
    end
end

% transform to array space
axons_orig = axons;
% make transformation
TA1F1_reverse = cp2tform(base_points,input_points,'lwm'); % these should be from the TA1F1 transformation
T=maketform('composite',TA1F1_reverse,fliptform(TA1));
% apply to each axon
for aa=1:num_axons
    if ~isempty(axons{aa})
        axons{aa} = tforminv(T,axons_orig{aa});
    end
end


% plot, for fun
if 1
    figure(1);clf;imagesc(f1);axis image; colormap gray; hold on; aa(1)=gca;
    for aa=1:num_axons
        if ~isempty(axons_orig{aa})
            plot(axons_orig{aa}(:,1),axons_orig{aa}(:,2),'Color',(rand(3,1)+0.2)/1.2,'LineWidth',0.2)
            text(axons_orig{aa}(1,1),axons_orig{aa}(1,2),num2str(aa),'Color',[0.5 0 0],'FontSize',1,...
                'HorizontalAlignment','Center','VerticalAlignment','Bottom')
        end
    end
end


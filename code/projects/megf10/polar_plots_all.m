function[] = polar_plots_all(U,V,num, mag)

%Input X and Y vectors of all cells for all stimuli for all speeds
%Plots polar plot for all speeds for all cells

% Sneha Ravi 
% Last revision: 12-18-2012

axes_handle = [];
ylim = zeros(1, length(U));
X = floor(sqrt(length(U)));
for i = 1:length(U)
    ax(i) = subplot(X,ceil(length(U)/X),i);
    set(ax(i), 'FontName', 'AvantGarde', 'FontSize', 18)
    h1 = compass(U{i,1},V{i,1});
    h = findall(gca, 'type', 'line'); % find all lines in the current figure
    h(h==ax(i)) = []; % remove the line you ploted from the list.
    set(h, 'LineStyle', '--');
    set(h1,'linestyle','-', 'LineWidth', 1); 
    axes_handle = [axes_handle ax(i)];
    ylim(i) = max(mag{i,1});
    titlechar = [num2str(num(i))];
    %titlechar = ['Vector Averages Plot for all cells of temporal period: ' num2str(num(i))];
    %title(titlechar, 'Color', 'k', 'FontWeight', 'bold' , 'FontSize', 18 ,'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Right');
    title(titlechar);
end

%linkaxes(axes_handle,'xy'); %Link axes to have same x and y limits
%[C, I] = max(ylim);
%set(ax(I),'xlimmode','auto');
%set(ax(I),'ylimmode','auto');

end


%Return vector averages for each speed for each cell to next function

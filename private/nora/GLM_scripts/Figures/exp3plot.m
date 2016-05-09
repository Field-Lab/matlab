function exp3plot(bps_exp, dist)
hold on
default_colors = get(gca,'ColorOrder');
ExpColors = [default_colors(2,:); default_colors(5,:); default_colors(1,:)];
TypeSymbols = {'+','o'};
TypeLines = {'--', ''};
% scatter
if iscell(dist)
    for i = 1:3
        for j = 1:2
            plot(bps_exp{i, j}, dist{i, j}, TypeSymbols{j}, 'Color', ExpColors(i,:))
        end
    end
    % plot distribution
elseif ~isvector(bps_exp{1,1})
    hold on
    for i =3:-1:1
        for j = 1:2
            %subplot(2,1,j); hold on;
            %plot(bps_exp{i,j}', 'Color', [ExpColors(i,:), 0.1])
            hold on;
            avg = mean(bps_exp{i,j});
            %plot(avg, 'Color', ExpColors(i,:), 'LineWidth', 2.5)
            plot(avg, TypeLines{j}, 'Color', ExpColors(i,:), 'LineWidth', j)
        end
    end
elseif dist
    for i = 1:3
        for j = 1:2
            plot((i+(j/2-0.5))*ones(length(bps_exp{i,j}),1),bps_exp{i,j},TypeSymbols{j},'Color', ExpColors(i,:))
        end
    end
    xlim([0.5 4])
    set(gca, 'XTick', 1:0.5:3.5)
    set(gca, 'XTickLabels', {'On', 'Off'})
    
    % plot lines
else
    for i = 1:3
        for j = 1:2
            plot(bps_exp{i,j},TypeLines{j}, 'Color', ExpColors(i,:), 'LineWidth', j)
        end
    end
end
title(inputname(1))
end
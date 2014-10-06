function [] = print_cell_plots(M, N, P)

	x0 = M(:, 2);
	y0 = M(:, 3);
	a = M(:, 4);
	b = M(:, 5);
	rad = M(:, 6);
	alfa = linspace(0,2 * pi, 100);
    K = [ 'y' 'm' 'c' 'r' 'g' 'b' ];
    %imdata = imread('okno.png');
    time = 1;
    axis(P.Axis);
    
    figure(time)
    %imagesc(P.XBackgroundLimits, P.YBackgroundLimits, imdata);
    hold on;
                
    for i = 1:size(N, 1) % czasy
		if (i ~= 1)
			if (N(i, 1) ~= N(i - 1, 1))
                close(f);
				time = time + 1;
                
                figure(time)
                %imagesc(P.XBackgroundLimits, P.YBackgroundLimits, imdata);
                hold on;
        
            end
        end

        f = gcf;
        set(f, 'Visible', 'off');
        g = gca;
        set(g, 'Visible', 'off');
        
		for j = 1:size(M, 1) % kola
			x = x0(j) + a(j) * cos(alfa);
			y = y0(j) + b(j) * sin(alfa);
				if (N(i, 2) == j)
					h = fill(x, y, K(N(i, 3)), 'EdgeColor', P.EdgeColor, 'LineWidth', P.LineWidth);
					hold on;
				else
					h = plot(x, y, 'LineWidth', P.LineWidth);
					hold on;
                end
			rotate(h, [0 0 1], rad(j), [M(j, 2), M(j, 3), 0]);
            alpha(0.7);
        end
		
        axis(P.Axis);
        set(g, 'ydir', 'normal');
		set(f, 'PaperUnits', 'inches');
		set(f, 'PaperSize', P.PaperSize);
		set(f,'PaperPosition', P.PaperPosition); % xLeft yTop xSize ySize
		name = [P.PathName P.FileName int2str(time)];
		print(f, P.FormatType, P.Resolution, name);		
    end
    close(f);
end

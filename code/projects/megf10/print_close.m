function [] = print_close(figure_n, paper_size, name)

set(figure(figure_n), 'paperpositionmode', 'auto');
set(gcf, 'PaperUnits', 'inch');
set(figure(figure_n), 'PaperSize', paper_size);
print(figure(figure_n), '-dpdf', name)
close
end

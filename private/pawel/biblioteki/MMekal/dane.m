%1) Wstaw macierz M, macierz N oraz ustaw parametry P.
%2) Wywolaj i przypisz pod macierz funkcje: Nc = convert_matrix(N)
 %  w celu konwersji macierzy N na trojkolumnowa Nc.
%3) Wywolaj funkcje: print_cell_plots(M, Nc, P)

M = [
	1 2 2 3 5 30
	2 4 4 1 6 45
	3 2 2 3 1 0
	4 10 10 3 4 80
	5 13 10 3 3 0
	6 1 1 2 1 40
	7 13 7 3 4 0
	8 0 0 3 4 0
	9 6 6 3 1 35
	10 15 15 5 2 80
	] % nr_kola x0 y0 a b kat

N = [
	0 2 2 1
	0 3 2 3
	0 6 1 4
	0 10 6 1
	1.1 9 4 5
	1.1 10 4 4
	1.5 1 6 1
	1.6 1 5 1
	1.6 2 5 1
	1.6 3 5 7
	1.6 5 3 5
	1.8 1 3 3
	1.9 10 2 2
	2.3 4 1 1
	2.3 7 1 2
	2.3 8 2 3
	2.3 9 3 4
	2.3 10 4 3
	3.6 4 5 2
	3.6 6 6 1 
	3.9 6 1 2
	3.9 7 2 2
	4.5 8 3 1
	4.5 10 4 1
	] % czas nr_kola kolor czas_zycia

cd C:\home\Pawel\nauka\praktyki2011;
P = struct('PaperSize', [8 6], 'PaperPosition', [0 0 8 6], 'PathName', '', 'FileName', {'wykres'}, 'EdgeColor', 'k', 'LineWidth', 0.5, 'Resolution', '-r100', 'FormatType', '-dtiff', 'Axis', [-5 20 -5 20], 'XBackgroundLimits', [-10 30], 'YBackgroundLimits', [-10 30])
Nc = convert_matrix(N)
print_cell_plots(M, Nc, P);
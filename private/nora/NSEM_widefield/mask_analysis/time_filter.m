function time_filter = time_filter(b)
% 1 fit scale 1
% 2 tau one
% 3 n one
% 4 fit scale 2
% 5 tau two
% 6 n two
    time_filter = b(1).* ((1:30)./b(2)).^b(3) .* exp(-b(3)*(((1:30)./b(2)) - 1))-b(4).* ((1:30)./b(5)).^b(6) .* exp(-b(6)*(((1:30)./b(5)) - 1)); 
end
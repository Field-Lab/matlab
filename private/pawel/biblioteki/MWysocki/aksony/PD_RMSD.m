function rms = PD_RMSD( x0, y0, line)
    y = polyval(line,x0);
    rms = sqrt(sum((y0-y).^2)/length(y0));
end


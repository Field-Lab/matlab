function result = isSomaSpike(timecourse)
    [minV, index] = min(timecourse);
    maxV = max(timecourse(1:index));
    result = 1;
    if abs(maxV/minV) > 0.25
        result = 0;
    end
end
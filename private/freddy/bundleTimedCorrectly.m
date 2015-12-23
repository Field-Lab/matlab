function correct = bundleTimedCorrectly(bundle)

timingInfo = bundle(:, 3);

decreasing = true;

correct = true;

for i = 2:size(timingInfo, 1)
    curr = timingInfo(i);
    prev = timingInfo(i-1);
    if decreasing
        if curr - prev > 2
            decreasing = false;
            continue;
        end
    else
        if curr - prev < -2
            correct = false;
            break;
        end
    end
end
    

end
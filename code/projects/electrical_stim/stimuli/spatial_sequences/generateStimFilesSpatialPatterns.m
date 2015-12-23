function generateStimFilesSpatialPatterns(electrodes, fullAmps, threshAmps, offsets, nOrders, periodInMs)

nElec = length(electrodes);
%nOffsets = length(offsets);

if max(offsets)*(nElec+1) >= periodInMs;
    h = msgbox('need a longer period for this number of electrodes and length of time offsets');
    boxPos = get(h, 'position');
    set(h, 'position', [400 500 boxPos(3) boxPos(4)])
    return
end

%% generates distinct random orders of pulses

orders = zeros(nOrders, nElec);

if nOrders <= factorial(nElec)
    for ii = 1:nOrders
        newOrder = false;
        while ~newOrder
            randOrder = randperm(nElec);
            newOrder = true;
            for jj = 1:ii-1
                if randOrder == orders(jj,:)
                    newOrder = false;
                    break
                end
            end
        end
        orders(ii,:) = randOrder;
    end
else
    errordlg(['can''t generate ' num2str(nOrders) ' distinct orders of ' num2str(nElec) ' electrodes'])
end


%% varying all electrodes together near 100% response

Array = generatePatternSpatialPatternsVaryAllElecs(electrodes, fullAmps);

MovieChunksFile = generateMovieSpatialPatternsVaryAllElecs(length(electrodes), periodInMs, offsets, orders);

keyboard

fid = fopen('spatial_short_electrodes','wb','ieee-le.l64')
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('spatial_short_patterns','wb','ieee-le.l64')
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('spatial_short_movie','wb','ieee-le.l64')
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);



%% varying each electrode around threshold
    
[Array variedElec] = generatePatternSpatialPatternsVary1Elec(electrodes, fullAmps, threshAmps);

MovieChunksFile = generateMovieSpatialPatternsVary1Elec(length(electrodes), periodInMs, offsets, orders, variedElec);

keyboard

fid = fopen('spatial_vary_thresh_electrodes','wb','ieee-le.l64');
fwrite(fid,electrodes,'int32');
fclose(fid);

fid = fopen('spatial_vary_thresh_patterns','wb','ieee-le.l64');
fwrite(fid,Array,'double');
fclose(fid);

fid = fopen('spatial_vary_thresh_movie','wb','ieee-le.l64');
fwrite(fid,MovieChunksFile,'int32');
fclose(fid);
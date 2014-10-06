function electrodeOrder = electrodeOrderGenerator(plotOn)

% Generates a pseudorandom order of electrodes on the 61 MEAs that never has sequential electrodes
% within 90 microns (for the 30 micron array) or 180 microns (for the 60 micron array) of eachother.
%
% arguments: use 1 if the step-through plot is desired, or 0 if not
%
% returns: a vector containing the electrode numbers in the pseudorandom order (doesn't include
% disconnected electrodes)
%
% Algorithm:
%   9 "frames" of electrodes are defined
%       each frame is comprised of the triangular lattice of electrodes with 90 or 180 micron edges
%
%   randperm randomizes the order of the frames to be used (all electrodes within a frame are added to
%       electrodeOrder before the next frame is used)
%
%   within each frame, randperm randomizes the order of the electrodes within the frame
%
%   when finishing with one frame and going to the next, the distance between the last electrode in
%       the previous frame and the first electrode in the current frame is determined
%
%   if the distance is less than or equal to 90 or 180 microns, a new random order of electrodes
%       within the current frame is chosen




%% info about electrode coordinates (xCoords, yCoords)
% Originally defined by Anastacia
% They were modified by Justin so that down was towards electrode 30 and groundpin.
% Hexagonal array.
% 64 electrodes.
% Electrodes 9, 25, and 57 are disconnected.
% Central electrode is #41.
% The position of electrode #41 is taken to be the origin (0,0).
% X+ axis: from electrode 41 towards electrode 13.
% Y+ axis: from electrode 41 towards electrode 62.
% If the array is positioned ground pin down, then X is left-to-right and Y is down-to-up.
% The units are such that electrode spacing = 2.
% for a 60 micron spaced array, scale factor should be 30 e.g. 30*2 = 60 microns
% Justin 2008-04-30

%electrodeNumbers = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 ...
%    30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 ...
%    61 62 63 64];


%%

if ~(plotOn==1||plotOn==0)
    error('The argument must be either 0 (to suppress plot) or 1 (to allow plot).')
end

xCoords = [1.7320508 3.4641016 3.4641016 1.7320508 5.196152 3.4641016 5.196152 6.928203 nan ...
    6.928203 5.196152 3.4641016 6.928203 5.196152 1.7320508 6.928203 3.4641016 5.196152 6.928203 ...
    5.196152 3.4641016 1.7320508 3.4641016 1.7320508 nan 0 1.7320508 0 0 0 -1.7320508 -1.7320508 ...
    -1.7320508 -3.4641016 -3.4641016 -1.7320508 -5.196152 -3.4641016 -5.196152 -6.928203 0 ...
    -6.928203 -5.196152 -3.4641016 -6.928203 -5.196152 -1.7320508 -6.928203 -3.4641016 -5.196152 ...
    -6.928203 -5.196152 -3.4641016 -1.7320508 -3.4641016 -1.7320508 nan 0 -1.7320508 0 0 0 ...
    1.7320508 1.7320508];

yCoords = [3 6 4 1 5 2 3 4 nan 2 1 0 0 -1 -1 -2 -2 -3 -4 -5 -4 -3 -6 -5 nan -2 -7 -4 -6 -8 -7 -5 ...
    -3 -6 -4 -1 -5 -2 -3 -4 0 -2 -1 0 0 1 1 2 2 3 4 5 4 3 6 5 nan 2 7 4 6 8 7 5];

electrodeOrder = zeros(1, 61);

frame1Elecs = [8 16 63 4 24 53 38];
frame2Elecs = [10 19 64 15 27 49 35];
frame3Elecs = [13 1 22 55 44 34];
frame4Elecs = [5 14 62 58 28 52 43];
frame5Elecs = [7 18 61 41 29 50 39];
frame6Elecs = [11 20 60 26 30 46 37];
frame7Elecs = [2 12 23 54 33 45];
frame8Elecs = [3 17 59 47 32 51 42];
frame9Elecs = [6 21 56 36 31 48 40];


frameOrder = randperm(9);
lastFilled = 0;

for i = 1:9
    switch frameOrder(i)
        case 1
            orderInFrame = randperm(length(frame1Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame1Elecs(orderInFrame(1))) yCoords(frame1Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame1Elecs));
                    nextPosition = [xCoords(frame1Elecs(orderInFrame(1))) yCoords(frame1Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame1Elecs)
                electrodeOrder(lastFilled+j) = frame1Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame1Elecs);
        case 2
            orderInFrame = randperm(length(frame2Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame2Elecs(orderInFrame(1))) yCoords(frame2Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame2Elecs));
                    nextPosition = [xCoords(frame2Elecs(orderInFrame(1))) yCoords(frame2Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame2Elecs)
                electrodeOrder(lastFilled+j) = frame2Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame2Elecs);
        case 3
            orderInFrame = randperm(length(frame3Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame3Elecs(orderInFrame(1))) yCoords(frame3Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame3Elecs));
                    nextPosition = [xCoords(frame3Elecs(orderInFrame(1))) yCoords(frame3Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame3Elecs)
                electrodeOrder(lastFilled+j) = frame3Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame3Elecs);
        case 4
            orderInFrame = randperm(length(frame4Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame4Elecs(orderInFrame(1))) yCoords(frame4Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame4Elecs));
                    nextPosition = [xCoords(frame4Elecs(orderInFrame(1))) yCoords(frame4Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame4Elecs)
                electrodeOrder(lastFilled+j) = frame4Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame4Elecs);
        case 5
            orderInFrame = randperm(length(frame5Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame5Elecs(orderInFrame(1))) yCoords(frame5Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame5Elecs));
                    nextPosition = [xCoords(frame5Elecs(orderInFrame(1))) yCoords(frame5Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame5Elecs)
                electrodeOrder(lastFilled+j) = frame5Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame5Elecs);
        case 6
            orderInFrame = randperm(length(frame6Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame6Elecs(orderInFrame(1))) yCoords(frame6Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame6Elecs));
                    nextPosition = [xCoords(frame6Elecs(orderInFrame(1))) yCoords(frame6Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame6Elecs)
                electrodeOrder(lastFilled+j) = frame6Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame6Elecs);
        case 7
            orderInFrame = randperm(length(frame7Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame7Elecs(orderInFrame(1))) yCoords(frame7Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame7Elecs));
                    nextPosition = [xCoords(frame7Elecs(orderInFrame(1))) yCoords(frame7Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame7Elecs)
                electrodeOrder(lastFilled+j) = frame7Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame7Elecs);
        case 8
            orderInFrame = randperm(length(frame8Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame8Elecs(orderInFrame(1))) yCoords(frame8Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame8Elecs));
                    nextPosition = [xCoords(frame8Elecs(orderInFrame(1))) yCoords(frame8Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame8Elecs)
                electrodeOrder(lastFilled+j) = frame8Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame8Elecs);
        case 9
            orderInFrame = randperm(length(frame9Elecs));
            if i ~= 1
                prevPosition = [xCoords(electrodeOrder(lastFilled)) yCoords(electrodeOrder(lastFilled))];
                nextPosition = [xCoords(frame9Elecs(orderInFrame(1))) yCoords(frame9Elecs(orderInFrame(1)))];
                while norm(prevPosition - nextPosition) < 6
                    orderInFrame = randperm(length(frame9Elecs));
                    nextPosition = [xCoords(frame9Elecs(orderInFrame(1))) yCoords(frame9Elecs(orderInFrame(1)))];
                end
            end
            for j = 1:length(frame9Elecs)
                electrodeOrder(lastFilled+j) = frame9Elecs(orderInFrame(j));
            end
            lastFilled = lastFilled + length(frame9Elecs);
    end
end



%% displays electrode order
% a red marker, normal size, indicates the previous electrode in the order
% a red marker, larger size, indicates the current electrode (the one just added to the figure)
% markers are added one at a time: hit any key to add the next marker

if plotOn
    figure
    axis([-11 11 -11 11])
    text(0, 9.5,'62','HorizontalAlignment','center','VerticalAlignment','middle')
    text(0, -9.5,'30','HorizontalAlignment','center','VerticalAlignment','middle')
    text(9,0,'13','HorizontalAlignment','center','VerticalAlignment','middle')
    text(-9,0,'45','HorizontalAlignment','center','VerticalAlignment','middle')
    hold on
    disp('Press any key to advance the plot.')
    for i = 1:61
        if i>1
            plot(xCoords(electrodeOrder(i-1)), yCoords(electrodeOrder(i-1)),'c.', 'MarkerSize',20)
            if i>2
                plot(xCoords(electrodeOrder(i-2)), yCoords(electrodeOrder(i-2)),'k.', 'MarkerSize',20)
            end
        end
        plot(xCoords(electrodeOrder(i)), yCoords(electrodeOrder(i)),'r.', 'MarkerSize',20)
        pause
    end
    hold off  
end



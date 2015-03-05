
clear
[w, screenRect] = Screen(0,'OpenWindow',[120 120 120],[100 100 1400 900]); % open window, get screen size
Screen('Preference', 'OverrideMultimediaEngine',1)
[movie a b c d]=Screen('OpenMovie', w, '/Users/alexth/Desktop/test.avi');


Screen('PlayMovie', movie, 1,[],0);
pause(5)
tic
i=1;
while ~KbCheck
    % Wait for next movie frame, retrieve texture handle to it
    tex = Screen('GetMovieImage', w, movie);

    % Valid texture returned? A negative value means end of movie reached:
    if tex<=0
        % We're done, break out of loop:
        break;
    end;

    % Draw the new texture immediately to screen:
    Screen('DrawTexture', w, tex);

    % Update display:
    Screen('Flip', w);
    timeindex(i) = Screen('GetMovieTimeIndex', movie);
    i=i+1;

    % Release texture:
    Screen('Close', tex);
end;
toc
[droppedframes]=Screen('PlayMovie', movie, 0);
Screen('CloseMovie', movie);

Screen('CloseAll');
plot(diff(timeindex))

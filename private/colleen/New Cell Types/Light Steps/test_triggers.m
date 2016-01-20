trigger_marks = [];
counter = 1;
while counter< length(x)-1
    trigger_marks = [trigger_marks x(counter)];
    while x(counter+1) == x(counter) + 1;
        counter = counter+1;
        if counter> length(x)-1
            break;
        end
    end
    counter = counter+1;
end

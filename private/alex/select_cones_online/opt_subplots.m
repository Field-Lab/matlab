function [rows cols]=opt_subplots(amount)
if amount<4
    rows=1;cols=amount;
elseif amount==4
    rows=2;cols=2;    
elseif amount<7
    rows=2;cols=3;    
elseif amount<=8
    rows=2;cols=4;    
elseif amount==9
    rows=3;cols=3;    
elseif amount==10
    rows=2;cols=5;    
elseif amount<13
    rows=3;cols=4;
elseif amount<16
    rows=3;cols=5;
elseif amount==16
    rows=4;cols=4;
elseif amount<19
    rows=3;cols=6;
else
    rows=floor(sqrt(amount));
    if amount==rows*rows
        cols=rows;
    else
        if rows*(rows+1)<amount
            rows=rows+1;
            cols=rows;
        else
            cols=rows+1;
        end
    end
end
if amount==48
    rows=6;
    cols=8;
end

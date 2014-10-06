Xsize=100;
Ysize=100;
Xblock=80;
Yblock=80;
Xshift=10;
Yshift=10;
%Xstep=5;
%YStep=5;

signal=rand(Xsize,Ysize);

signal_to_move=signal(Xshift+1:Xsize-Xshift,Yshift+1:Ysize-Yshift);
w=zeros(2*Xshift-1,2*Yshift-1);
for xi=1%:2*Xshift
    xi
    for yi=1%:2*Yshift
        s=signal(xi:xi+Xblock-1,yi:yi+Yblock-1);
        w(xi,yi)=sum(sum(s*signal_to_move));
    end
end
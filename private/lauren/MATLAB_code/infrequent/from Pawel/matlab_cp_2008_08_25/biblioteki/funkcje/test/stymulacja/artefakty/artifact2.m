function w=artifact2(model,Vdc,t0,t1,A1,t2,A2,t3,A3,t_delay,t_disch,T,dt);
%Impuls pradowy trzyfazowy plus rozladowanie.
%model - model elektrody
%t0 - moment zaczecia impuslu (wzgledem poczatku symulacji)
%t1 - czas trwania pierwszej czesci impulsu
%A1 - amplituda pierwszej czesci
%t2 - czas trwanai drugiej czesci
%A2 - amplituda drugiej czesci
%T - czas symulacji
%dt - krok czasowy

Ce=model.Ce;
beta=model.beta;
Re=model.Re;
Rs=model.Rs;
Idc=Vdc/Re;

N=floor(T/dt);
n0=floor(t0/dt);
n1=floor(t1/dt);
n2=floor(t2/dt);
n3=floor(t3/dt);
n_delay=floor(t_delay/dt);
n_disch=floor(t_disch/dt);

Ucpe=zeros(1,N);
Ire=zeros(1,N);
Icpe=Ire;
Irs=Ire;
U=Ire;
%gamma=1.5;
czas=[1:N]*dt;
tran_cpe=1./(czas.^(1-beta));
tran_cpe=tran_cpe*Ce;
%plot(tran_cpe);
cos=1;
%pierwsza czesc
for i=n0+1:n0+n1
    Ire(i)=Ucpe(i-1)/Re;
    Icpe(i)=A1-Ire(i);
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i)+A1*Rs;
end
i;
%druga czesc
for i=n0+n1+1:n0+n1+n2
    Ire(i)=Ucpe(i-1)/Re;
    Icpe(i)=A2-Ire(i);
    %Icpe(i)
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i)+A2*Rs;
end
i;
%trzeca czesc
for i=n0+n1+n2+1:n0+n1+n2+n3
    Ire(i)=Ucpe(i-1)/Re;
    Icpe(i)=A3-Ire(i);
    %Icpe(i)
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i)+A3*Rs;
end
i;
%przerwa miedzy impulsem i rozladowywaniem
for i=n0+n1+n2+n3+1:n0+n1+n2+n3+n_delay
    Ire(i)=Ucpe(i-1)/Re;
    Icpe(i)=-Ire(i);
    %Icpe(i)
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i);
end
i;
%rozladowywanie
for i=n0+n1+n2+n3+n_delay+1:n0+n1+n2+n3+n_delay+n_disch
    Ire(i)=Ucpe(i-1)/Re;
    Irs(i)=(Vdc-Ucpe(i-1))/Rs;    
    Icpe(i)=Irs(i)-Ire(i);
    %Icpe(i)
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i)+Irs(i)*Rs;
end
i;

for i=n0+n1+n2+n3+n_delay+n_disch+1:N
    Ire(i)=Ucpe(i-1)/Re;
    Icpe(i)=-Ire(i);
    %Icpe(i)
    u0=0;
    u1=0;
    for k=0:i-cos
        u0=u0+Icpe(i-k)*dt*tran_cpe(k+1);  %aktualne napiecie, uwzgledniajace rowinez 'zastrzyk' pradu z aktualnego kroku czasowego                                        
        %u1=u1+Icpe(i-k)*tran_cpe(k+2);  
    end
    Ucpe(i)=u0;
    U(i)=Ucpe(i);
end

w=Ucpe;
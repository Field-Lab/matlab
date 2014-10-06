function model=fit_el_mod1(pomiar,model_pocz);
Ce_krok=2e+8;
beta_krok=0.005;
Re_krok=100000;
warburg=1;
Rs_krok=200;

r=pomiar(:,2).*cos(pomiar(:,3)/90*pi/2);
ur=pomiar(:,2).*sin(pomiar(:,3)/90*pi/2);

imp=(r+i*ur)*1000;

f=pomiar(:,1);

Ce=model_pocz.Ce;
beta=model_pocz.beta
Re=model_pocz.Re;
Rs=model_pocz.Rs;
figure(3)

for i=1:200
    %Ce_krok=Ce_krok*0.99;
    %Ce_krok=Ce_krok*1.01;
    %beta_krok=beta_krok*0.99;
    %Re_krok=Re_krok*0.99;
    %Rs_krok=Rs_krok*0.99;
    
    Ce0=[Ce-Ce_krok Ce Ce+Ce_krok];
    beta0=[beta-beta_krok beta beta+beta_krok];
    Re0=[Re-Re_krok Re Re+Re_krok];
    Rs0=[Rs-Rs_krok Rs Rs+Rs_krok];        
    
    for j=0:80;
        iCe=floor(j/27);
        ibeta=floor((j-iCe*27)/9);
        iRe=floor((j-iCe*27-ibeta*9)/3);
        iRs=floor((j-iCe*27-ibeta*9-iRe*3));
        
        z=el_imp_mod1(Ce0(iCe+1),beta0(ibeta+1),Re0(iRe+1),warburg,Rs0(iRs+1),f);      
        %warburg na razie nieuzywany w funkcji el_imp_mod1, do pozniejszej
        %implmentacji
        blad0=abs(z-imp).^2;
        blad1=abs(imp).^2;                
        blad(j+1)=sum(blad0./blad1);
        
        blad0r=real(z-imp).^2;
        blad1r=real(imp).^2;
        bladr=sum(blad0r./blad1r);        
        blad0i=imag(z-imp).^2;
        blad1i=imag(imp).^2;
        bladi=sum(blad0i./blad1i);
        blad(j+1)=sqrt(bladr^2+bladi^2);
        
        %loglog(f,abs(imp),f,abs(z))
        %semilogx(angle(z))
    end
    blad
    indeks=find(blad==min(blad))-1
    iCe=floor(indeks/27);
    ibeta=floor((indeks-iCe*27)/9);
    iRe=floor((indeks-iCe*27-ibeta*9)/3);
    iRs=floor((indeks-iCe*27-ibeta*9-iRe*3));
    
    Ce=Ce0(iCe+1);
    beta=beta0(ibeta+1)
    Re=Re0(iRe+1);
    Rs=Rs0(iRs+1);
    z=el_imp_mod1(Ce,beta,Re,warburg,Rs,f);
    subplot(1,2,1);
    s=loglog(f,abs(z),f,abs(imp),'g*');
    grid on;
    subplot(1,2,2);
    semilogx(f,angle(z)*57,f,angle(imp)*57,'g*');
    grid on
    pause(0.01);
    blad(indeks+1);
end   
grid on;
Ce
beta
Re
Rs 
model=blad;
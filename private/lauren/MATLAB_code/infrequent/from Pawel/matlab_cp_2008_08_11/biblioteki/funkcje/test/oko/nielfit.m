function [beta,przebieg] = nielfit(t,x,model,beta0,delta)
%NLINFIT Nonlinear least-squares data fitting by the Gauss-Newton method.
%   NLINFIT(X,Y,'MODEL',BETA0) finds the coefficients of the nonlinear 
%   function described in MODEL. MODEL is a user supplied function having 
%   the form y = f(beta,x). That is MODEL returns the predicted values of y
%   given initial parameter estimates, beta, and the independent variable, X.   
%   [BETA,R,J] = NLINFIT(X,Y,'MODEL',BETA0) returns the fitted coefficients
%   BETA the residuals, R, and the Jacobian, J, for use with NLINTOOL to
%   produce error estimates on predictions.

%   B.A. Jones 12-06-94.
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.12 $  $Date: 1998/09/09 19:39:39 $

if min(size(t)) ~= 1
    error('Wektor "t" musi byc tablica jednowymiarowa');
end
t=t(:)';

if min(size(x)) ~= 1
    error('Wektor "t" musi byc tablica jednowymiarowa');
end
x=x(:)';

if min(size(beta0)) ~= 1
    error('Wektor "beta0" musi byc tablica jednowymiarowa');
end
beta0=beta0(:)';

%if min(size(delta)) ~= 1
%    error('Wektor "delta" musi byc tablica jednowymiarowa');
%end
delta=delta(:)';


if size(t)~=size(x)
    error('Rozmiary tablic "t" oraz "x" musza byc identyczne');
end

if size(beta0)~=size(delta)
    error('Rozmiary tablic "t" oraz "x" musza byc identyczne');
end

beta=beta0;
beta2=beta0;
gradient=zeros(size(beta));
y=feval(model,beta,t);

for i=1:5000   % NA RAZIE !!!!
    %y=feval(model,beta,t);
    norma=mean((y-x).^2);
    gradient_old=gradient;
    %beta2=beta;
    %delta=delta*0.998;
    for j=1:length(beta)
        beta2=beta;
        beta2(j)=beta2(j)+delta(j);
        %yplus = feval(model,beta+delta,X);
        yprim=feval(model,beta2,t);
        normaprim=mean((yprim-x).^2);
        beta3=beta2;
        beta3(j)=beta(j)-delta(j);
        ybis=feval(model,beta3,t);
        normabis=mean((ybis-x).^2);
        gradient(1,j)=sign(round((normabis-normaprim)/norma*10000));
    end;
    
    %beta2
    %gradient
    %plot(t,x,t,feval(model,beta,t));
    beta=beta+delta.*gradient;
    przebieg(i,:)=beta;
    y=feval(model,beta,t);
    if gradient_old == -gradient
        i
        break;
    end
end


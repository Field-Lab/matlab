function out = cacRegisterCertificate()
%cacRegisterCertificate() - associates a MyProxy Certificate with a users cac
%user name.  The user will need a valid CAC account in order to complete
%this process.  This process ONLY needs to be done once for a user!
%
% See also cacsched, cacscheduler
%
% Copyright 2009-2010 Cornell Center for Advanced Computing

%Ensure we clear all CertManager stuff first.
CM = edu.cornell.cac.tuc.cacscheduler.globus.CertManager.getInstance();
CM.clearCreds();

%No register certificate
CR = edu.cornell.cac.cacscheduler.tools.CertReg();
CR.registerCert();
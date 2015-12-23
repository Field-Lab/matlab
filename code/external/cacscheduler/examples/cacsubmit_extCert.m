function j = cacsubmit_extCert(sched)
%CACSUBMIT_EXTCERT(SCHED) - part of the CAC Scheduler test suite.  This tests the ability to supply an already-available 
% x509 certificate from either a path location or from an path stored in an environment variable.

%Setup the CM to use a stored certificate
cm = edu.cornell.cac.tuc.cacscheduler.globus.CertManager.getInstance();
%cm.setExternalCredential('C:\DOCUME~1\naw47\LOCALS~1\Temp\x509up_u500');
%or
cm.setCredentialFromEnv('x509_home');

j = createJob(sched);
createTask(j,@rand,1,{3,3});
createTask(j,@rand,1,{3,3});
submit(j);
waitForState(j);
a = getAllOutputArguments(j);
a{1}

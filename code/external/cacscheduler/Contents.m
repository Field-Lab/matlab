% MATLAB
%
% Files
%   cacCancelJob                  - CACCANCELJOB() - cancels a job on the TUC cluster at CAC.  The end user
%   cacCancelTask                 - CACCANCELTASK() - cancels a job on the TUC cluster at CAC.  The end user
%   cacColSum                     - cacColSum() - example parallel job worker function.  See examples for
%   cacCopyJobFilesIfFinished     - CACCOPYJOBFILESIFFINISHED() Copy the job files of a finished job to the
%   cacDestroyJob                 - CACDESTROYJOB() - deletes all files related to the provided job locally 
%   cacGetJobState                - CACGetJobState() - Retrieves the state of a job running on TUC.  The end user
%   cacNonSharedParallelDecodeFcn - cacNONSHAREDPARALLELDECODEFCN Prepares a worker to run a MATLAB task.
%   cacNonSharedParallelSubmitFcn - CACNONSHAREDPARALLELSUBMITFCN() - PCT interface function to submit a
%   cacNonSharedSimpleDecodeFcn   - Prepares a worker to run a MATLAB task.
%   cacNonSharedSimpleSubmitFcn   - CACNONSHAREDSIMPLESUBMITFCN() - PCT interface function to submit a
%   cacRegisterCertificate        - cacRegisterCertificate() - associates a MyProxy Certificate with a users cac
%   cacsched                      - CACSCHED() - construct a scheduler object to interface with TUC located at
%   copyDataFromCAC               - COPYDATAFROMCAC() - PCT interface function to move data from TUC to the
%   copyDataToCAC                 - COPYDATATOCAC() - PCT interface function to move data from the local machine
%   getErrors                     - GETERRORS(job) - Pretty print or return the errors associated with a job.
%   getOutput                     - GETOUTPUT(job) - Pretty print or return the command window output for all
%   cacscheduler                    - MATLAB
%   runtests                      - runTests(sched) - test function to test cacscheduler Install.
%   sendDirToCAC                  - SENDDIRTOCAC(localdir, remoteName) - send the directory at localdir to the
%   sendFileToCAC                 - SENDFILETOCAC(localfile) - send a single local file to TUC storage system

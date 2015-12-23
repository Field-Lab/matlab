function cac_sendmail(fromaddr,toaddr,subject,message,attachments)
%CAC_SENDMAIL Send e-mail.
%   CAC_SENDMAIL(FROMADDR,TOADDR,SUBJECT,MESSAGE,ATTACHMENTS) sends an e-mail.
% sends mail using Cornell SMTP server, be sure to specify fromaddr
% see SENDMAIL for documentation


% Argument parsing.
error(nargchk(3,5,nargin,'struct'));
if ischar(toaddr)
    toaddr = {toaddr};
end
if (nargin < 4)
    message = '';
end
if (nargin < 5) 
    attachments = [];
elseif ischar(attachments)
    attachments = {attachments};
end
setpref('Internet','E_mail',fromaddr);
setpref('Internet','SMTP_Server','appsmtp.mail.cornell.edu');
message = [message sprintf('\nMATLAB: %s\nHPC Job ID#: %s\n',getenv('MDCE_TASK_LOCATION'),getenv('CCP_JOBID'))];
sendmail(toaddr,subject,message,attachments);
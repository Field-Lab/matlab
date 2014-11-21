% made a dummy email account to
% example usage: gmail('crhoades227@gmail.com', sprintf('Done with %s %s_data_run_%02d_config_%d',run_opt.data_set, run_opt.cell_type, run_opt.data_run, run_opt.config_num))

function []=gmail(address,subject,message)

% Define these variables appropriately:
mail = 'crhoadesDA@gmail.com'; %Your gmail email address (DA = dummy account)
password = 'dummyaccount'; %gmail password

% Then this code will set up the preferences properly:
setpref('Internet','E_mail',mail);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',mail);
setpref('Internet','SMTP_Password',password);
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');

% Send the email
if nargin == 3
    sendmail(address,subject,message)
elseif nargin == 2
    sendmail(address,subject)
end

end
function getCACCerts()

SSLCert = 'https://myproxy.cac.cornell.edu/CA/fb4341f6.0';
SSLCertSign = 'https://myproxy.cac.cornell.edu/CA/fb4341f6.signing_policy';
TRUSTED_CERT_PATH =  fullfile('.globus','certificates');

fprintf('Making new cert directory.\n');
f = java.io.File(java.lang.System.getProperty('user.home'), TRUSTED_CERT_PATH);
f.mkdirs();
path = char(f.getAbsolutePath());
fprintf('Using %s for storing certs.\n', path);

cert = urlread(SSLCert);
sign = urlread(SSLCertSign);

[fid,message] = fopen(fullfile(path,'fb4341f6.0'),'w');
if fid < 0
    fprintf('Failed opening cert file:%s',message);
end
fwrite(fid,cert);
fclose(fid);

fid = fopen(fullfile(path,'fb4341f6.signing_policy'),'w');
if fid < 0
    fprintf('Failed opening cert file:%s',message);
end
fwrite(fid,sign);
fclose(fid);

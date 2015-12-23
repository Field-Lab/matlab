function cac_addpath(override)
CACHome = fileparts(which('runtests'));
examples = fullfile(CACHome,'examples');
pct = fullfile(CACHome,'PCT');
contrib = fullfile(CACHome,'contrib');
sf = java.io.File('startup.m');
addpath(CACHome);
addpath(pct);
addpath(contrib);
if sf.exists() | nargin == 1
    if ~nargin==1
    fprintf(['startup.m exists\nPlease add to startup.m\n %%START\n\' ...
             'n']);
    fprintf('java.lang.System.setProperty(''sun.security.ssl.allowUnsafeRenegotiation'',''true'');\n');
    fprintf('addpath(''%s'');\n',CACHome);
    fprintf('addpath(''%s'');\n',pct);
    fprintf('addpath(''%s'');\n\n %%END',contrib);
    fprintf('\n\nto edit run >> edit startup\n');
    input('press any key to continue\n','s');
    end
else
    fprintf('creating startup.m\n');
    fout = fopen('startup.m','w');
    fprintf(fout,'java.lang.System.setProperty(''sun.security.ssl.allowUnsafeRenegotiation'',''true'');\n');
    fprintf(fout,'addpath(''%s'');\n',CACHome);
    fprintf(fout,'addpath(''%s'');\n',pct);
    fprintf(fout,'addpath(''%s'');\n',contrib);
    fclose(fout);
end
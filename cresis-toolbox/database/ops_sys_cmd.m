% OPS SYSTEM CONNECTION SETTINGS
%
% profile_cmd: (logical) should the python profilier be run.
% sys_url: (string) url to current system.
%
% Author: Kyle W. Purdon

% USER INPUT
% ==========================================================
profile_cmd = false;

sys_url = 'https://ops.cresis.ku.edu/';
% sys_url = 'http://192.168.111.222/';

% AUTOMATIC SECTION (DONT MODIFY)
% ============================================================

db_user = 'doesnt';
db_pswd = 'doanything';

sys_name = 'ops';
server_url = strcat(sys_url,sys_name,'/');
geoserver_url = strcat(sys_url,'geoserver/');

if profile_cmd
    web(strcat(sys_url,'profile-logs/'));
end

clear sys_url
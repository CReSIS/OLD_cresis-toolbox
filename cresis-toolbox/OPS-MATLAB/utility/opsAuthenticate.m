function [opsParam,opsAuth,opsProfile] = opsAuthenticate(param,forceLogin)
% input param from any ops function
% output param with user info added
% forces login if not logged in
%
% more docs coming soon
%
% Author: Kyle W. Purdon

opsParam = param;

global gRadar;
opsCmd;

if nargin < 2
  forceLogin =  false;
end

if forceLogin
  % FORCE USER TO LOGIN
  [status,loginNotice] = opsLoginUser();
  if status == 1
    opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
    opsProfile = load(fullfile(gRadar.tmp_path,'ops.profile.mat'));
    opsParam.properties.userName = opsAuth.userName;
    opsParam.properties.isAuthenticated = opsAuth.isAuthenticated;
    opsParam.properties.mat = true;
  elseif status ~= 1
    warning(loginNotice);
  end
else
  if exist(fullfile(gRadar.tmp_path,'ops.mat'),'file') && exist(fullfile(gRadar.tmp_path,'ops.profile.mat'),'file')
    % A USER HAS ALREADY LOGGED IN
    opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
    if opsAuth.isAuthenticated
      % THE USER IS ALREADY AUTHENTICATED
      opsProfile = load(fullfile(gRadar.tmp_path,'ops.profile.mat'));
      opsParam.properties.userName = opsAuth.userName;
      opsParam.properties.isAuthenticated = opsAuth.isAuthenticated;
      opsParam.properties.mat = true;
    else
      % THE USER IS NOT AUTHENTICATED
      [status,loginNotice] = opsLoginUser();
      if status == 1
        opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
        opsProfile = load(fullfile(gRadar.tmp_path,'ops.profile.mat'));
        opsParam.properties.userName = opsAuth.userName;
        opsParam.properties.isAuthenticated = opsAuth.isAuthenticated;
        opsParam.properties.mat = true;
      elseif status ~= 1
        warning(loginNotice);
      end
    end
  else
    % NO USER IS LOGGED IN, LOG A USER IN
    [status,loginNotice] = opsLoginUser();
    if status == 1
      opsAuth = load(fullfile(gRadar.tmp_path,'ops.mat'));
      opsProfile = load(fullfile(gRadar.tmp_path,'ops.profile.mat'));
      opsParam.properties.userName = opsAuth.userName;
      opsParam.properties.isAuthenticated = opsAuth.isAuthenticated;
      opsParam.properties.mat = true;
    elseif status ~= 1
      warning(loginNotice);
    end
  end
end

end
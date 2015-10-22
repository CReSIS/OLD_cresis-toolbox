% =========================================================================
% OPS TEST UTILITY
%
% Tests that MATLAB can communicate with the OPS system.
%
% A user must exist or be created to run all tests.
%
% Author: Kyle W. Purdon, Trey Stafford
%
% =========================================================================
opsCmd;

warn = false;
%Make sure MATLAB can query the database:
[status,data] = opsQuery('select 9;');
if status == 1
    if data{1} == 9
        fprintf('Database successfully queried.\n')
    else
        warning('%s s% \n','Something went wrong with querying the database. Expected 9 but got ', data{1})
        warn = true;
    end
else
    warning('%s %s \n','Something went wrong with querying the database! Got status ',status)
    warn = true;
end
clear status data

% Stop process to see if a user should be created.
if ~any(strcmp(gOps.serverUrl,{'https://ops.cresis.ku.edu/ops/','http://ops2.cresis.ku.edu/ops/'}))
    confirmParams = {'Do you wish to create a new user?','','Note: for these tests to continue a user must be defined.'};
    confirmButton = questdlg(confirmParams,'Create User?','YES: Create User','NO: Use Existing User','NO: Quit Tests','YES: Create User');
    switch confirmButton
        case 'NO: Quit Tests'
            error('TESTS ENDED BY USER.');
        case 'YES: Create User'
            opsCreateUser();            
    end
end
       
%Make sure MATLAB can handle a django view response
[status,~] = opsGetSystemInfo();
if status == 1
    fprintf('opsGetSystemInfo() successful! Django is responsive. \n')
else
    warning('%s \n','Something went wrong with opsGetSystemInfo()!')
    warn = true;
end
clear status 

%Make sure the database has a surface layer defined for rds (should be
%added automatically by django fixture)
[status,data] =  opsGetLayers('rds');
if status == 1
    if any(ismember(data.properties.lyr_id,1)) && any(ismember(data.properties.lyr_id,2))
        fprintf('opsGetLayers() successfully shows that there is a surface and bottom layer. \n');
    else
        warning('opsGetLayers() seemed to have worked, but there is no surface and bottom layer! Something may have gone wrong with syncdb.')
        warn = true;
    end
else
    warning('Something went wrong with opsGetLayers()!')
    warn = true;
end
clear status data

%Logout the user.
[status,logoutNotice] = opsLogoutUser();
if status ~= 1
  warning(logoutNotice);
  warn = true;
else
  fprintf('%s\n',logoutNotice)
end
clear status logoutNotice

%Let the user know the tests have completed. 
if ~warn
    fprintf('opsTests have completed. No Errors were encountered. \n');
else
    warning('opsTests have completed but errors were encountered!');
end
clear warn
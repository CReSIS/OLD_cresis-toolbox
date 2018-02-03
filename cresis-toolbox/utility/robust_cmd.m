function [success] = robust_cmd(cmd,num_attempts)

success = false;
attempts = 0;
while attempts < num_attempts
  try
    evalin('caller',cmd);
    success = true;
    return
  catch ME
    ME.getReport
    cmd
    warning('Function failed, pausing 30 sec: %s', cmd);
    pause(30);
  end
end

function [success] = robust_cmd(cmd,num_attempts)

success = false;
attempt = 1;
while attempt <= num_attempts
  try
    evalin('caller',cmd);
    success = true;
    return
  catch ME
    ME.getReport
    cmd
    warning('cmd failed on attempt %d of %d: %s', attempt, num_attempts, cmd);
    pause(10);
  end
  attempt = attempt + 1;
end

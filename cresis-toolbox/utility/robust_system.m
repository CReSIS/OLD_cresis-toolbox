function [status,result] = robust_system(cmd)

status = -1;
attempts = 0;
while status ~= 0
  try
    [status,result] = system(cmd);
  catch ME
    ME.getReport
    cmd
    warning('system call failed');
    keyboard
  end
  if status ~= 0
    warning('%s\n  FAILED %d %s', cmd, status, result);
    attempts = attempts + 1;
    delay_period = 3*2^(attempts-1);
    fprintf('  Delaying %d seconds\n', delay_period)
    pause(delay_period);
    if attempts > 10
      % There is potentially something wrong with the command if it
      % fails this many times. Look into it before running "dbcont".
      keyboard;
    end
  end
end

function cmd = peak(obj)
% Peaks at the current command

if obj.pointer > 0
  cmd = obj.stack(obj.pointer);
else
  cmd = [];
end

end

function name_enc = ni_xml_name_encode(name)
% name_enc = ni_xml_name_encode(name)
%
% Letters and Numbers pass unencoded
% All other letters are encoded as:
%   ZHH
% where HH is the hexidecimal representation of the ASCII value for
% that letter.
%
% Author: John Paden
%
% See also: print_ni_xml_object.m, ni_xml_name_decode.m,
%   ni_xml_name_encode.m, read_ni_xml_object.m, write_ni_xml_object.m

name_enc = '';
out_idx = 0;
for idx = 1:length(name)
  if (name(idx) >= 'A' && name(idx) <= 'Y') ...
      || (name(idx) >= 'a' && name(idx) <= 'z')
    out_idx = out_idx + 1;
    name_enc(out_idx) = name(idx);
  else
    out_idx = out_idx + 1;
    name_enc(out_idx:out_idx+2) = sprintf('Z%s', dec2hex(name(idx)));
    out_idx = out_idx + 2;
  end
end

end

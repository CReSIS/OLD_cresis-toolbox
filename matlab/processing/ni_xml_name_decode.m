function name = ni_xml_name_decode(name_enc)
% name = ni_xml_name_decode(name_enc)
%
% Letters and Numbers pass undecoded EXCEPT 'Z'
% All other letters are decoded as:
%   ZHH
% where HH is the hexidecimal representation of the ASCII value for
% that letter. Therefore "Z20" gets decoded into a space character " ".
%
% Author: John Paden
%
% See also: print_ni_xml_object.m, ni_xml_name_decode.m,
%   ni_xml_name_encode.m, read_ni_xml_object.m, write_ni_xml_object.m

name = '';
idx = 1;
out_idx = 0;
while idx <= length(name_enc)
  out_idx = out_idx + 1;
  if (name_enc(idx) >= 'A' && name_enc(idx) <= 'Y') ...
      || (name_enc(idx) >= 'a' && name_enc(idx) <= 'z')
    name(out_idx) = name_enc(idx);
    idx = idx + 1;
  else
    name(out_idx) = sprintf('%s', char(hex2dec(name_enc(idx+(1:2)))));
    idx = idx + 3;
  end
end

end

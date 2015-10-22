function [freq,sparam] = create_snp_from_s2p(fn_base)
% [freq,sparam] = create_snp_from_s2p(fn_base);
%
% Takes S2P files saved with filenames [fn_base][port1][port2].s2p
% and converts them to an SNP file. For example, if
% fn_base = 'A', then filenames should be A12.s2p, A13.s2p, A23.s2p,
% etc.
%
% sparam(port1 axes,port2 axes,freq axes)

[fn_base_dir,fn_base_name] = fileparts(fn_base);
fns = get_filenames(fn_base_dir,fn_base_name,'','.s2p');

sparam = [];
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  [~,fn_name] = fileparts(fn);
  
  if length(fn_base_name) ~= length(fn_name)-2
    % File name format does not match
    continue;
  end
  
  %fprintf('Reading %s\n', fn);
  
  ports(1) = str2double(fn_name(length(fn_base_name)+1));
  ports(2) = str2double(fn_name(length(fn_base_name)+2));
  
  if 0
    % Use cresis toolbox reader
    [freq,sparam_s2p] = SXPParse(fn);
  else
    % Use Matlab toolbox reader
    a1 = read(rfdata.data,fn);
    freq = a1.Freq;
    sparam_s2p = a1.S_Parameters;
  end
  
  % Write the S2P params into the correct part of the SNP matrix
  for port1_s2p=1:2
    for port2_s2p=1:2
      sparam(ports(port1_s2p),ports(port2_s2p),:) = sparam_s2p(port1_s2p,port2_s2p,:);
    end
  end
end

end

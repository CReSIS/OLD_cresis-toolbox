function [xo_table_tag, idx_xo, reuse_loc, xo_hdr] = de_ident(ident)

global gRadar;

xo_table_tag = extractBefore(ident,'_xo_');
idx_xo = str2double(extractAfter(ident,'_xo_'));
reuse_loc = fullfile(gRadar.ct_tmp_path, 'radcal_mat', xo_table_tag);
xo_hdr = sprintf('(xo %d) ', idx_xo);

end
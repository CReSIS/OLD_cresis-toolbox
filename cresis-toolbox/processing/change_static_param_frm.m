function change_static_param_frm(old_filename, new_filename, frm_str)

if ~ischar(frm_str)
    error('Frames Input has to be a character with comma separated values. Please try again!');
end

frms = cellfun(@str2num,strsplit(frm_str,','));

load(old_filename,'static_param');
static_param.cmd.frms = frms;
save(new_filename,'static_param');
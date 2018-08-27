function [success] = csarp_task_ollie(static_param_file_name,dynamic_param_file_name,frm,chunk,wf,adc)
% [success] = csarp_task_ollie(static_param_file_name,dynamic_param_file_name,frm,chunk,wf,adc)

% Test
if (ischar(frm))
  frm=str2num(frm);
end
if (ischar(chunk))
  chunk=str2num(chunk);
end
if (ischar(wf))
  wf=str2num(wf);
end
if (ischar(adc))
  adc=str2num(adc);
end
  
load(static_param_file_name,'static_param');
param=static_param;
load(dynamic_param_file_name,'dynamic_param');
param.load.frm = frm;
param.load.chunk_idx = chunk;
param.load.num_chunks = dynamic_param.frms.(['frm',num2str(frm)]).num_chunks;
param.load.recs = dynamic_param.frms.(['frm',num2str(frm)]).chunks.(['chunk',num2str(chunk)]).recs;
param.load.imgs = {[wf adc]};
param.load.sub_band_idx = 1;

clear frm chunk wf adc;
[success] = csarp_task(param);
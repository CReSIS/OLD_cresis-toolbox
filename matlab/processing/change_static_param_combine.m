function change_static_param_combine(old_filename,new_filename,new_method,out_path,ch1,ch2,ch3,ch4,ch5,ch6,ch7,ch8,ch9,ch10,ch11,ch12,ch13,ch14,ch15,ch16,ch17,ch18,ch19,ch20,ch21,ch22,ch23,ch24)

if(~isequal(new_method,'standard') && ~isequal(new_method,'mvdr'))
    error('Method selection is not valid. Choose "standard" or "mvdr"');
end

if (ischar(ch1))
  ch1=str2num(ch1);
end
if (ischar(ch2))
  ch2=str2num(ch2);
end
if (ischar(ch3))
  ch3=str2num(ch3);
end
if (ischar(ch4))
  ch4=str2num(ch4);
end
if (ischar(ch5))
  ch5=str2num(ch5);
end
if (ischar(ch6))
  ch6=str2num(ch6);
end
if (ischar(ch7))
  ch7=str2num(ch7);
end
if (ischar(ch8))
  ch8=str2num(ch8);
end
if (ischar(ch9))
  ch9=str2num(ch9);
end
if (ischar(ch10))
  ch10=str2num(ch10);
end
if (ischar(ch11))
  ch11=str2num(ch11);
end
if (ischar(ch12))
  ch12=str2num(ch12);
end
if (ischar(ch13))
  ch13=str2num(ch13);
end
if (ischar(ch14))
  ch14=str2num(ch14);
end
if (ischar(ch15))
  ch15=str2num(ch15);
end
if (ischar(ch16))
  ch16=str2num(ch16);
end
if (ischar(ch17))
  ch17=str2num(ch17);
end
if (ischar(ch18))
  ch18=str2num(ch18);
end
if (ischar(ch19))
  ch19=str2num(ch19);
end
if (ischar(ch20))
  ch20=str2num(ch20);
end
if (ischar(ch21))
  ch21=str2num(ch21);
end
if (ischar(ch22))
  ch22=str2num(ch22);
end
if (ischar(ch23))
  ch23=str2num(ch23);
end
if (ischar(ch24))
  ch24=str2num(ch24);
end

channel_selection = [];
if(ch1==1) channel_selection = [channel_selection 1]; end
if(ch2==1) channel_selection = [channel_selection 2]; end
if(ch3==1) channel_selection = [channel_selection 3]; end
if(ch4==1) channel_selection = [channel_selection 4]; end
if(ch5==1) channel_selection = [channel_selection 5]; end
if(ch6==1) channel_selection = [channel_selection 6]; end
if(ch7==1) channel_selection = [channel_selection 7]; end
if(ch8==1) channel_selection = [channel_selection 8]; end
if(ch9==1) channel_selection = [channel_selection 9]; end
if(ch10==1) channel_selection = [channel_selection 10]; end
if(ch11==1) channel_selection = [channel_selection 11]; end
if(ch12==1) channel_selection = [channel_selection 12]; end
if(ch13==1) channel_selection = [channel_selection 13]; end
if(ch14==1) channel_selection = [channel_selection 14]; end
if(ch15==1) channel_selection = [channel_selection 15]; end
if(ch16==1) channel_selection = [channel_selection 16]; end
if(ch17==1) channel_selection = [channel_selection 17]; end
if(ch18==1) channel_selection = [channel_selection 18]; end
if(ch19==1) channel_selection = [channel_selection 19]; end
if(ch20==1) channel_selection = [channel_selection 20]; end
if(ch21==1) channel_selection = [channel_selection 21]; end
if(ch22==1) channel_selection = [channel_selection 22]; end
if(ch23==1) channel_selection = [channel_selection 23]; end
if(ch24==1) channel_selection = [channel_selection 24]; end

load(old_filename,'static_param');

static_param.combine.imgs = {[ones([size(channel_selection,2) 1]),channel_selection.'],[2*ones([size(channel_selection,2) 1]),channel_selection.'],[3*ones([size(channel_selection,2) 1]),channel_selection.']};
static_param.combine.method = new_method;
static_param.combine.out_path = out_path;

save(new_filename,'static_param');

return;
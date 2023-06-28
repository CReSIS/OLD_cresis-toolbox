% this function is used to correct and interpolate unreasonable time
% sequence---qishi
function time_after=interpolation(time_before,previous_mark,previous_time)
time=time_before;
if previous_mark
  second_valid_idx=find(time>time(1),1,'first');
  time(1:second_valid_idx)=interp1([0  second_valid_idx],[previous_time time(second_valid_idx)],[1:1:second_valid_idx]);
end

while ~isempty(find(diff(time)<0))% deal with decrease time points
  idx=find(diff(time)<0,1,'last');
  time(idx)=time(idx+1);
end

if (any(diff(time)==0))&&( ~all( (time-time(1))==0))% deal with those "flat" time points,and there is a rare situation that all time points are equal
    jump_idx=find(diff(time)>0);
    jump_idx=jump_idx+1;
    if jump_idx(end)~=length(time)%we need extrapolation to interp last flat time sequence
      last_mark=1;
    else
      last_mark=0;
    end
    jump_idx=[1,jump_idx];% valid time point
    for jump_idx_idx=2:length(jump_idx)
        lseg=jump_idx(jump_idx_idx)-jump_idx(jump_idx_idx-1);%length between 2 valid point
        if lseg>1% means there is a "flat" time (i.e. repeated time)
            time(jump_idx(jump_idx_idx-1):jump_idx(jump_idx_idx))=interp1([jump_idx(jump_idx_idx-1),jump_idx(jump_idx_idx)],[time(jump_idx(jump_idx_idx-1)),time(jump_idx(jump_idx_idx))],jump_idx(jump_idx_idx-1):1:jump_idx(jump_idx_idx),'linear');       
        end
    end
    if last_mark %interp last flat time sequence
      time(jump_idx(end-1):end) = interp1( [1:jump_idx(end)-jump_idx(end-1)+1],time(jump_idx(end-1):jump_idx(end)),[1:length(time)-jump_idx(end-1)+1],'linear','extrap' );
    end      
end  
time_after=time;
end


     

    
    
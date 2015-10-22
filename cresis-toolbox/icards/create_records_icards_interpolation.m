function time_after=interpolation(time_before)
time=time_before;
if any(diff(time)<=0)
    jump_idx=find(diff(time)>0);
    jump_idx=jump_idx+1;
    jump_idx=[1,jump_idx];
    for jump_idx_idx=2:length(jump_idx)
        lseg=jump_idx(jump_idx_idx)-jump_idx(jump_idx_idx-1);
        if lseg>1
            time(jump_idx(jump_idx_idx-1):jump_idx(jump_idx_idx))=interp1([jump_idx(jump_idx_idx-1),jump_idx(jump_idx_idx)],[time(jump_idx(jump_idx_idx-1)),time(jump_idx(jump_idx_idx))],jump_idx(jump_idx_idx-1):1:jump_idx(jump_idx_idx),'linear');
        
        end
    end    
    if jump_idx(end)~=length(time)
        interpolation_average=(time(jump_idx(1))-time(jump_idx(end)))/(jump_idx(1)-jump_idx(end));
        for ii=1:length(time)-jump_idx(end)
            time(jump_idx(end)+ii)=time(jump_idx(end)+ii)+interpolation_average*ii;
        end
    end       
end
time_after=time;
end


     

    
    
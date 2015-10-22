%NOTES:This function merges TRAJ and NMEA data and generate a .csv file for
%a certain day as the output. TRAJ file is the main data resources and NMEA
%is used as complemtary. But actually, TRAJ file is not available for most
%of those years from 1993 to 2002, in this case, only NMEA is used to generate 
%the final .csv file. As I didn't find any multiple TRAJ files of one cetain
%day, this function does not take the situation of merging multiple TRAJ files 
%into consideration.
icards_gps_loaddata; %Call function to load data
out_fn='X:\metadata\';
TRAJ_MAX_TIMEGAP=1;%set up traj file time gap value 

if exist('TRAJ')%we have both of traj and nmea files of a certain day

 if NMEA.time(1)<TRAJ.gps_time(1)%add NMEA data points to the left of TRAJ data
   [~,start_idx]=find(NMEA.time<TRAJ.gps_time(1));
   TRAJ.gps_time=[NMEA.time(start_idx),TRAJ.gps_time];
   TRAJ.lon=[NMEA.lon(start_idx),TRAJ.lon];
   TRAJ.lat=[NMEA.lat(start_idx),TRAJ.lat];
   TRAJ.elev=[NMEA.elev(start_idx),TRAJ.elev];
   TRAJ.roll=[NMEA.roll(start_idx),TRAJ.roll];
   TRAJ.pitch=[NMEA.pitch(start_idx),TRAJ.pitch];
   TRAJ.heading=[NMEA.heading(start_idx),TRAJ.heading];
 end
 
 if TRAJ.gps_time(end)<NMEA.time(end) %add NMEA data points to the right of TRAJ data
   [~,end_idx]=find(NMEA.time>TRAJ.gps_time(end));
   TRAJ.gps_time=[TRAJ.gps_time,NMEA.time(end_idx)];
   TRAJ.lon=[TRAJ.lon,NMEA.lon(end_idx)];
   TRAJ.lat=[TRAJ.lat,NMEA.lat(end_idx)];
   TRAJ.elev=[TRAJ.elev,NMEA.elev(end_idx)];
   TRAJ.roll=[TRAJ.roll,NMEA.roll(end_idx)];
   TRAJ.pitch=[TRAJ.pitch,NMEA.pitch(end_idx)];
   TRAJ.heading=[TRAJ.heading,NMEA.heading(end_idx)];
 end
 
  traj_time_gaps_idx = find(abs(diff(TRAJ.gps_time)) > TRAJ_MAX_TIMEGAP);%find out those gaps in TRAJ file
  length_plus=0;
  
  for ii=1:length(traj_time_gaps_idx)%fill these gaps with NMEA gap by gap
    left_time=TRAJ.gps_time(traj_time_gaps_idx(ii)+length_plus);%the start point of a gap
    right_time=TRAJ.gps_time(traj_time_gaps_idx(ii)+1+length_plus);%the end point of a gap
    patch=intersect(find(NMEA.time>left_time),find(NMEA.time<right_time));%"patch" is a block of time sequence from NMEA time to fill the gap
    
      if isempty(patch)%sometimes there's no time block in NMEA time can be used to fill the gap, in this case I leave the gap there without doing any thing       
        TRAJ.gps_time=TRAJ.gps_time;
        TRAJ.lon=TRAJ.lon;
        TRAJ.lat=TRAJ.lat;
        TRAJ.elev=TRAJ.elev;
        TRAJ.roll=TRAJ.roll;
        TRAJ.pitch=TRAJ.pitch;
        TRAJ.heading=TRAJ.heading;
      else             %if the "patch" exists, then fill the data gap  
        TRAJ.gps_time=[TRAJ.gps_time(1:traj_time_gaps_idx(ii)+length_plus),NMEA.time(patch),TRAJ.gps_time(traj_time_gaps_idx(ii)+1+length_plus:end)];
          if  TRAJ.gps_time(traj_time_gaps_idx(ii)+1+length_plus)-TRAJ.gps_time(traj_time_gaps_idx(ii)+length_plus)>20   % if the time gap is longer than 20 secs           
            [~,idx1]=find(NMEA.time(patch)<TRAJ.gps_time(traj_time_gaps_idx(ii)+length_plus)+10);
            head=patch(1:idx1(end));
            [~,idx2]=find(NMEA.time(patch)>TRAJ.gps_time(traj_time_gaps_idx(ii)+1+length_plus)-10);
            tail=patch(idx2(1):end);
            lon_patch_head=linspace(TRAJ.lon(traj_time_gaps_idx(ii)+length_plus)-NMEA.lon(head(1)),0,length(head))+NMEA.lon(head);
            lon_patch_tail=linspace(0,TRAJ.lon(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.lon(tail(end)),length(tail))+NMEA.lon(tail);
            lon_patch=[lon_patch_head,NMEA.lon(patch(length(head)+1:length(patch)-length(tail))),lon_patch_tail];
            
            lat_patch_head=linspace(TRAJ.lat(traj_time_gaps_idx(ii)+length_plus)-NMEA.lat(head(1)),0,length(head))+NMEA.lat(head);
            lat_patch_tail=linspace(0,TRAJ.lat(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.lat(tail(end)),length(tail))+NMEA.lat(tail);
            lat_patch=[lat_patch_head,NMEA.lat(patch(length(head)+1:length(patch)-length(tail))),lat_patch_tail];
            
            elev_patch_head=linspace(TRAJ.elev(traj_time_gaps_idx(ii)+length_plus)-NMEA.elev(head(1)),0,length(head))+NMEA.elev(head);
            elev_patch_tail=linspace(0,TRAJ.elev(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.elev(tail(end)),length(tail))+NMEA.elev(tail);
            elev_patch=[elev_patch_head,NMEA.elev(patch(length(head)+1:length(patch)-length(tail))),elev_patch_tail];
            
            roll_patch_head=linspace(TRAJ.roll(traj_time_gaps_idx(ii)+length_plus)-NMEA.roll(head(1)),0,length(head))+NMEA.roll(head);
            roll_patch_tail=linspace(0,TRAJ.roll(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.roll(tail(end)),length(tail))+NMEA.roll(tail);
            roll_patch=[roll_patch_head,NMEA.roll(patch(length(head)+1:length(patch)-length(tail))),roll_patch_tail];
            
            pitch_patch_head=linspace(TRAJ.pitch(traj_time_gaps_idx(ii)+length_plus)-NMEA.pitch(head(1)),0,length(head))+NMEA.pitch(head);
            pitch_patch_tail=linspace(0,TRAJ.pitch(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.pitch(tail(end)),length(tail))+NMEA.pitch(tail);
            pitch_patch=[pitch_patch_head,NMEA.pitch(patch(length(head)+1:length(patch)-length(tail))),pitch_patch_tail];
            
            heading_patch_head=linspace(TRAJ.heading(traj_time_gaps_idx(ii)+length_plus)-NMEA.heading(head(1)),0,length(head))+NMEA.heading(head);
            heading_patch_tail=linspace(0,TRAJ.heading(traj_time_gaps_idx(ii)+1+length_plus)-NMEA.heading(tail(end)),length(tail))+NMEA.heading(tail);
            heading_patch=[heading_patch_head,NMEA.heading(patch(length(head)+1:length(patch)-length(tail))),heading_patch_tail];                                 
          else % if the time gap is shorter than 20 secs
            time_patch=NMEA.time(patch);
            time_temp=[NMEA.time(patch-5:patch-1),left_time,NMEA.time(patch),right_time,NMEA.time(patch+1:patch+5)]; 
            
            lon_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.lon(patch-5:patch+5),time_temp);
            lon_interp_left=lon_interp(6);
            lon_interp_right=lon_interp(end-5);                       
            lon_patch=linspace(TRAJ.lon(traj_time_gaps_idx(ii)+length_plus)-lon_interp_left,TRAJ.lon(traj_time_gaps_idx(ii)+1+length_plus)-lon_interp_right,length(time_patch))+NMEA.lon(patch);
            
            lat_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.lat(patch-5:patch+5),time_temp);
            lat_interp_left=lat_interp(6);
            lat_interp_right=lat_interp(end-5);                       
            lat_patch=linspace(TRAJ.lat(traj_time_gaps_idx(ii)+length_plus)-lat_interp_left,TRAJ.lat(traj_time_gaps_idx(ii)+1+length_plus)-lat_interp_right,length(time_patch))+NMEA.lat(patch);
            
            elev_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.elev(patch-5:patch+5),time_temp);
            elev_interp_left=elev_interp(6);
            elev_interp_right=elev_interp(end-5);                       
            elev_patch=linspace(TRAJ.elev(traj_time_gaps_idx(ii)+length_plus)-elev_interp_left,TRAJ.elev(traj_time_gaps_idx(ii)+1+length_plus)-elev_interp_right,length(time_patch))+NMEA.elev(patch);
            
            roll_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.roll(patch-5:patch+5),time_temp);
            roll_interp_left=roll_interp(6);
            roll_interp_right=roll_interp(end-5);                       
            roll_patch=linspace(TRAJ.roll(traj_time_gaps_idx(ii)+length_plus)-roll_interp_left,TRAJ.roll(traj_time_gaps_idx(ii)+1+length_plus)-roll_interp_right,length(time_patch))+NMEA.roll(patch);
            
            pitch_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.pitch(patch-5:patch+5),time_temp);
            pitch_interp_left=pitch_interp(6);
            pitch_interp_right=pitch_interp(end-5);                       
            pitch_patch=linspace(TRAJ.pitch(traj_time_gaps_idx(ii)+length_plus)-pitch_interp_left,TRAJ.pitch(traj_time_gaps_idx(ii)+1+length_plus)-pitch_interp_right,length(time_patch))+NMEA.pitch(patch);
            
            heading_interp=interp1([NMEA.time(patch-5:patch+5)],NMEA.heading(patch-5:patch+5),time_temp);
            heading_interp_left=heading_interp(6);
            heading_interp_right=heading_interp(end-5);                       
            heading_patch=linspace(TRAJ.heading(traj_time_gaps_idx(ii)+length_plus)-heading_interp_left,TRAJ.heading(traj_time_gaps_idx(ii)+1+length_plus)-heading_interp_right,length(time_patch))+NMEA.heading(patch);
          end 
        TRAJ.lon=[TRAJ.lon(1:traj_time_gaps_idx(ii)+length_plus),lon_patch,TRAJ.lon(traj_time_gaps_idx(ii)+1+length_plus:end)];
        TRAJ.lon=mod(TRAJ.lon+180, 360)-180;
        TRAJ.lat=[TRAJ.lat(1:traj_time_gaps_idx(ii)+length_plus),lat_patch,TRAJ.lat(traj_time_gaps_idx(ii)+1+length_plus:end)];
        TRAJ.elev=[TRAJ.elev(1:traj_time_gaps_idx(ii)+length_plus),elev_patch,TRAJ.elev(traj_time_gaps_idx(ii)+1+length_plus:end)];
        TRAJ.roll=[TRAJ.roll(1:traj_time_gaps_idx(ii)+length_plus),roll_patch,TRAJ.roll(traj_time_gaps_idx(ii)+1+length_plus:end)];
        TRAJ.pitch=[TRAJ.pitch(1:traj_time_gaps_idx(ii)+length_plus),pitch_patch,TRAJ.pitch(traj_time_gaps_idx(ii)+1+length_plus:end)];
        TRAJ.heading=[TRAJ.heading(1:traj_time_gaps_idx(ii)+length_plus),heading_patch,TRAJ.heading(traj_time_gaps_idx(ii)+1+length_plus:end)];
        length_plus=length_plus+length(patch);
      end
  end
  FINAL.time=TRAJ.gps_time;%a struct called "FINAL" will be printed to a .csv file at last
  FINAL.elev=TRAJ.elev;
  FINAL.lat=TRAJ.lat;
  FINAL.lon=TRAJ.lon;
  FINAL.roll=TRAJ.roll;
  FINAL.pitch=TRAJ.pitch;
  FINAL.heading=TRAJ.heading;  
else %we only have nmea files   

  FINAL.time=NMEA.time;
  FINAL.elev=NMEA.elev;
  FINAL.lat=NMEA.lat;
  FINAL.lon=NMEA.lon;
  FINAL.roll=NMEA.roll;
  FINAL.pitch=NMEA.pitch;
  FINAL.heading=NMEA.heading;
end
fix_minor_idx=find(diff(FINAL.time)<=1e-6);
FINAL.time(fix_minor_idx+1)=FINAL.time(fix_minor_idx+1)+1e-6;
%% write data and generate a .csv file per day
fid = fopen(full_file_out,'w');%%%%%%%write csv file
% % % % fid = fopen(strcat(out_dir,param1.date,'_nmea'),'w');%%%%%%%write csv file
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t\n',...
'YR',',','Day',',','Time',',','Latitude',',','Longitude',',','Elevation',',','Roll',',','Pitch',',','Heading');

  for n=1:length(FINAL.time)
    fprintf(fid,'%s\t%s\t%f\t%s\t%.6f\t%s\t%s\t%s\t%s\t%s\t%f\t%s\t%f\t%s\t%f\t%s\t%f\t\n',...
     param1.year(param1.year~='0'),',', day,',', FINAL.time(n),',',FINAL.lat(n),',',FINAL.lon(n),',',FINAL.elev(n),',',FINAL.roll(n),',',FINAL.pitch(n),',',FINAL.heading(n));
  
  end

fclose(fid);

fprintf('Successful\n')



function [traj_mark_out, season_name, traj_folder, traj_filename,gps_source] = icards_data_traj_list( day_seg,traj_mark )
% This function is for the convenience to load trajactory files of some
% certain days. Only 1999-2002 has this kind of file as a supplements to
% nmea file
  if traj_mark
    day=day_seg(1:8);
    switch day
      case '19990507'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990507_aa_l12_jgs_itrf96_28jun99';gps_source=strcat('atm-final_',day);
      case '19990510'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990510_aa_l12_jgs_itrf96_09jun99';gps_source=strcat('atm-final_',day);
      case '19990511'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990511_aa_l12_jgs_itrf96_15jun99';gps_source=strcat('atm-final_',day);
      case '19990512'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990512_aa_l12_jgs_itrf96_10jun99';gps_source=strcat('atm-final_',day);
      case '19990513'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990513_aa_l12_jgs_itrf96_15jun99';gps_source=strcat('atm-final_',day); 
      case '19990514'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990514_aa_l12_jgs_itrf96_11jun99';gps_source=strcat('atm-final_',day);
      case '19990517'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990517_aa_l12_jgs_itrf96_09jun99';gps_source=strcat('atm-final_',day);  
      case '19990518'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990518_aa_l12_jgs_itrf96_21jun99';gps_source=strcat('atm-final_',day);
      case '19990519'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990519_aa_l12_jgs_itrf96_11jun99';gps_source=strcat('atm-final_',day);  
      case '19990521'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990521_aa_l12_jgs_itrf96_24jun99';gps_source=strcat('atm-final_',day);
      case '19990523'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990523_aa_l12_jgs_itrf96_25jun99';gps_source=strcat('atm-final_',day);  
      case '19990524'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990524_aa_l12_jgs_itrf96_10jun99';gps_source=strcat('atm-final_',day);
      case '19990525'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='990525_aa_l12_jgs_itrf96_08jun99';gps_source=strcat('atm-final_',day);
      case '20010520'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='010520_aa_l12_jgs_itrf97_09jul01';gps_source=strcat('atm-final_',day);
      case '20010521'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='010521_aa_l12_jgs_itrf97_19jul01';gps_source=strcat('atm-final_',day);
      case '20010523'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='010523_aa_l12_jgs_itrf97_09aug01';gps_source=strcat('atm-final_',day);
      case '20010524'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='010524_aa_l12_jgs_itrf97_06jul01';gps_source=strcat('atm-final_',day);
      case '20010527'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='010527_aa_l12_cfm_itrf97_02aug01';gps_source=strcat('atm-final_',day);
      case '20020518'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020518_aa_l12_jgs_itrf00_05jul02_b898';gps_source=strcat('atm-final_',day); 
      case '20020524'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020524_aa_l12_jgs_itrf00_17jul02_npm';gps_source=strcat('atm-final_',day);
      case '20020528'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020528_aa_l12_jgs_itrf00_16jul02_b898';gps_source=strcat('atm-final_',day);
      case '20020529'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020529_aa_l12_cfm_itrf00_16jul02_thule';gps_source=strcat('atm-final_',day);
      case '20020530'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020530_aa_l12_cfm_itrf00_01jul02_2sta';gps_source=strcat('atm-final_',day);
      case '20020531'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020531_aa_l12_jgs_itrf00_06aug02_6138';gps_source=strcat('atm-final_',day); 
      case '20020601'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020601_aa_l12_jgs_itrf00_05aug02_6138';gps_source=strcat('atm-final_',day);
      case '20020604'
        traj_mark_out=1;season_name=strcat(day(1:4),'_Greenland_P3');traj_folder=strcat(season_name,'_traj');
        traj_filename='020604_aa_l12_jgs_itrf00_09jul02_6138';gps_source=strcat('atm-final_',day);
      otherwise
        traj_mark_out=0; season_name=[]; traj_folder=[]; traj_filename=[]; gps_source=[];
%         fprintf('No trajactory file for %s\n',day);
        
    end
  else
    traj_mark_out=0;season_name=[];traj_folder=[];traj_filename=[];gps_source=[];
  end
    
    
  
end


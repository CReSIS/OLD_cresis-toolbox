function temp_cham_rate_DO_NOT_USE(temp_celsius,celsius_per_min,plot_en)
% PLEASE VISIT https://wiki.cresis.ku.edu/cresis/Temperature_Chamber FOR
% MORE DETAILS
%This program caused a presure switch to trip and the option for rate of
%change was removed. So, Please use temp_cham.m
%
%This is only for reference and please do not use this program.
%>>function temp_cham(temp_celsius,celsius_per_min)        
%>>Sets the Chamber to the temperature specified.(and Rate of change)
%>>Usage: temp_cham(-30)
%         temp_cham(-20,1/3)
         
 if ~exist('plot_en','var')
     plot_en=0;
 end
if(temp_celsius<-30 || temp_celsius>177)  
          fprintf('Temperature out of testing Limits(-30 to 177)\n');
          return;
else
    fprintf('Temperature with in the specified Range\n');
end
switch nargin
    case 1
        fprintf('Required Temperature = %d\n', temp_celsius);
        celsius_per_min=0;
        fprintf('case 1');
    case {2,3}
        fprintf('Required Temperature = %d\n', temp_celsius);
        fprintf('Required Rate = %d\n', celsius_per_min);
        
    otherwise
        fprintf('NO Input Arguments\n');
        fprintf('case no');
        return;
end                     

%% Initializing the chamber
% Set the chamber to INTERFACE: RS 232
                     %MODE: COMMAND
                     %COMMAND END OF LINE: CR/LF  (changeable)
                     %SERIAL CONFIGURATION----> #following
%% CREATE PORT AND CONFIGURE CHAMBER

    s=serial('COM10'); 
        %Define the COM# before interfacing
        %Check the values on the Controller screen OR Set the following default values
    set(s, 'BaudRate', 19200); 
        % Chamber can talk at 310, 600, 1200, 2400, 4800, 9600, 19200
    set(s,'DataBits', 8); 
        % 8 or 7 bits
    set(s,'StopBits', 1); 
        % 1 or 2 bits
    set(s,'Parity', 'none'); 
        % EVEN or ODD or dont care--->none
    set(s,'Terminator', 'CR/LF'); 
        % Check for same End of Command Line. Possibility of Time-out
        % Chamber Delimiter toggles between either of the 3----> CR, LF, CR/LF
    set(s,'FlowControl', 'software'); 
        % Xon/Xoff(software) or RTS/CTS(hard wired)
set(s,'TimeOut',3);
%% CONNECT TO CHAMBER
fopen(s);
fprintf('*****SYSTEM CONNECTED*****');
%% TALK TO CHAMBER

%fprintf(s,'start system');
    fprintf(s,'set lights = true'); % Switch on the chamber lights //NOT operational
    temp=fscanf(s);
    fprintf(s,'read auxs'); %//This program uses no AUXS and events
    auxs=fscanf(s);
    temp=fscanf(s);
    fprintf(s,'read events');
    events=fscanf(s);
    temp=fscanf(s);
%% Set the req value

    fprintf(s,sprintf('Set setpoint 1 = %d',temp_celsius));
    temp=fscanf(s);
    
    fprintf(s,'read setpoint 1');
    set_point=str2double(fscanf(s));
    temp=fscanf(s);
    fprintf('Set Point = %d\n',set_point);

    fprintf(s,'read pv 1');
    pro_var=str2double(fscanf(s)); % Process Variable
    temp=fscanf(s);
    fprintf('Process Variable = %d\n',pro_var);

    fprintf(s,'read deviation 1');
    deviation=str2double(fscanf(s));
    temp=fscanf(s);
    fprintf('Deviation = %d\n',deviation);
    
%% Set the step time ADDED FOR RATE OF CHANGE
if ~(celsius_per_min==0)
    time_s=ceil(abs(deviation)*60/celsius_per_min);
    
    set_hh = floor(time_s / 3600);
    tem = time_s - set_hh * 3600;
    set_mm = floor(tem / 60);
    set_ss = tem - set_mm * 60;
    fprintf('Please enter the step time as %2d:%2d:%2d\n',set_hh,set_mm,set_ss);
    while(1)
        response=input('Enter "y" when it is set(y)\n','s');
        if strcmpi(response,'y')
            break;
         end
    end
    fprintf('Step Time is Set\n');
end
%% start

    fprintf(s,'start system');
    temp=fscanf(s);
    fprintf('*****SYSTEM STARTED\n');
%% LIVE
    
    fprintf('Temperature is being set to --->   %d...\n',temp_celsius);
    fprintf(s,'read pv 1');
    live_temp=str2double(fscanf(s));
    temp=fscanf(s);
    x=live_temp;
    time=0;
    stop_time=0;

  
    while(live_temp~=temp_celsius)  
        fprintf(s,'read pv 1');
        live_temp=str2double(fscanf(s));
        temp=fscanf(s);

        x=[x live_temp];  
        if(plot_en)
        plot(x);
        axis auto;
        grid on;
        hold on;
        end
        time=time+1;
        pause(1);
        if(live_temp<-31)       
             fprintf(s,'stop system');
             temp=fscanf(s);
             stop_time=time;
             fprintf('Temperature Exceeded the Minimum Limits\n');
        end
        if(live_temp>177)       
             fprintf(s,'stop system');
             temp=fscanf(s);
             stop_time=time;
             fprintf('Temperature Exceeded the Maximum Limits\n');
        end
    end
    fprintf('System is at %d in %d seconds\n',live_temp,time);
    
 %% stop

      response=input('Stop the system? y\n','s');
        if strcmpi(response,'y')
            fprintf(s,'stop system');
            temp=fscanf(s);
            fprintf('*****SYSTEM STOPPED*****\n');
        end
         if strcmpi(response,'n')
            fprintf('*****SYSTEM STILL RUNNING*****\n');
         end
     
%% OPTIONAL:
dir_name=datestr(clock,'yyyy_mmm_dd');
if ~exist(dir_name,'dir')
    mkdir(dir_name);
end
save(datestr(clock,'yyyy_mmm_dd\\HH_MM_SS'));

    
%% DISCONNECT AND CLEAN
    
    fclose(s);
    fprintf('*****DISCONNECTED*****');
    delete(s);
    clear s;
end
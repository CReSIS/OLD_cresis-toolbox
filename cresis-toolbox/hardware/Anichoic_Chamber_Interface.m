   function varargout = Anichoic_Chamber_Interface(varargin)
% ANICHOIC_CHAMBER_INTERFACE MATLAB code for Anichoic_Chamber_Interface.fig
%      ANICHOIC_CHAMBER_INTERFACE, by itself, creates a new ANICHOIC_CHAMBER_INTERFACE or raises the existing
%      singleton*.
%
%      H = ANICHOIC_CHAMBER_INTERFACE returns the handle to a new ANICHOIC_CHAMBER_INTERFACE or the handle to
%      the existing singleton*.
%
%      ANICHOIC_CHAMBER_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANICHOIC_CHAMBER_INTERFACE.M with the given input arguments.
%
%      ANICHOIC_CHAMBER_INTERFACE('Property','Value',...) creates a new ANICHOIC_CHAMBER_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Anichoic_Chamber_Interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Anichoic_Chamber_Interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Anichoic_Chamber_Interface

% Last Modified by GUIDE v2.5 17-Mar-2013 09:29:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Anichoic_Chamber_Interface_OpeningFcn, ...
                   'gui_OutputFcn',  @Anichoic_Chamber_Interface_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before Anichoic_Chamber_Interface is made visible.
function Anichoic_Chamber_Interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Anichoic_Chamber_Interface (see VARARGIN)

% Choose default command line output for Anichoic_Chamber_Interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Delete All GPIB interfaces
interfaces=instrfind;
if ~isempty(interfaces)
    fclose(interfaces);
end

% UIWAIT makes Anichoic_Chamber_Interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Anichoic_Chamber_Interface_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in Inst_Board_Index.
function Inst_Board_Index_Callback(hObject, eventdata, handles)
% hObject    handle to Inst_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Inst_Board_Index contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Inst_Board_Index


% --- Executes during object creation, after setting all properties.
function Inst_Board_Index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inst_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Inst_Prime_Address.
function Inst_Prime_Address_Callback(hObject, eventdata, handles)
% hObject    handle to Inst_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Inst_Prime_Address contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Inst_Prime_Address


% --- Executes during object creation, after setting all properties.
function Inst_Prime_Address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Inst_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Turntable_Board_Index.
function Turntable_Board_Index_Callback(hObject, eventdata, handles)
% hObject    handle to Turntable_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Turntable_Board_Index contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Turntable_Board_Index


% --- Executes during object creation, after setting all properties.
function Turntable_Board_Index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Turntable_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Turntable_Prime_Address.
function Turntable_Prime_Address_Callback(hObject, eventdata, handles)
% hObject    handle to Turntable_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Turntable_Prime_Address contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Turntable_Prime_Address


% --- Executes during object creation, after setting all properties.
function Turntable_Prime_Address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Turntable_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Phi_Board_Index.
function Phi_Board_Index_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Phi_Board_Index contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Phi_Board_Index


% --- Executes during object creation, after setting all properties.
function Phi_Board_Index_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi_Board_Index (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Phi_Prime_Address.
function Phi_Prime_Address_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Phi_Prime_Address contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Phi_Prime_Address


% --- Executes during object creation, after setting all properties.
function Phi_Prime_Address_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi_Prime_Address (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Inst_Connect.
function Inst_Connect_Callback(hObject, eventdata, handles)
% hObject    handle to Inst_Connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Inst_board_Index = get(handles.Inst_Board_Index,'Value')-1;
Inst_Prim_Address = get(handles.Inst_Prime_Address,'Value')-1;

try         % Try and create gpib object, update status text on error
    Inst_gpib = gpib('ni',Inst_board_Index,Inst_Prim_Address);
catch err
    if (strcmp(err.identifier,'instrument:gpib:invalidInterface'))
        set(handles.Inst_Status_Text,'String','Invalid Interface')
    end
end

if exist('Inst_gpib')
    try
        fopen(Inst_gpib);    %Open connection to instrument
    catch err
        if (strcmp(err.identifier,'instrument:fopen:opfailed'))
            set(handles.Inst_Status_Text,'String','Unable to connect')
            return
        end
    end
    
    % Check if connection opened correctly
    if get(Inst_gpib, 'Status') ~= 'open' % Error out if not open
        set(handles.Inst_Status_Text,'String','Unable to Connect')
        return
    else % Get instrument name if connected
        fprintf(Inst_gpib,'*IDN?')
        Inst_IDN=fscanf(Inst_gpib);
        if strncmpi(Inst_IDN,'Agilent Technologies,N5230C',27)
            set(handles.Inst_Status_Text,'String','PNA')
        elseif strncmpi(Inst_IDN,'Agilent Technologies, E4446A',28)
            set(handles.Inst_Status_Text,'String','SA')
        elseif strncmpi(Inst_IDN,'TEKTRONIX,DPO70404C',19)
            set(handles.Inst_Status_Text,'String','Oscope')
        else
            set(handles.Inst_Status_Text,'String',Inst_IDN)
        end
        fclose(Inst_gpib)
    end
    
    
end

% --- Executes on button press in Theta_Connect.
function Theta_Connect_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_Connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Theta_board_Index = get(handles.Turntable_Board_Index,'Value')-1;
Theta_Prim_Address = get(handles.Turntable_Prime_Address,'Value')-1;

try
    Theta_gpib = gpib('ni',Theta_board_Index,Theta_Prim_Address);
catch err
    if (strcmp(err.identifier,'instrument:gpib:invalidInterface'))
        set(handles.Theta_Status_Text,'String','Invalid Interface')
    end
end

if exist('Theta_gpib')
    try
        fopen(Theta_gpib);    %Open connection to Theta
    catch err
        if (strcmp(err.identifier,'instrument:fopen:opfailed'))
            set(handles.Theta_Status_Text,'String','Unable to connect')
            return
        end
    end
    
    % Check if connection opened correctly
    if get(Theta_gpib, 'Status') ~= 'open' % Error out if not open
        set(handles.Theta_Status_Text,'String','Unable to Connect')
        return
    else % Get instrument name if connected
        fprintf(Theta_gpib,'*IDN?')
        Theta_IDN=fscanf(Theta_gpib);
        if strncmpi(Theta_IDN,'EMCO,2090-TT',12)
            set(handles.Theta_Status_Text,'String','Connected')
        else
            set(handles.Theta_Status_Text,'String',Theta_IDN)
        end
        fclose(Theta_gpib)
    end
    
end

% --- Executes on button press in Phi_Connect.
function Phi_Connect_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Connect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Phi_board_Index = get(handles.Phi_Board_Index,'Value')-1;
Phi_Prim_Address = get(handles.Phi_Prime_Address,'Value')-1;

try
    Phi_gpib = gpib('ni',Phi_board_Index,Phi_Prim_Address);
catch err
    if (strcmp(err.identifier,'instrument:gpib:invalidInterface'))
        set(handles.Phi_Status_Text,'String','Invalid Interface')
    end
end

if exist('Phi_gpib')
    try
        fopen(Phi_gpib);    %Open connection to Phi
    catch err
        if (strcmp(err.identifier,'instrument:fopen:opfailed'))
            set(handles.Phi_Status_Text,'String','Unable to connect')
            return
        end
    end
    
    % Check if connection opened correctly
    if get(Phi_gpib, 'Status') ~= 'open' % Error out if not open
        set(handles.Phi_Status_Text,'String','Unable to Connect')
        return
    else % Get instrument name if connected
        fprintf(Phi_gpib,'*IDN?')
        Phi_IDN=fscanf(Phi_gpib);
        if strncmpi(Phi_IDN,'EMCO,2090-TT',12)
            set(handles.Phi_Status_Text,'String','Connected')
        else
            set(handles.Phi_Status_Text,'String',Phi_IDN)
        end
        fclose(Phi_gpib)
    end
end

function Phi_Start_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Phi_Start as text
%        str2double(get(hObject,'String')) returns contents of Phi_Start as a double


% --- Executes during object creation, after setting all properties.
function Phi_Start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Phi_Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Phi_Stop as text
%        str2double(get(hObject,'String')) returns contents of Phi_Stop as a double


% --- Executes during object creation, after setting all properties.
function Phi_Stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Phi_Step_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Phi_Step as text
%        str2double(get(hObject,'String')) returns contents of Phi_Step as a double


% --- Executes during object creation, after setting all properties.
function Phi_Step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Phi_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_Start_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_Start as text
%        str2double(get(hObject,'String')) returns contents of Theta_Start as a double


% --- Executes during object creation, after setting all properties.
function Theta_Start_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_Start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_Stop_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_Stop as text
%        str2double(get(hObject,'String')) returns contents of Theta_Stop as a double


% --- Executes during object creation, after setting all properties.
function Theta_Stop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_Stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Theta_Step_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Theta_Step as text
%        str2double(get(hObject,'String')) returns contents of Theta_Step as a double


% --- Executes during object creation, after setting all properties.
function Theta_Step_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Theta_Step (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Output_Path_Text_Callback(hObject, eventdata, handles)
% hObject    handle to Output_Path_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Output_Path_Text as text
%        str2double(get(hObject,'String')) returns contents of Output_Path_Text as a double


% --- Executes during object creation, after setting all properties.
function Output_Path_Text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Output_Path_Text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Output_Path_Button.
function Output_Path_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Output_Path_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[Output_File, Output_Path] = uiputfile;
set(handles.Output_Path_Text,'String',[Output_Path, Output_File])


% --- Executes on button press in Test_Data_Capture_Button.
function Test_Data_Capture_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Test_Data_Capture_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Inst_board_Index = get(handles.Inst_Board_Index,'Value')-1;
Inst_Prim_Address = get(handles.Inst_Prime_Address,'Value')-1;

try         % Try and create gpib object, display error in msgbox on error
    Inst_gpib = gpib('ni',Inst_board_Index,Inst_Prim_Address);
catch err
    msgbox(['Error: ' err.identifier])
end



if exist('Inst_gpib')
    
    Buffer_size = 32001 * 8;    % Max buffer for 32001 8-bit data points
    set(Inst_gpib,'InputBufferSize',Buffer_size)
    
    try
        fopen(Inst_gpib);    %Open connection to instrument
    catch err
        error('Unable to Connect')
        return
    end
    
    % Check if connection opened correctly
    if get(Inst_gpib, 'Status') ~= 'open' % Error out if not open
        msgbox('Unable to Connect')
        return
    else
        % Identify Instrument
        fprintf(Inst_gpib,'*IDN?')
        Inst_IDN=fscanf(Inst_gpib);
        
        % Capture A Single Data Set
            %%%%%%%
            % PNA %
            %%%%%%%
            if strncmpi(Inst_IDN,'Agilent Technologies,N5230C',27)
                fprintf(Inst_gpib,'FORM:BORD SWAP')         % We need reverse GPIB bitorder to read data 
                fprintf(Inst_gpib,'format:data real,64')    % Set Data type to double presission
                                                            % (this doesn't  really increase the
                                                            % presission so much as speed up data
                                                            % transfer.)
                
                fprintf(Inst_gpib,'CALC:PAR:CAT?')
                Meas_Name=fscanf(Inst_gpib);
                [matchstart,matchend]=regexpi(Meas_Name,'(.)+,');
                fprintf(Inst_gpib,'CALC:PAR:SEL ''%s''',Meas_Name(matchstart+1:matchend-1))
                
                % Get Data
                fprintf(Inst_gpib,'SENS:SWE:POIN?');
                num_of_points=str2num(fscanf(Inst_gpib));
                
                fprintf(Inst_gpib,'CALC:DATA? fdata');

                start_Data_Char = fscanf(Inst_gpib,'%c',1);     % Data block always starts with a
                                                                % '#' char that we need to clear.
                byte_count_size = fscanf(Inst_gpib,'%c',1);     % Found out how many characters the
                                                                % byte count is.
                byte_count = fscanf(Inst_gpib,'%c',str2num(byte_count_size));   % Clear out the 
                                                                                % byte size block

                Test_Data = fread(Inst_gpib,num_of_points,'double');    % Read Data
              
                % Get X-Axis Values
                fprintf(Inst_gpib,'CALC:X?');

                start_Data_Char = fscanf(Inst_gpib,'%c',1);     % Data block always starts with a
                                                                % '#' char that we need to clear.
                byte_count_size = fscanf(Inst_gpib,'%c',1);     % Found out how many characters the
                                                                % byte count is.
                byte_count = fscanf(Inst_gpib,'%c',str2num(byte_count_size));   % Clear out the 
                                                                                % byte size block

                Test_Data_X_Axis = fread(Inst_gpib,num_of_points,'double');    % Read Data
                
            end
            %%%%%%%%%%%%%%%%%%%%%
            % Spectrum Analyzer %
            %%%%%%%%%%%%%%%%%%%%%
            if strncmpi(Inst_IDN,'Agilent Technologies, E4446A',28)
                
            end
            %%%%%%%%%%%
            % O-Scope %
            %%%%%%%%%%%
            if strncmpi(Inst_IDN,'TEKTRONIX,DPO70404C',19)     
                %Acquires a single acquisition
                fprintf(Inst_gpib,'ACQ:STATE ON');
                pause(0.2)

                
                % Set Channel 1 as source
                fprintf(Inst_gpib,'DATA:SOURCE CH1');
                pause(0.01)
                fprintf(Inst_gpib,'HORizontal:ACQLENGTH?')
                % Setbuffer size
                numSamples=str2num(fscanf(Inst_gpib));
                fprintf(Inst_gpib,'WFMOutpre:BYT_Nr?')
                bytes_per_sample = str2num(fscanf(Inst_gpib));
                fclose(Inst_gpib);
                Buffer_size = bytes_per_Sample*numSamples;
                set(Inst_gpib,'InputBufferSize',Buffer_size)
                fopen(Inst_gpib)
                %Configure byte order
                fprintf(Inst_gpib,'WFMOutpre:BYT_Or LSB')
                % Setup Data Transfer to Send the Entire Record Length
                fprintf(Inst_gpib,'Data:Start 1')
                fprintf(Inst_gpib,['Data:Stop ' num2str(numSamples)])
                % Request data
                fwrite(Inst_gpib,'CURVe?');
                pause(0.01)
                % Clear first bit
                fread(Inst_gpib,1);
                % Read length of header
                a = char(fread(Inst_gpib,1));
                % Read data length
                bytes = char(fread(Inst_gpib,str2double(a))');
                % Read Test data
                Test_Data = fread(Inst_gpib,numSamples,['int' num2str(bytes_per_sample*8)]);
                % Get Multiplication Factor
                fprintf(Inst_gpib,'WFMOUTPRE:YMULT?');
                pause(0.1);
                mult=str2num(fscanf(Inst_gpib));
                Test_Data = mult*Test_Data;
                
                % Generate Test Data axis
                fprintf(Inst_gpib,'HOR:ACQDURATION?');
                Acq_duration = str2num(fscanf(Inst_gpib));
                Test_Data_X_Axis = linspace(0,Acq_duration,numSamples)';

            end
    end
            
        % Display Data
        plot(Test_Data_X_Axis,Test_Data)
        
        %Cleanup
        fclose(Inst_gpib)
end

% --- Executes on button press in Start_Test_Button.
function Start_Test_Button_Callback(hObject, eventdata, handles)
% hObject    handle to Start_Test_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check For Required Fields
    % Output Path
    Output_Path_Text = get(handles.Output_Path_Text,'String');
    if isempty(Output_Path_Text)
        msgbox('Output Path Not Set')
        return
    end
    
    % At least 1 active axis
    Theta_Active = get(handles.Theta_Active,'Value');
    Phi_Active = get(handles.Phi_Active,'Value');
    if ~(Theta_Active || Phi_Active)
        msgbox('No Active Axis Set')
        return
    end
        
    % Valid active axis start/stop/step values
    Theta_Start = str2num(get(handles.Theta_Start,'String'));
    Theta_Stop = str2num(get(handles.Theta_Stop,'String'));
    Theta_Step = str2num(get(handles.Theta_Step,'String'));
    
    Phi_Start = str2num(get(handles.Phi_Start,'String'));
    Phi_Stop = str2num(get(handles.Phi_Stop,'String'));
    Phi_Step = str2num(get(handles.Phi_Step,'String'));
   
    if Theta_Active
        if (Theta_Stop < Theta_Start) || Theta_Step < 0
            msgbox('Invalid Theta Settings')
            return
        end
    end
    
    if Phi_Active
        if (Phi_Stop < Phi_Start) || Phi_Step < 0
            msgbox('Invalid Phi Settings')
            return
        end
    end
    
% Check connection status of required devices
    Inst_board_Index = get(handles.Inst_Board_Index,'Value')-1;
    Inst_Prim_Address = get(handles.Inst_Prime_Address,'Value')-1;
    Phi_board_Index = get(handles.Phi_Board_Index,'Value')-1;
    Phi_Prim_Address = get(handles.Phi_Prime_Address,'Value')-1;
    Theta_board_Index = get(handles.Turntable_Board_Index,'Value')-1;
    Theta_Prim_Address = get(handles.Turntable_Prime_Address,'Value')-1;


    try        
        Inst_gpib = gpib('ni',Inst_board_Index,Inst_Prim_Address);
    catch err
        msgbox(['Error: ' err.identifier])
    end

    if Phi_Active
        try
            Phi_gpib = gpib('ni',Phi_board_Index,Phi_Prim_Address);
        catch err
            msgbox(['Error: ' err.identifier])
        end
    end

    if Theta_Active
        try
            Theta_gpib = gpib('ni',Theta_board_Index,Theta_Prim_Address);
        catch err
            msgbox(['Error: ' err.identifier])
        end
    end

% Create Measurement Position List
    if Theta_Active
        Theta_Position_List = [Theta_Start:Theta_Step:Theta_Stop];
    else
        Theta_Position_List = Theta_Start;
    end
    
    if Phi_Active
        Phi_Position_List = [Phi_Start:Phi_Step:Phi_Stop];
    else
        Phi_Position_List = Phi_Start;
    end
    
% Identify Instrument
    if exist('Inst_gpib')
        fopen(Inst_gpib)
        % Check if connection opened correctly
        if get(Inst_gpib, 'Status') ~= 'open' % Error out if not open
            msgbox('Unable to Connect')
            return
        else % Get instrument name if connected
            fprintf(Inst_gpib,'*IDN?')
            Inst_IDN=fscanf(Inst_gpib);
            if strncmpi(Inst_IDN,'Agilent Technologies,N5230C',27)
                Test_Inst = 'PNA';
            elseif strncmpi(Inst_IDN,'Agilent Technologies, E4446A',28)
                Test_Inst = 'SA';
            elseif strncmpi(Inst_IDN,'TEKTRONIX,DPO70404C',19)
                Test_Inst = 'Oscope';
            else
                msgbox('Unknown Instrument')
                return
            end
            fclose(Inst_gpib)
        end
    end

    
% Measurement Loop
    switch Test_Inst
        case 'PNA'
            % Open Connections
                Buffer_size = 2*32001 * 8;    % Max buffer for 32001 8-bit data points
                set(Inst_gpib,'InputBufferSize',Buffer_size)
                fopen(Inst_gpib)

                if Theta_Active
                    fopen(Theta_gpib)
                end

                if Phi_Active
                    fopen(Phi_gpib)
                end
            
            % Move Axis to inital position
                if Theta_Active
                    fprintf(Theta_gpib,'N2')
                    fprintf(Theta_gpib,'SK %3.1f',Theta_Position_List(1))
                end

                if Phi_Active
                    fprintf(Phi_gpib,'N2')
                    fprintf(Phi_gpib,'SK %3.1f',Phi_Position_List(1))
                end
            
            % Setup PNA Data Transfer Parameters
                
                fprintf(Inst_gpib,'FORM:BORD SWAP')         % We need reverse GPIB bitorder to read data 
                fprintf(Inst_gpib,'format:data real,64')    % Set Data type to double presission
                                                            % (this doesn't  really increase the
                                                            % presission so much as speed up data
                                                            % transfer.)
                fprintf(Inst_gpib,'SENS:SWE:POIN?');
                num_of_points=str2num(fscanf(Inst_gpib));
                
            % Get X-Axis Scale
                fprintf(Inst_gpib,'CALC:PAR:CAT?')
                Meas_Name=fscanf(Inst_gpib);
                [matchstart,matchend]=regexpi(Meas_Name,'(.)+,');
                fprintf(Inst_gpib,'CALC:PAR:SEL ''%s''',Meas_Name(matchstart+1:matchend-1))
                
                fprintf(Inst_gpib,'CALC:X?');

                start_Data_Char = fscanf(Inst_gpib,'%c',1);     % Data block always starts with a
                                                                % '#' char that we need to clear.
                byte_count_size = fscanf(Inst_gpib,'%c',1);     % Found out how many characters the
                                                                % byte count is.
                byte_count = fscanf(Inst_gpib,'%c',str2num(byte_count_size));   % Clear out the 
                                                                                % byte size block

                Test_Data_X_Axis = fread(Inst_gpib,num_of_points,'double');    % Read Axis Data
                
                % Check Measurement Format and expect extra points for polar
                % and smith formats
                    fprintf(Inst_gpib,'CALC:FORM?');
                    Data_Format = fscanf(Inst_gpib);
                fprintf(Inst_gpib,'SENS:SWE:POIN?');
                num_of_points=str2num(fscanf(Inst_gpib));
                % In Polar and Smith mode we get 2 data points for every
                % frequency instead of just 1.
                if( strcmp(Data_Format(1:3),'POL') || strcmp(Data_Format(1:3),'SMI') )
                   Extra_Data = 2;
                else
                    Extra_Data =1;
                end
                
            % Setup Primary & Secondary Position Lists
                Primary_Axis_Theata = get(handles.Primary_Axis_Theta,'Value');
                
                if Primary_Axis_Theata
                    Primary_Axis = 'Theta';
                    Primary_Axis_Position_List = Theta_Position_List;
                    Secondary_Axis_Position_List = Phi_Position_List;
                    Primary_Axis_gpib = Theta_gpib;
                    if Phi_Active
                        Secondary_Axis_gpib = Phi_gpib;
                    end
                    
                else
                    Primary_Axis = 'Phi';
                    Primary_Axis_Position_List = Phi_Position_List;
                    Secondary_Axis_Position_List = Theta_Position_List;
                    Primary_Axis_gpib = Phi_gpib;
                    if Phi_Active
                        Secondary_Axis_gpib = Theta_gpib;
                    end
                end
                
             % Setup Data Variables
             num_of_points = Extra_Data*num_of_points;
%              Test_Data = zeros(num_of_points,length(Primary_Axis_Position_List),length(Secondary_Axis_Position_List));
             Primary_Position_Index = repmat(Primary_Axis_Position_List',1,length(Secondary_Axis_Position_List));
             Secondary_Position_Index = repmat(Secondary_Axis_Position_List,length(Primary_Axis_Position_List),1);

             % Setup Output File
            Output_File_Path = get(handles.Output_Path_Text,'String');
            warning('off','all')
            delete(Output_File_Path)
            warning('on','all')
            Output_File_Obj = matfile(Output_File_Path,'Writable',true);
            Output_File_Obj.Test_Data(1:num_of_points,length(Primary_Axis_Position_List),length(Secondary_Axis_Position_List)) = zeros(num_of_points,1,1);

             % Primary Axis Loop
                h = waitbar(0,'Test In Progress','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
                setappdata(h,'canceling',0)
                for i = 1:length(Primary_Axis_Position_List)
                    if getappdata(h,'canceling')
                            break;
                    end
                    
                    waitbar(i/length(Primary_Axis_Position_List))
                    % Position Axis
                    fprintf(Primary_Axis_gpib,'SK %3.1f',Primary_Axis_Position_List(i));
                    in_motion = true;
                    pause(2);
                    
                    % Wait for movement to stop
                    while in_motion
                        if getappdata(h,'canceling')
                            break;
                        end
                        fprintf(Primary_Axis_gpib,'CP?')
                        current_location = str2num(fscanf(Primary_Axis_gpib));
                        
                        if (current_location < Primary_Axis_Position_List(i)+0.9) && (current_location > Primary_Axis_Position_List(i)-0.9)
                            in_motion = false;
                        else
                            pause(1)
                        end
                    end
                    
                    % Secondary Axis Loop
                    for j = 1:length(Secondary_Axis_Position_List)
                        if getappdata(h,'canceling')
                            break;
                        end
                        % Position Axis
                        if length(Secondary_Axis_Position_List) > 1;
                            fprintf(Secondary_Axis_gpib,'SK %3.1f',Secondary_Axis_Position_List(j));
                            in_motion = true;
                            pause(2);

                            % Wait for movement to stop
                            while in_motion
                                if getappdata(h,'canceling')
                                    break;
                                end
                                fprintf(Secondary_Axis_gpib,'CP?')
                                current_location = str2num(fscanf(Secondary_Axis_gpib));

                                if (current_location < Secondary_Axis_Position_List(j)+0.9) && (current_location > Secondary_Axis_Position_List(j)-0.9)
                                    in_motion = false;
                                else
                                    pause(1)
                                end
                            end
                        end
                        
                        % Take Measurement
                        fprintf(Inst_gpib,'CALC:DATA? fdata');

                        start_Data_Char = fscanf(Inst_gpib,'%c',1);     % Data block always starts with a
                                                                        % '#' char that we need to clear.
                        byte_count_size = fscanf(Inst_gpib,'%c',1);     % Found out how many characters the
                                                                        % byte count is.
                        byte_count = fscanf(Inst_gpib,'%c',str2num(byte_count_size));   % Clear out the 
                                                                                        % byte size block

                        Test_Data = fread(Inst_gpib,num_of_points,'double');    % Read Data
                        if length(Secondary_Axis_Position_List) == 1
                            Output_File_Obj.Test_Data(1:num_of_points,i) = Test_Data;
                        else
                            Output_File_Obj.Test_Data(1:num_of_points,i,j) = Test_Data;
                        end
                    end
                end
            delete(h)
            % Save Data File
            Output_File_Obj.Test_Data_X_Axis = Test_Data_X_Axis;
            Output_File_Obj.Primary_Position_Index = Primary_Position_Index;
            Output_File_Obj.Secondary_Position_Index = Secondary_Position_Index;
            Output_File_Obj.Primary_Axis = Primary_Axis;

            % Cleanup
                
                fclose(Inst_gpib)
                
                 if Theta_Active
                    fprintf(Theta_gpib,' SK 0')
                    fclose(Theta_gpib)
                end

                if Phi_Active
                    fprintf(Phi_gpib,' SK 0')
                    fclose(Phi_gpib)
                end
        case 'SA'
            
        case 'Oscope'
            % Open Connections
                fopen(Inst_gpib)
                % Set buffer size
                fprintf(Inst_gpib,'HORizontal:ACQLENGTH?')  % Ask scope for buffer size
                numSamples=str2num(fscanf(Inst_gpib));
                fprintf(Inst_gpib,'WFMOutpre:BYT_Nr?')
                bytes_per_sample = str2num(fscanf(Inst_gpib));
                fclose(Inst_gpib);                          % Close connection and  reconfigure input buffer
                Buffer_size = bytes_per_sample*numSamples;
                set(Inst_gpib,'InputBufferSize',Buffer_size)
                fopen(Inst_gpib)                            % Reopen Inst connection
                
                % Configure byte order
                fprintf(Inst_gpib,'WFMOutpre:BYT_Or LSB')  
                
                % Setup Data Transfer to Send the Entire Record Length
                fprintf(Inst_gpib,'Data:Start 1');
                fprintf(Inst_gpib,['Data:Stop ' num2str(numSamples)])
                
                % Setup
                if Theta_Active
                    fopen(Theta_gpib)
                end

                if Phi_Active
                    fopen(Phi_gpib)
                end
            
            % Move Axis to inital position
                if Theta_Active
                    fprintf(Theta_gpib,'SK %3.1f',Theta_Position_List(1))
                end

                if Phi_Active
                    fprintf(Phi_gpib,'SK %3.1f',Phi_Position_List(1))
                end
            
          % Get X-Axis Scale
              fprintf(Inst_gpib,'HOR:ACQDURATION?');
              Acq_duration = str2num(fscanf(Inst_gpib));
              Test_Data_X_Axis = linspace(0,Acq_duration,numSamples)';
          
          % Setup Primary & Secondary Position Lists
                Primary_Axis_Theata = get(handles.Primary_Axis_Theta,'Value');
                
                if Primary_Axis_Theata
                    Primary_Axis = 'Theta';
                    Primary_Axis_Position_List = Theta_Position_List;
                    Secondary_Axis_Position_List = Phi_Position_List;
                    Primary_Axis_gpib = Theta_gpib;
                    if Phi_Active
                        Secondary_Axis_gpib = Phi_gpib;
                    end
                    
                else
                    Primary_Axis = 'Phi';
                    Primary_Axis_Position_List = Phi_Position_List;
                    Secondary_Axis_Position_List = Theta_Position_List;
                    Primary_Axis_gpib = Phi_gpib;
                    if Phi_Active
                        Secondary_Axis_gpib = Theta_gpib;
                    end
                end
                
        % Setup Data Variables
%             Test_Data = zeros(numSamples,length(Primary_Axis_Position_List),length(Secondary_Axis_Position_List),'int8');
            Primary_Position_Index = repmat(Primary_Axis_Position_List',1,length(Secondary_Axis_Position_List));
            Secondary_Position_Index = repmat(Secondary_Axis_Position_List,length(Primary_Axis_Position_List),1);
            
        % Setup Output File
            Output_File_Path = get(handles.Output_Path_Text,'String');
            warning('off','all')
            delete(Output_File_Path)
            warning('on','all')
            Output_File_Obj = matfile(Output_File_Path,'Writable',true);
            Output_File_Obj.Test_Data(1:numSamples,length(Primary_Axis_Position_List),length(Secondary_Axis_Position_List)) = zeros(numSamples,1,1);
            
        % Get data muliplication factor
            fprintf(Inst_gpib,'WFMOUTPRE:YMULT?');
            Y_Axis_Multiplication_Factor=str2num(fscanf(Inst_gpib));
            
        % Primary Axis Loop
            h = waitbar(0,'Test In Progress','CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0)
            for i = 1:length(Primary_Axis_Position_List)
                if getappdata(h,'canceling')
                        break;
                end

                waitbar(i/length(Primary_Axis_Position_List))
                % Position Axis
                fprintf(Primary_Axis_gpib,'SK %3.1f',Primary_Axis_Position_List(i));
                in_motion = true;
                pause(2);

                % Wait for movement to stop
                while in_motion
                    if getappdata(h,'canceling')
                        break;
                    end
                    fprintf(Primary_Axis_gpib,'CP?')
                    current_location = str2num(fscanf(Primary_Axis_gpib));

                    if (current_location < Primary_Axis_Position_List(i)+0.9) && (current_location > Primary_Axis_Position_List(i)-0.9)
                        in_motion = false;
                    else
                        pause(1)
                    end
                end

                % Secondary Axis Loop
                for j = 1:length(Secondary_Axis_Position_List)
                    if getappdata(h,'canceling')
                        break;
                    end
                    % Position Axis
                    if length(Secondary_Axis_Position_List) > 1;
                        fprintf(Secondary_Axis_gpib,'SK %3.1f',Secondary_Axis_Position_List(j));
                        in_motion = true;
                        pause(2);

                        % Wait for movement to stop
                        while in_motion
                            if getappdata(h,'canceling')
                                break;
                            end
                            fprintf(Secondary_Axis_gpib,'CP?')
                            current_location = str2num(fscanf(Secondary_Axis_gpib));

                            if (current_location < Secondary_Axis_Position_List(j)+0.9) && (current_location > Secondary_Axis_Position_List(j)-0.9)
                                in_motion = false;
                            else
                                pause(1)
                            end
                        end
                    end

                    % Take Measurement
                        fwrite(Inst_gpib,'CURVe?');
                        pause(0.01)
                        % Clear first bit
                        fread(Inst_gpib,1);
                        % Read length of header
                        a = char(fread(Inst_gpib,1));
                        % Read data length
                        bytes = str2num(char(fread(Inst_gpib,str2double(a))'));
                        % Read Test data
                        Test_Data = fread(Inst_gpib,numSamples,['int' num2str(bytes_per_sample*8)]);
                        if length(Secondary_Axis_Position_List) == 1
                            Output_File_Obj.Test_Data(1:numSamples,i) = Test_Data;
                        else
                            Output_File_Obj.Test_Data(1:numSamples,i,j) = Test_Data;
                        end
%                       
                end
            end
        delete(h)
        % Save Data File
        fprintf(Inst_gpib,'WFMOutpre?')
        Waveform_Preamble = fscanf(Inst_gpib);
        
        Output_File_Obj.Waveform_Preamble = Waveform_Preamble;
        Output_File_Obj.Test_Data_X_Axis = Test_Data_X_Axis;
        Output_File_Obj.Primary_Position_Index = Primary_Position_Index;
        Output_File_Obj.Secondary_Position_Index = Secondary_Position_Index;
        Output_File_Obj.Primary_Axis = Primary_Axis;
        Output_File_Obj.Y_Axis_Multiplication_Factor = Y_Axis_Multiplication_Factor;
        
%             % Email Notify Brian of Test Status
            myaddress = 'rsl.emc.lab@gmail.com';
            mypassword = 'rslemclab9wrd';

            setpref('Internet','E_mail',myaddress);
            setpref('Internet','SMTP_Server','smtp.gmail.com');
            setpref('Internet','SMTP_Username',myaddress);
            setpref('Internet','SMTP_Password',mypassword);

            props = java.lang.System.getProperties;
            props.setProperty('mail.smtp.auth','true');
            props.setProperty('mail.smtp.socketFactory.class', ...
                              'javax.net.ssl.SSLSocketFactory');
            props.setProperty('mail.smtp.socketFactory.port','465');
%             sendmail('bdcordill@gmail.com','Anachoic Chamber Test Complete','Chamber measurements are complete')

        % Cleanup

            fclose(Inst_gpib)

             if Theta_Active
                fprintf(Theta_gpib,' SK 0')
                fclose(Theta_gpib)
            end

            if Phi_Active
                fprintf(Phi_gpib,' SK 0')
                fclose(Phi_gpib)
            end

    end



% --- Executes on button press in Theta_Active.
function Theta_Active_Callback(hObject, eventdata, handles)
% hObject    handle to Theta_Active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Theta_Active



% --- Executes on button press in Phi_Active.
function Phi_Active_Callback(hObject, eventdata, handles)
% hObject    handle to Phi_Active (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Phi_Active

function Waitbar_Cancel_Callback()
    fclose(Inst_gpib)
    
    if Theta_Active
        fprintf(Theta_gpib,' SK 0')
        fclose(Theta_gpib)
    end

    if Phi_Active
        fprintf(Phi_gpib,' SK 0')
        fclose(Phi_gpib)
    end
    
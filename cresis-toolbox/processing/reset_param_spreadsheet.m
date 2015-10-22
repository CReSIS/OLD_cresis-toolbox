function reset_param_spreadsheet(param_fn)
%% Opens a the specified Param sheet (param_fn) with columns:
%   Process
%   Vectors
%   Records
%   Frames
%   Get_Heights
%   Set to blank.
%
% INPUT:
%   param_fn:   String, absolute filepath + filename to MS Excel 
%               params spreadsheet(.xls).  
%
% OUTPUT: 
%   
%
% EXAMPLE:
%   param_fn = 'H:\scripts\params\mcords_param_2010_Greenland_DC8.xls';
%   
%   reset_param_spreadsheet(param_fn);
%   
%   
%%
% Create connection to Excel
Excel = actxserver('Excel.Application');

% Make the Excel window visible
Excel.visible = 1;

% Select sheet to operate on
Workbook = Excel.Workbooks.Open(param_fn);
Command_sheet = Workbook.Worksheets.Item('command');

% Specify range
range = ['C6',':','G1000'];

% Set specified range to ''
Command_sheet.Range(range).Value = '';

    
    
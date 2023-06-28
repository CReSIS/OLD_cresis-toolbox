function [] = cmd_line(dateorstruct,seg_id,mission_name,notes)

%Takes the inputs and prints a line of the cmd worksheet
    %date needs to be a string following format YYYYMMDD
    %seg_id should be an integer
    %mission_name  Name of this flight. 
        %For OIB, each day there is a single flight with a specific name 
        %and this column matches that. For other missions, without specific
        %names, a basic regional description should be used 
            %(e.g. Byrd, Jakobshavn, Helheim, etc.)
        %Don't include as field in structure if unused.
    %notes This field must be text and is used for notes about the 
        %corresponding segment. If the segment should not be processed the 
        %phrase "do not process" should be contained in the notes
        %browse_ni called this settings_fn
        %Don't include as field in structure if unused.
if nargin < 2
    S = dateorstruct;
    date = S.date;
    seg_id = S.seg_id;
    if ~isfield(S,'mission_name')
        mission_name = '';
    else
        mission_name = S.mission_name;
    end
    if ~isfield(S,'notes')
        notes = '';
    else
        notes = S.notes;
    end
elseif nargin < 4
    date = dateorstruct;
    notes = '';
    if nargin < 3
        mission_name = '';
    end
else
    date = dateorstruct;
end
fprintf('%s\t%02d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        date,seg_id,'','','','','','','',mission_name,notes);
end
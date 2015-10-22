classdef (HandleCompatible = true) undo_stack < handle
% undo_stack class
%
% A generic class for handling undo across multiple documents (i.e. each
% document can hold overlapping information so when the overlapping
% information in one document changes it needs to change in all the
% documents that are looking at the same information.
%
% undo_stack(param)
% param.id = unique id of this stack stored in "unique_id" property
%   (not used by undo_stack, helpful to uniquely identifying your stack)
%   Can be of any type.
  
  properties
    % unique_id: unique identifier for this stack, can be any type, set
    % when stack is created (not used by undo_stack)
    unique_id
    % doc_list: cell vector list of documents, use attach_document and
    % remove_document to add documents to the list (otherwise not used by
    % undo_stack)
    doc_list
    % stack: cell vector of commands
    stack
    % pointer: location in the stack that represents the current state of
    % the documents (required since we allow "redo" operation)
    pointer
  end
  
  properties (SetAccess = private, GetAccess = private)
    last_pointer
  end
  
  events
    synchronize_event
    save_event
  end
  
  methods
    function obj = undo_stack(param)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.unique_id = param.id;
      obj.doc_list = {};
      obj.stack = {};
      obj.pointer = 0;
      obj.last_pointer = NaN;
    end
    
    % Attach document to undo stack and return the set of commands that
    % need to be run to synchronize that document with the stack
    cmds_list = attach_document(obj, h_doc);
    
    % Remove document from undo stack
    remove_document(obj, h_doc);

    % Get the list set of commands that need to be executed because a
    % synchronize event has occurred.
    [cmds_list,cmds_direction] = get_synchronize_cmds(obj);

    % Get save commands
    cmds_list = get_save_cmds(obj,remove_cmds);

    % Push a set of commands to the undo stack (causes synchronize event)
    push(obj,cmds);
    
    % Causes a set of popped commands to be repushed (causes synchronize event)
    redo(obj);
    
    % Pop a set of commands off the undo stack (this does not cause
    %   these commands to be deleted permanently until another push is
    %   done. Until that happens, they can be repushed with redo.
    %   This command is equivalent to undo. (causes synchronize event)
    pop(obj,cmds);
    
    % Peak at current command (only returns the command, but does not move
    % the pointer)
    cmd = peak(obj);
    
    % Cause all commands that have been pushed and not popped to be
    % committed (i.e. cannot be popped). (causes save event)
    save(obj);

    % Returns logical of whether or not any commands have been pushed
    % and not popped (i.e. current document state is modified).
    modified = ismodified(obj);
      
  end
  
end


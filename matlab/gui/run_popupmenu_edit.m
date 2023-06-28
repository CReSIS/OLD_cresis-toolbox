% script run_popupmenu_edit
%
% Script for testing popupmenu_edit class
%
% See also: popupmenu_edit, run_popupmenu_edit
%
% Author: John Paden

% For debugging purposes, it is nice to have the handle (gpopupmenu_edit) available
global gpopupmenu_edit;
delete(gpopupmenu_edit)
h_fig = figure(1);
% gpopupmenu_edit = popupmenu_edit(h_fig,{'layerData','CSARP_post/layerData'},2);
% gpopupmenu_edit = popupmenu_edit(h_fig,{'layerData','CSARP_post/layerData'});
% gpopupmenu_edit = popupmenu_edit(h_fig,{'layerData','CSARP_post/layerData'},[]);
% gpopupmenu_edit = popupmenu_edit(h_fig,{'layerData','CSARP_post/layerData'},23);
% gpopupmenu_edit = popupmenu_edit(h_fig,{'layerData','CSARP_post/layerData'},NaN);
% gpopupmenu_edit = popupmenu_edit(h_fig,{},NaN);
gpopupmenu_edit = popupmenu_edit(h_fig,{});
% gpopupmenu_edit = popupmenu_edit(h_fig);

return

% Test commands

gpopupmenu_edit.get_value
gpopupmenu_edit.get_selected_string
gpopupmenu_edit.get_list

gpopupmenu_edit.set_value(2)

gpopupmenu_edit.set_selected_string('1')
gpopupmenu_edit.set_selected_string('4')
gpopupmenu_edit.set_list({'layerData','CSARP_post/layerData'})

gpopupmenu_edit.set_enable(false)
gpopupmenu_edit.set_enable(true)

function modePM_callback(obj,hObj,event)

% modePM_callback(obj,hObj,event)
%
% Called when the snake tool's mode select menu has its value changed (from
% basic snake to crandall mode or vice versa). Changes the size of the
% window and hides options that are inapplicable to the current mode.
%

tool_idx = get(obj.top_panel.tool_PM,'Value');
% 1 is basic snake
% 2 is crandall pick
% 3 isn't used

switch tool_idx
  case 3 % not implemented
    return;
  case 1 % basic snake -- hide crandall parameters
    if tool_idx ~= obj.cur_mode
      % backup tool values
      obj.in_rng_sv = str2double(get(obj.top_panel.insert_range_TE, 'String'));
      obj.sn_rng_sv = str2double(get(obj.top_panel.snake_range_TE, 'String'));
      obj.top_sm_sv = get(obj.bottom_panel.topSmoothSlider, 'Value');
      obj.top_pk_sv = get(obj.bottom_panel.topPeakRatioSlider, 'Value');
      obj.bot_sm_sv = get(obj.bottom_panel.bottomSmoothSlider, 'Value');
      obj.bot_pk_sv = get(obj.bottom_panel.bottomPeakRatioSlider, 'Value');
      obj.rep_sv = get(obj.bottom_panel.repulseSlider, 'Value');
      set(obj.h_fig,'Units','Pixels');
      % get position info
      cur_pos = get(obj.h_fig,'Position');
      % set position (window height is basic snake height)
      obj.create_ui_basic(cur_pos(1),cur_pos(2)-(obj.h - cur_pos(4)));
    end
  case 2 % crandall picker -- show crandall parameters
    if tool_idx ~= obj.cur_mode
      % backup tool values
      obj.in_rng_sv = str2double(get(obj.top_panel.insert_range_TE, 'String'));
      obj.sn_rng_sv = str2double(get(obj.top_panel.snake_range_TE, 'String'));
      set(obj.h_fig,'Units','Pixels');
      % get position info
      cur_pos = get(obj.h_fig,'Position');
      % set position (window height is basic snake height)
      obj.create_ui_crandall(cur_pos(1),cur_pos(2)-(obj.crandall_h - cur_pos(4)));
    end
end

obj.cur_mode = tool_idx;

return;
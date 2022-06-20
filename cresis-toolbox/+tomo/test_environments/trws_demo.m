global Nt Nsv Nx step;

close all;
h_fig = figure;
set(h_fig,'WindowKeyPressFcn',@keyboard_input);

Nt  = 10;
Nsv = 10;
Nx  = 10;

step = 45;

draw();

% left, backward, right, forward
% left/right along-track Nx
% backward/forward cross-track Nsv
function draw()
  global step Nx Nsv Nt;
  clf;
  hold on;
  dir = '';
  dir_ind = mod(step, 1);
  switch dir_ind
    case 0
      dir = 'Right';
    case .25
      dir = 'Backward';
    case .5
      dir = 'Left';
    case .75
      dir = 'Forward';
  end
  txt = sprintf('%d, %s', floor(step), dir);
  title(txt);
 
  for w_idx = 1:Nx
      for h_idx = 1:Nsv
          color = [.7 .7 .7];
          cur_col = floor(step);
          [cur_x, cur_sv] = coords(cur_col);
          iter_col = (w_idx-1)*Nsv + h_idx;
          [iter_x, iter_sv] = coords(iter_col);
          
          switch mod(step, 1)
            case 0 % left
              if iter_sv == cur_sv && iter_x == cur_x - 1
                color = [0 1 0];
              end
              
              if iter_x == cur_x && iter_sv == cur_sv - 1
                color = [0 0 1];
              end
              if iter_sv == cur_sv && iter_x == cur_x + 1
                color = [0 0 1];
              end
              if iter_x == cur_x && iter_sv == cur_sv + 1
                color = [0 0 1];
              end
            case .25 % backward
              if iter_x == cur_x && iter_sv == cur_sv - 1
                color = [0 1 0];
              end
              
              if iter_sv == cur_sv && iter_x == cur_x - 1
                color = [0 0 1];
              end
              if iter_sv == cur_sv && iter_x == cur_x + 1
                color = [0 0 1];
              end
              if iter_x == cur_x && iter_sv == cur_sv + 1
                color = [0 0 1];
              end
            case .5 % right
              if iter_sv == cur_sv && iter_x == cur_x + 1
                color = [0 1 0];
              end
              
              if iter_sv == cur_sv && iter_x == cur_x - 1
                color = [0 0 1];
              end
              if iter_x == cur_x && iter_sv == cur_sv - 1
                color = [0 0 1];
              end
              if iter_x == cur_x && iter_sv == cur_sv + 1
                color = [0 0 1];
              end
            case .75 % forward
              if iter_x == cur_x && iter_sv == cur_sv + 1
                color = [0 1 0];
              end
              
              if iter_sv == cur_sv && iter_x == cur_x - 1
                color = [0 0 1];
              end
              if iter_x == cur_x && iter_sv == cur_sv - 1
                color = [0 0 1];
              end
              if iter_sv == cur_sv && iter_x == cur_x + 1
                color = [0 0 1];
              end
          end
          if Nsv*(w_idx-1) + h_idx == step
            color = [.8 0 0];
          end
          if Nsv*(w_idx-1) + h_idx == cur_col
            % Current column
            color = [.8 0 0];
          end
       plotcube([.5 .5 Nx], [w_idx - .25, h_idx - .25, 0], .8, color);
%           plot(w_idx, h_idx, 'Marker', 's', 'MarkerFaceColor', color, 'Color', color);
      end
  end
  

  xlim([.5 Nx + .5]);
  %     xticks(1:Nx);
  xlabel('Along-Track');
%   set(gca, 'XColor', 'b');
  set(gca, 'xdir', 'reverse');

  ylim([.5 Nsv + .5]);
  %     yticks(1:Nsv);
  ylabel('Cross-Track');
%   set(gca, 'YColor', 'g');
%   set(gca, 'ydir', 'reverse');


   zlim([0 Nt]);
%   %     zticks(0:Nt);
   zlabel('Fast-Time');
%    set(gca, 'ZColor', 'r');
   set(gca, 'zdir', 'reverse');

%   view([-45 90 90]);
  view([90 90]);
  cameratoolbar('SetCoordSys', 'none');
  cameratoolbar('SetMode', 'nomode');
end


function [col_x, col_sv] = coords(ind)
  global Nsv;
  col_sv = ceil(mod(ind - 1, Nsv)) + 1;
  col_x = ceil(ind / Nsv);
end

function keyboard_input(~,event)
  global step Nsv;
  if ~isempty(event.Key)
      step_size = .25;
      if any(strcmpi(event.Modifier,'shift'))
         step_size = 1;
      end
      if any(strcmpi(event.Modifier,'control'))
         step_size = Nsv;
      end
      switch event.Key
      
        case 'rightarrow' % Right arrow
          step = step + step_size;
        case 'leftarrow'
          step = step - step_size;
      end
    draw();
  end
end



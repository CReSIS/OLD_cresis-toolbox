function [path,lower,upper]=run_ice(fname, top_smooth, bottom_smooth, ...
    top_peak, bottom_peak, repulse, pts1, pts2, alg, plot, png, savepath)

% function [path,lower,upper]=run_ice(fname, top_smooth, bottom_smooth, ...
%    top_peak, bottom_peak, repulse, pts1, pts2, alg, plot, png, savepath)
%
%  Call MCMC or HMM algorithm with desired parameters
%  Usually called from MCMC_HMM_automated_frames.m
%
%  See also: MCMC_HMM_automated_frames.m, run_MCMC_HMM_automated_frames.m
%  Actual layer tracking algorithms: RJ_MCMC.cpp, stereo.cpp
%
%  Author: Mingze Xu, Victor Berger

% LOAD .mat file for segment
Original = load(fname);
myData = Original.Data;
LogData = 10*log10(abs(myData));

% SCALE between 0-255:
Scaled = 255  * ((LogData - min(LogData(:))) ./ (max(LogData(:) - min(LogData(:)))));

% COPY into RGB matrix with the same values
A(:,:,1) = Scaled;
A(:,:,2) = Scaled;
A(:,:,3) = Scaled;

% CONVERT to uint8 format
A = uint8(A);

warning('off', 'Images:initSize:adjustingMag');

while 1
    
    if alg
        
        % MCMC
        if(numel(pts1) + numel(pts2) > 0 )
            [path, lower, upper]=RJ_MCMC(double(A(:,:,1)), pts1, pts2);
        else
            [path, lower, upper]=RJ_MCMC(double(A(:,:,1)));
        end
        
        if plot
            figure
        else
            figure('visible', 'off');
        end
        
        imshow(A, 'Border', 'tight', 'InitialMagnification', 100);
        title('MCMC')
        line(1:size(A,2), path(1,:), 'Color', [1 0 0], 'LineWidth', 3);
        line(1:size(A,2), path(2,:), 'Color', [0 1 0], 'LineWidth', 3);
        
        if(numel(pts1) > 0)
            line(pts1(1,:), pts1(2,:), 'LineStyle', 'none', 'Marker', '*', ...
                'Color', [1 0 0]);
        end
        if(numel(pts2) > 0)
            line(pts2(1,:), pts2(2,:), 'LineStyle', 'none', 'Marker', '*', ...
                'Color', [0 1 0]);
        end
        
        if (png && ~isempty(savepath))
            F = getframe(gcf);
            imwrite(F.cdata, savepath);
        end
        
        
    else
        
        % HMM
        opts = [top_smooth bottom_smooth top_peak bottom_peak repulse];
        if(numel(pts1) + numel(pts2) > 0 )
            [~, path]=stereo(1, double(A(:,:,1)), opts, pts1, pts2);
        else
            [~, path]=stereo(1, double(A(:,:,1)), opts);
        end
        
        if plot
            figure
        else
            figure('visible', 'off');
        end
        
        imshow(A, 'Border', 'tight', 'InitialMagnification', 100);
        title('HMM')
        
        line(1:size(A,2), path(1,:), 'Color', [1 0 0], 'LineWidth', 3);
        line(1:size(A,2), path(2,:), 'Color', [0 1 0], 'LineWidth', 3);
        
        upper = 1;
        lower = 1;
        
        if(numel(pts1) > 0)
            line(pts1(1,:), pts1(2,:), 'LineStyle', 'none', 'Marker', '*', ...
                'Color', [1 0 0]);
        end
        if(numel(pts2) > 0)
            line(pts2(1,:), pts2(2,:), 'LineStyle', 'none', 'Marker', '*', ...
                'Color', [0 1 0]);
        end
        
        if (png && ~isempty(savepath))
            F = getframe(gcf);
            imwrite(F.cdata, savepath);
        end
        
    end
    
    x = [];
    if(numel(x) == 0)
        break
    end
    
    if(b == 3)
        pts1 = [pts1 double([x y]')];
    else
        pts2 = [pts2 double([x y]')];
    end
end

warning('on', 'Images:initSize:adjustingMag');

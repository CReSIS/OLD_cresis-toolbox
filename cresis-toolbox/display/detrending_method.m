% clear; close all; clc; warning off;

function [detrended_data, record_y] = detrending_method(data, seg_len, fit_order)

if size(data, 1) < size(data, 2)
    data = data';
end

% seg_len = 1001; % odd num
nonoverlap_len = (seg_len - 1)/2;

% load brod2005.txt
% data = brod2005;
data_len = length(data);

% calculate the coefficient, given a window size and fitting order
[coeff_output, A] = detrending_coeff(seg_len, fit_order);
A_coeff = A * coeff_output;

for seg_index = 1
    % left trend
    % seg_index = 1;
    xi = 1 + (seg_index - 1) * (seg_len - 1) : seg_index * (seg_len - 1) + 1;
    xi_left = xi;
    seg_data = data(xi);
    %slope = polyfit(xi, seg_data', fit_order);
    %left_trend1 = polyval(slope, xi);    
    % to be replaced by pre-calculated coefficients
    left_trend = (A_coeff * seg_data)';
    
    % mid trend
    if seg_index * (seg_len - 1) + 1 + nonoverlap_len > data_len
        xi = 1 + (seg_index - 1) * (seg_len - 1) + nonoverlap_len : length(data);
        xi_mid = xi;
        seg_data = data(xi);
        %slope = polyfit(xi, seg_data', fit_order);
        %mid_trend = polyval(slope, xi);
        if (length(seg_data) < seg_len) 
            [coeff_output1, A1] = detrending_coeff(length(seg_data), fit_order);
            A_coeff1 = A1 * coeff_output1;
            mid_trend = (A_coeff1 * seg_data)';
        else
            mid_trend = (A_coeff * seg_data)';
        end
        
        xx1 = left_trend((seg_len + 1)/2:seg_len);
        xx2 = mid_trend(1:(seg_len + 1)/2);
        w = (0 : nonoverlap_len)/nonoverlap_len;
        xx_left = xx1 .* (1 - w) + xx2 .* w;
        
        record_x = xi_left(1:nonoverlap_len); 
        record_y = left_trend(1:nonoverlap_len);
        
        mid_start_index = find(xi_mid == (xi_left(end) + 1));
        if (length(mid_start_index) == 0)
            record_x = [record_x xi_left((end+3)/2 : end)];
            record_y = [record_y xx_left(2 : end)];
        else
            record_x = [record_x xi_left((end+1)/2 : end) xi_mid(mid_start_index : end)];
            record_y = [record_y xx_left(1 : end) mid_trend((seg_len + 3)/2:end)];
        end
        
        detrended_data = data - record_y';
        return;
    else
        xi = 1 + (seg_index - 1) * (seg_len - 1) + nonoverlap_len : seg_index * (seg_len - 1) + 1 + nonoverlap_len;
        xi_mid = xi;
        seg_data = data(xi);
        %slope = polyfit(xi, seg_data', fit_order);
        %mid_trend = polyval(slope, xi);
        mid_trend = (A_coeff * seg_data)';
        
        % right trend    
        if (seg_index + 1) * (seg_len - 1) + 1 > data_len  % not enough data points available in the right part
            xi = seg_index * (seg_len - 1) + 1 : data_len;
            xi_right = xi;
            seg_data = data(xi);
            %slope = polyfit(xi, seg_data', fit_order);
            %right_trend = polyval(slope, xi);
            if (length(seg_data) < seg_len)
                % need to re-calculate coefficient and A
                [coeff_output1, A1] = detrending_coeff(length(seg_data), fit_order);
                A_coeff1 = A1 * coeff_output1;
                right_trend = (A_coeff1 * seg_data)';
            else
                right_trend = (A_coeff * seg_data)';
            end
        
            xx1 = left_trend((seg_len + 1)/2:seg_len);
            xx2 = mid_trend(1:(seg_len + 1)/2);
            w = (0 : nonoverlap_len)/nonoverlap_len;
            xx_left = xx1 .* (1 - w) + xx2 .* w;  
        
            xx1 = mid_trend((seg_len + 1)/2:seg_len);
            xx2 = right_trend(1:(seg_len + 1)/2);
            w = (0 : nonoverlap_len)/nonoverlap_len;
            xx_right = xx1 .* (1 - w) + xx2 .* w;
        
            record_x = xi_left(1:nonoverlap_len); 
            record_y = left_trend(1:nonoverlap_len);
    
            record_x = [record_x xi_left((end+1)/2 : end) xi_mid((end+1)/2 + 1 : end)];
            record_y = [record_y xx_left(1 : end) xx_right(2:end)];
        
            right_start_index = find(xi_right == xi_mid(end) + 1);    
            record_x = [record_x xi_right(right_start_index:end)];
            record_y = [record_y right_trend(right_start_index:end)];
            detrended_data = data - record_y';
            return;
        else
            xi = seg_index * (seg_len - 1) + 1 : (seg_index + 1) * (seg_len - 1) + 1;
            xi_right = xi;
            seg_data = data(xi);
            %slope = polyfit(xi, seg_data', fit_order);
            %right_trend = polyval(slope, xi);
            right_trend = (A * coeff_output * seg_data)';
    
            xx1 = left_trend((seg_len + 1)/2:seg_len);
            xx2 = mid_trend(1:(seg_len + 1)/2);
            w = (0 : nonoverlap_len)/nonoverlap_len;
            xx_left = xx1 .* (1 - w) + xx2 .* w;
            % plot(xi_left(end/2 : end), xx_left, 'r', 'LineWidth', 2)

            xx1 = mid_trend((seg_len + 1)/2:seg_len);
            xx2 = right_trend(1:(seg_len + 1)/2);
            w = (0 : nonoverlap_len)/nonoverlap_len;
            xx_right = xx1 .* (1 - w) + xx2 .* w;
            % plot(xi_mid(end/2:end), xx_right, 'r', 'LineWidth', 2)
    
            record_x = xi_left(1:nonoverlap_len); 
            record_y = left_trend(1:nonoverlap_len);

            record_x = [record_x xi_left((end+1)/2 : end) xi_mid((end+1)/2 + 1 : end)];
            record_y = [record_y xx_left(1 : end) xx_right(2:end)];
        end
    end
end

for seg_index = 2 : floor((data_len-1)/(seg_len-1)) - 1
    % left trend
    % seg_index = 1;
    xi = 1 + (seg_index - 1) * (seg_len - 1) : seg_index * (seg_len - 1) + 1;
    xi_left = xi;
    seg_data = data(xi);
    %slope = polyfit(xi_left, seg_data', fit_order);
    %left_trend = polyval(slope, xi);
    left_trend = (A_coeff * seg_data)';
    
    % mid trend
    xi = 1 + (seg_index - 1) * (seg_len - 1) + nonoverlap_len : seg_index * (seg_len - 1) + 1 + nonoverlap_len;
    xi_mid = xi;
    seg_data = data(xi);
    %slope = polyfit(xi, seg_data', fit_order);
    %mid_trend = polyval(slope, xi);
    mid_trend = (A_coeff * seg_data)';
    
    % right trend
    xi = seg_index * (seg_len - 1) + 1 : (seg_index + 1) * (seg_len - 1) + 1;
    xi_right = xi;
    seg_data = data(xi);
    %slope = polyfit(xi, seg_data', fit_order);
    %right_trend = polyval(slope, xi);
    right_trend = (A_coeff * seg_data)';
    
    xx1 = left_trend((seg_len + 1)/2:seg_len);
    xx2 = mid_trend(1:(seg_len + 1)/2);
    w = (0 : nonoverlap_len)/nonoverlap_len;
    xx_left = xx1 .* (1 - w) + xx2 .* w;
%     plot(xi_left(end/2 : end), xx_left, 'r', 'LineWidth', 2)

    xx1 = mid_trend((seg_len + 1)/2:seg_len);
    xx2 = right_trend(1:(seg_len + 1)/2);
    w = (0 : nonoverlap_len)/nonoverlap_len;
    xx_right = xx1 .* (1 - w) + xx2 .* w;
%     plot(xi_mid(end/2:end), xx_right, 'r', 'LineWidth', 2)
    
    record_x = [record_x xi_left((end+3)/2 : end) xi_mid((end+1)/2 + 1 : end)];
    record_y = [record_y xx_left(2 : end) xx_right(2:end)];
end

% last part of data
for seg_index = floor((data_len-1)/(seg_len-1))
    % left trend
    % seg_index = 1;
    xi = 1 + (seg_index - 1) * (seg_len - 1) : seg_index * (seg_len - 1) + 1;
    xi_left = xi;
    seg_data = data(xi);
    %slope = polyfit(xi_left, seg_data', fit_order);
    %left_trend = polyval(slope, xi);
    left_trend = (A_coeff * seg_data)';
    
    % mid trend
    if (seg_index * (seg_len - 1) + 1 + nonoverlap_len > length(data)) 
        xi = 1 + (seg_index - 1) * (seg_len - 1) + nonoverlap_len : length(data);
        xi_mid = xi;
        seg_data = data(xi);
        %slope = polyfit(xi, seg_data', fit_order);
        %mid_trend = polyval(slope, xi);

        % may need to re-calcualte coefficient and A, since the length of
        % 'seg_data' may be smaller than seg_len
        if (length(seg_data) < seg_len)
            % need to re-calculate coefficient and A
            [coeff_output1, A1] = detrending_coeff(length(seg_data), fit_order);
            A_coeff1 = A1 * coeff_output1;
            mid_trend = (A_coeff1 * seg_data)';
        else
            mid_trend = (A_coeff * seg_data)';
        end
        
        xx1 = left_trend((seg_len + 1)/2:seg_len);
        xx2 = mid_trend(1:(seg_len + 1)/2);
        w = (0 : nonoverlap_len)/nonoverlap_len;
        xx_left = xx1 .* (1 - w) + xx2 .* w;
        
        mid_start_index = find(xi_mid == (xi_left(end) + 1));
        if (length(mid_start_index) == 0)
            record_x = [record_x xi_left((end+3)/2 : end)];
            record_y = [record_y xx_left(2 : end)];
        else
            record_x = [record_x xi_left((end+3)/2 : end) xi_mid(mid_start_index : end)];
            record_y = [record_y xx_left(2 : end) mid_trend((seg_len + 3)/2:end)];
        end
    
        break;
    else
        xi = 1 + (seg_index - 1) * (seg_len - 1) + nonoverlap_len : seg_index * (seg_len - 1) + 1 + nonoverlap_len;
        xi_mid = xi;
        seg_data = data(xi);
        %slope = polyfit(xi, seg_data', fit_order);
        %mid_trend = polyval(slope, xi);
        mid_trend = (A_coeff * seg_data)';
    end
%     plot(xi, mid_trend);
%     hold on;
    
    % right trend
    xi = seg_index * (seg_len - 1) + 1 : data_len;
    xi_right = xi;
    seg_data = data(xi);
    %slope = polyfit(xi, seg_data', fit_order);
    %right_trend = polyval(slope, xi);
    
    if (length(seg_data) < seg_len)
        % need to re-calculate coefficient and A
        [coeff_output1, A1] = detrending_coeff(length(seg_data), fit_order);
        A_coeff1 = A1 * coeff_output1;
        right_trend = (A_coeff1 * seg_data)';
    else
        right_trend = (A_coeff * seg_data)';
    end
    
    xx1 = left_trend((seg_len + 1)/2:seg_len);
    xx2 = mid_trend(1:(seg_len + 1)/2);
    w = (0 : nonoverlap_len)/nonoverlap_len;
    xx_left = xx1 .* (1 - w) + xx2 .* w;
%     plot(xi_left(end/2 : end), xx_left, 'r', 'LineWidth', 2)

    xx1 = mid_trend((seg_len + 1)/2:seg_len);
    xx2 = right_trend(1:(seg_len + 1)/2);
    w = (0 : nonoverlap_len)/nonoverlap_len;
    xx_right = xx1 .* (1 - w) + xx2 .* w;
%     plot(xi_mid(end/2:end), xx_right, 'r', 'LineWidth', 2)
    
    record_x = [record_x xi_left((end+3)/2 : end) xi_mid((end+1)/2 + 1 : end)];
    record_y = [record_y xx_left(2 : end) xx_right(2:end)];
    
    right_start_index = find(xi_right == xi_mid(end) + 1);    
    record_x = [record_x xi_right(right_start_index:end)];
    record_y = [record_y right_trend(right_start_index:end)];
end

detrended_data = data - record_y';
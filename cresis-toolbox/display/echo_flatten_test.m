GENERATE = false;
USE_ELEVATION = false;

physical_constants;
vel_air = c/2;

if ~GENERATE
    param = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140516_01');
    mdata = echo_load(param,'CSARP_post/standard',31);
    
    mdata.Data = real(10*log10(mdata.Data));
    
    max_val = max(max(mdata.Data));
    mean_val = mean(mean(mdata.Data));
    min_val = min(min(mdata.Data));
    
    if USE_ELEVATION
        slope_type = 'Elevation';
        elevation = mdata.Elevation;
        time = mdata.Time;

        dt = time(2) - time(1);
        drange = dt * vel_air;
        bins = elevation / drange;
    else
        slope_type = 'Surface';
        bins = interp1(mdata.Time, 1:length(mdata.Time), mdata.Surface);
    end
    

    % Insert test layer which mirrors elevation
    for i = 1:length(bins)
        mdata.Data(round(max(bins) - bins(i) + min(bins)), i) = max_val;
    end
    
    figure(1);
    clf;
    hold on;
    imagesc(mdata.Data);
    colormap(1-gray);
    plot(bins);
    legend({slope_type});
    title('Echo before shift');
    axis ij;
    
    flattened = echo_flatten(mdata, bins);
    
    figure(2);
    clf;
    hold on;
    imagesc(flattened);
    colormap(1-gray);
    plot(bins);
    legend({sprintf('Original %s', slope_type)});
    title('Echo after shift');
    axis ij;
    
else
    rows = 20;
    cols = 20;
    shift_matrix = .5;
    shift_elev = .5;
    flatness = 100;
    layer_start = 0;
    down_shift = rows;

    % Create default matrix
    matrix = magic(max([rows cols]));
    matrix = matrix(1:rows, 1:cols);
    matrix = matrix ./ max(max(matrix)) .* 5;

    % Create borders
    matrix(1, : ) = 10;
    matrix(end, : ) = 10;

    % Create elevation model by which to shift the layer
    layer_elev = ((down_shift - (rows:-1:1) + 1) * shift_matrix + 5) ./ 1;
    elev_range = layer_elev - min(layer_elev);
    elev_range = interp_finite(elev_range);

    % Create layer sine wave
    layer_mat = zeros(size(matrix, 1) + floor(max(elev_range)*2), size(matrix, 2));
    for i = 1:cols 
      layer_mat(max(floor(abs(sin(i)/flatness)*rows) + 1, 1) + round(size(layer_mat, 1)/2), i) = 30;
    end

    % Skew the layer sine wave based on shift_matrix
    skewed_layer_mat = layer_mat;
    for i = 1:cols
      skewed_layer_mat(:, i) = interp1(1:size(skewed_layer_mat, 1), skewed_layer_mat(:, i), (1:size(matrix, 1) + floor(max(elev_range)*2)) + elev_range(i)); 
    end

    % remove nans
    skewed_layer_mat(isnan(skewed_layer_mat)) = 0;

    % Crop skewed layer matrix
    [r, ~] = ind2sub([size(skewed_layer_mat, 1), size(skewed_layer_mat, 2)], find(skewed_layer_mat));
    skewed_layer_mat = skewed_layer_mat(min(r):max(r), :);

    % Apply layer to matrix
    for c = 1:cols
       for r = 1:size(skewed_layer_mat, 1)
           row = r + layer_start;
           if skewed_layer_mat(r, c) > 0 && row <= size(matrix, 1)
               matrix(row, c) = skewed_layer_mat(r, c);
           end
       end
    end

    % Create elevation model by which to shift the matrix back
    elevation = ((down_shift - (rows:-1:1) + 1) * shift_elev + 5) ./ 1;

    flattened = echo_flatten(matrix, elevation);

    % Crop layer matrix
    [r, ~] = ind2sub([size(layer_mat, 1), size(layer_mat, 2)], find(layer_mat));
    layer_mat = layer_mat(min(r)-5:max(r)+5, :);
    figure(1);
    clf;
    hold on;
    imagesc(layer_mat);
    colormap(1-gray);
    title('Layer before adding elevation shift');

    figure(2);
    clf;
    hold on;
    imagesc(matrix);
    colormap(1-gray);
    plot(elevation);
    legend({'Plane elevation'});
    title('Elevation-shifted layer in echogram');

    figure(3);
    clf;
    hold on;
    imagesc(flattened);
    colormap(1-gray);
    title('Echogram shifted to account for elevation');
end



function [ pc ] = pcsynthetic( block )
    % Michael Christoffersen
    % May 2018
    % Pulse compress with a synthetic chirp, generated from block file info

    data = block.ch0;

    %% Subtract along track mean
    m = zeros(length(data(:,1)),1);
    for i=1:length(data(1,:))
        m = m+data(:,i);
    end
    m=m/length(data(1,:));

    for i=1:length(data(1,:))
       data(:,i) = data(:,i) - m;
    end

    %% Filtering

    %for i=1:length(data(1,:))
    %    data(:,i) = filter(Hbp,data(:,i));
    %end

    %% Generating synthetic chirp
    % This assumes a block file is loaded
    T = [0:block.dt:block.chirp.len]; % Time window for synthetic chirp

    c = chirp(T,block.chirp.cf-(block.chirp.cf*.5 *block.chirp.bw/100),block.chirp.len,block.chirp.cf+(block.chirp.cf*.5*block.chirp.bw/100));
    c(length(c)+1:length(block.ch0(:,1))) = 0;
    %w = gausswin(251); % Windows
    %w = hanning(length(c));
    %c = c.*w';  % Multiply chirp by window

    pc = zeros(size(block.ch0));

    %% Pulse compression
    C = fft(c)';
    for i=1:length(block.ch0(1,:))
        S = fft(data(:,i));
        pc(:,i) = ifft(C.*S); 
    end

    pcabs = abs(log(pc));

    %% Smoothing the speckle
    k = 0.1*ones(3);
    pc = conv2(pcabs,k,'same');
end

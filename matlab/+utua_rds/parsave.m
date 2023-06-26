function [] = parsave(path,block)
    % Michael Christoffersen
    % August 2018
    % Wrapper for saving a block file, allows use of matlab parfor in
    % syncnav_loop and tdms2block_loop

    %Saves a block file to a the given file name
    save(path,'block','-v7.3');
end


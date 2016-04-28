function [ data_out ] = Sort_NACA_Data( data )
%%%
% Sorts data points within each entry of the NACA data structure, such that they progress
% around the circumference of the deformed wing. This is done by computing the camber line
% of the deformed airfoil, and identifying nodes above and below it. These two sets of
% nodes are then sorted by increasing and decreasing x-coordinate, respectively.
%%%
    
    %%%
    % Load the remaining '.csv' files and build the output data structure.
    %%%
    
    try
        load_str = 'Sorting nodal point data...';
        hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(files)), ...
                        'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
        setappdata(hwait,'cancel_loading',0);
        for i = 1:length(files)
            if getappdata(hwait,'cancel_loading')
                delete(hwait);
                error('Nodal point sorting canceled by user. Program stopped.');
            end
            waitbar(i/length(files),hwait,sprintf('%s %i / %i',load_str,i,length(files)));
            [data_out(i).IDstr, ...
             data_out(i).x, ...
             data_out(i).xnorm, ...
             data_out(i).y, ...
             data_out(i).p, ...
             data_out(i).cp, ...
             data_out(i).ref] = Load_NACA_Data(directory_path, files(i).name); %#ok<*AGROW>
            data_out(i).filename = files(i).name;
        end
    catch
        delete(hwait);
    end
    delete(hwait);

end

    
    
    
    
    
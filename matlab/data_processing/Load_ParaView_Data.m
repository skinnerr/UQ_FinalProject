function [ varargout ] = Load_ParaView_Data( data_path, requested_fields )
%%%
% Loads data saved from ParaView.
%
% Usage:
%   [T, p, EV] = Load_ParaView_Data('data.csv',{'T','p','EV'})
%%%

    % Assert that requested_fields is a nonempty string or a nonempty cell array.
    validateattributes(requested_fields,{'char','cell'},{'nonempty'});

    %%%
    % Open file, read first line, and close file
    %%%
    
    file_ID = fopen(data_path);
    if file_ID == -1
        error('File does not exist: %s',data_path);
    end
    try
        first_line = fgetl(file_ID);
        if ~ischar(first_line)
            error('File ended before expected, unable to read field names: %s',data_path);
        end
    catch exception
        fclose(file_ID);
        throw(exception);
    end
    fclose(file_ID);
        
    %%%
    % Parse field names and perform a number of checks.
    %%%

    % Parse the line into variable names.
    field_names = regexp(first_line,'[^",]+','match');

    requested_field_indices = zeros(length(requested_fields),1);
    
    % If we're loading 'all' fields...
    if strcmp(requested_fields,'all')
        requested_field_indices = 1:length(field_names);

    % If we're loading from a list of desired fields...
    else
        % Loop through the requested field names...
        for i = 1:length(requested_fields)
            field = char(requested_fields(i));
            % Check for the requested field in the available variable names.
            matches = strcmp(field_names,field);
            if ~any(matches)
                error('Requested field (''%s'') was not found in the file.', field);
            elseif sum(matches) > 1
                error('Field (''%s'') should not be requested more than once.', field);
            else
                match_index = find(matches);
                requested_field_indices(i) = match_index;
            end
        end
    end

    %%%
    % Read subsequent lines and build varargout output.
    %%%
    
    data = csvread(data_path,1,0);

    varargout = cell(length(requested_field_indices),1);
    for i = 1:length(varargout)
        varargout{i} = data(:,requested_field_indices(i));
    end

end
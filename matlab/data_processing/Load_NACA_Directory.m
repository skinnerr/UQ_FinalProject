function [ NACA_data ] = Load_NACA_Directory( varargin )
%%%
% Loads all the .csv files from a directory.
%
% Input: Load_NACA_Data(directory_path)
%        Load_NACA_Data(directory_path, filter_string)
%         (where this function ignores all filenames not containing filter_string)
%
% Output: NACA_Data is an array containing a struct for each file loaded in.
%%%

    %%%
    % Sanitize inputs.
    %%%

    if nargin < 0 || nargin > 3
        error('Invalid number of arguments passed in.');
    end
    
    directory_path = varargin{1};
    if directory_path(end) ~= '/'
        directory_path(end+1) = '/';
    end
    
    if nargin < 2
        file_filter_string = '';
    else
        file_filter_string = varargin{2};
    end
    
    if nargin < 3
        n_files = 0;
    else
        n_files = varargin{3};
    end

    % Assert that inputs are strings.
    validateattributes(directory_path,    {'char'},{'nonempty'});
    validateattributes(file_filter_string,{'char'},{});
    validateattributes(n_files,           {'numeric'},{'>=',0});
    
    %%%
    % Obtain a list of all '.csv' files in the directory.
    %%%
    
    % List all files.
    files = dir(directory_path);
    % Retain only files that are not directories.
    files = files(find( ~[files.isdir] )); %#ok<*FNDSB>
    % Retain only files that end in '.csv'.
    files = files(find( ~cellfun(@isempty, regexp({files.name},'.csv$')) ));
    % Retain only files that contain the filter string, if it is non-empty.
    if ~isempty(file_filter_string)
        files = files(find( ...
                    ~cellfun(@isempty, regexp({files.name},file_filter_string)) ));
    end
    
    % Truncate file list if needed.
    if n_files > 0 && n_files < length(files)
        files = files(1:n_files);
    end
    
    %%%
    % Load the remaining '.csv' files and build the output data structure.
    %%%
    
    try
        load_str = 'Loading data...';
        hwait = waitbar(0,sprintf('%s %i / %i',load_str,1,length(files)), ...
                        'CreateCancelBtn','setappdata(gcbf,''cancel_loading'',1)');
        setappdata(hwait,'cancel_loading',0);
        for i = 1:length(files)
            if getappdata(hwait,'cancel_loading')
                delete(hwait);
                error('Directory loading canceled by user. Program stopped.');
            end
            waitbar(i/length(files),hwait,sprintf('%s %i / %i',load_str,i,length(files)));
            [NACA_data(i).IDstr, ...
             NACA_data(i).x, ...
             NACA_data(i).xnorm, ...
             NACA_data(i).y, ...
             NACA_data(i).p, ...
             NACA_data(i).cp, ...
             NACA_data(i).ref] = Load_NACA_Data(directory_path, files(i).name); %#ok<*AGROW>
            NACA_data(i).filename = files(i).name;
        end
    catch
        delete(hwait);
    end
    delete(hwait);

end

    
    
    
    
    
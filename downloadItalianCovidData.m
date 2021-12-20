function fname = downloadItalianCovidData(varargin)
    
    %
    % downloadItalianCovidData()
    %
    %
    %
    % downloadItalianCovidData()
    % -- Download the data of today
    %
    % downloadItalianCovidData(day)
    % -- Download the data of the specified day. day should be a string in
    % -- the form 'yyyymmdd'
    % -- day can also be set to 'today'
    %
    % downloadItalianCovidData(day, dest_folder)
    % -- Store the downloaded data in dest_folder.
    %
    % downloadItalianCovidData(day, dest_folder, 'force')
    % -- Force the download of the data, even if the corresponding file is 
    % -- already on the disk
    %
    %
    % Returns the name of the filename

    
    %% Process input arguments
    if (nargin == 0)
        date_string = datestr(now, 'yyyymmdd');
        dest_folder = '.';
        forceDownload = false;
    elseif (nargin == 1)
        date_string = varargin(1);
        dest_folder = '.';
        forceDownload = false;     
    elseif (nargin == 2)
        date_string = varargin(1);
        dest_folder = varargin(2);
        forceDownload = false; 
    elseif (nargin == 3)
        date_string = varargin(1);
        dest_folder = varargin(2);
        force_string = varargin(3);
        
        if (strcmp(force_string, 'force'))
            forceDownload = true;
        else
            s = sprintf('If specified, second argument should be ''force''');
            error(s)
        end
    end

    if (iscell(date_string))
        date_string = date_string{1};
    end
    
    if (iscell(dest_folder))
        dest_folder = dest_folder{1};
    end    
    
    if (strcmp(date_string, 'today'))
        date_string = datestr(now, 'yyyymmdd');
    else
        if (numel(date_string) ~= 8 || ...
                (min(double(date_string)) < 48 || max(double(date_string)) > 57))
            s = sprintf('First argument should be a string in the form ''yyyymmdd'' or ''today''');
            error(s)            
        end
    end    
    
    if (dest_folder(end) ~= '/')
        dest_folder = [dest_folder, '/'];
    end
    
    %% Download

    url = sprintf('https://raw.githubusercontent.com/pcm-dpc/COVID-19/master/dati-regioni/dpc-covid19-ita-regioni-%s.csv', date_string);    
    fname = sprintf('dpc-covid19-ita-regioni-%s.csv', date_string);
    fname = [dest_folder fname];
    
    try
        file_exists = isfile(fname);
    catch ex
        file_exists = (exist(fname, 'file') == 2); % matlab < R2017
    end

    try
        folder_exists = isfolder(dest_folder);
    catch ex
        folder_exists = false; %TODO matlab < R2017
    end    
    
    if (~folder_exists)
        s = sprintf('> Destionation folder %s does not exist\n', dest_folder);
        error(s)        
    end
    
    fprintf('> Downloading file %s...\n', url);
    if (forceDownload || ~file_exists)    
        try
            outfilename = websave(fname, url);
        catch ex
            s = sprintf('>> Problem downloading the file from url %s\n', url);
            error(s)
        end
    else
        fprintf('>> The file %s is already in the folder %s!\n', fname, dest_folder);
    end

end
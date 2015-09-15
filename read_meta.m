function varargout          = read_meta(name_path, name_file, file_type)
% READ_META
% 
% Read in metadata text file.
% 
% David Stillman, Joe MacGregor
% Last updated: 02/06/15

fid                         = fopen([name_path name_file '.txt'], 'r'); % open file

meta_test                   = fscanf(fid, '%s', 1); % read in string to see if it's new or not
if isnan(str2double(meta_test(1)))
    meta_type               = 'new';
    if mod((length(strfind(meta_type, ',')) - 6), 4) % if there are any commas left, then this must be the new type (which includes sample diameter)
        meta_type           = 'new+';
    end
else
    meta_type               = 'old';
end

fseek(fid, 0, -1); % return to start of file

% load metadata
switch file_type
    case 'data'
        switch meta_type
            case 'new+' % also includes sample diameter
                meta_all    = textscan(fid, '%s%s%f%f%f%d%d%d', 1, 'delimiter', ',');
                [name, name_os, thick, dia_sample, dia_elec, num_files, num_files_tot, num_loops] ...
                            = deal(meta_all{:});
            case 'new'
                meta_all    = textscan(fid, '%s%s%f%f%d%d%d', 1, 'delimiter', ',');
                [name, name_os, thick, dia_elec, num_files, num_files_tot, num_loops] ...
                            = deal(meta_all{:});
                dia_sample  = 55; % default sample diameter assumed (in mm)
            case 'old'
                meta_all    = textscan(fid, '%d%s%d%s%f%f%d%d%d', 1, 'delimiter', ',');
                [~, name, ~, name_os, thick, dia_elec, num_files, num_files_tot, num_loops] ...
                            = deal(meta_all{:});
                dia_sample  = 55; % default sample diameter assumed (in mm)                        
        end;
        name                = name{1};
        name_os             = name_os{1};
        thick               = thick * 0.001; % convert sample thickness to m
        dia_sample          = dia_sample * 0.001; % convert sample diameter to m
        dia_elec            = dia_elec * 0.001; % convert electrode diameter to m
        area_elec           = pi * ((dia_elec / 2) ^ 2); % lower electrode area in m^2
        num_files_tot       = (num_files_tot * 2) - 1;
    case 'os'
        switch meta_type
            case 'new'
                meta_all    = textscan(fid, '%s%f%f%d', 1, 'delimiter', ',');
                [~, ~, ~, num_loops] ...
                            = deal(meta_all{:});
            case 'old'
                meta_all    = textscan(fid, '%d%s%f%f%d', 1, 'delimiter', ',');
                [~, ~, ~, ~, num_loops] ...
                            = deal(meta_all{:});
        end
end

% scan loop parameters
loop_scan                   = textscan(fid, '%f', (num_loops * 4), 'delimiter', ','); 
loop_scan                   = loop_scan{1};
fclose(fid);

% load "loop" parameters for different sampling densities across different bands
[freq_max, freq_min, num_pts_decade, stack, num_pts_stack] ...
                            = deal(zeros(1, num_loops));
for ii = 1:num_loops
    freq_max(ii)            = loop_scan(1 + (4 * (ii - 1))); % max frequency for this loop
    freq_min(ii)            = loop_scan(2 + (4 * (ii - 1))); % min frequency for this loop
    num_pts_decade(ii)      = loop_scan(3 + (4 * (ii - 1))); % number of pts per frequency decade in this loop
    stack(ii)               = loop_scan(4 + (4 * (ii - 1))); % number of stack (repeats)
    num_pts_stack(ii)       = ceil((log10(freq_max(ii)) - log10(freq_min(ii))) * num_pts_decade(ii)) + 1; % number of points in stack
end
num_pts                     = sum(num_pts_stack);

% output depends on meta data type
switch file_type
    case 'data'
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, varargout{5}, varargout{6}, varargout{7}, varargout{8}, varargout{9}, varargout{10}, varargout{11}, varargout{12}, varargout{13}, varargout{14}, varargout{15}] ...
                            = deal(name, name_os, thick, dia_sample, dia_elec, area_elec, num_files, num_files_tot, num_loops, freq_max, freq_min, num_pts_decade, stack, num_pts_stack, num_pts);        
    case 'os'
        [varargout{1}, varargout{2}] ...
                            = deal(stack, num_pts);
end
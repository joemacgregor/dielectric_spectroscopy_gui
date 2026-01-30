function [time, freq, Z_abs, Z_phase, Z_real, Z_imag, Z_comp] ...
                            = read_solartron_csv(name_path, name_file, stack, num_pt, num_file)
% READ_SOLARTRON_CSV Read in Solartron .csv data file.
% 
% [TIME,FREQ,Z_ABS,Z_PHASE,Z_REAL,Z_IMAG,Z_COMP] = READ_SOLARTRON_CSV(NAME_PATH,NAME_FILE,STACK,NUM_PT,NUM_FILE)
% reads in a Solartron .csv data file. NAME_PATH and NAME_FILE are the
% path and filename, respectively, STACK is the number of stack, NUM_PT is
% the number of spectral samples and NUM_FILE is the number of files that
% compose the spectrum. 
% 
% The outputs are:
% 
% TIME: Measurement time
% FREQ: Frequency (Hz)
% Z_ABS/PHASE/REAL/IMAG/COMP: Absolute, phase, real part, imaginary part
% and complex impedance, respectively.
% 
% David Stillman, Joe MacGregor
% Last updated: 09/21/15

if (nargin ~= 5)
    error('read_solartron_csv:nargin', 'Number of arguments must be 5.')
end
if ~ischar(name_path)
    error('read_solartron_csv:name_path', 'NAME_PATH is not a string.')
elseif (~isempty(name_path) && ~exist(name_path, 'dir'))
    error('read_solartron_csv:name_path_dir', [name_path ' is not an existing directory.'])
end
if ~ischar(name_file)
    error('read_solartron_csv:name_file', 'NAME_FILE is not a string.')
elseif ~exist([name_path name_file '.csv'], 'file')
    error('read_solartron_csv:name_path_dir', [name_path name_file '.csv is not an existing file.'])
end
if (~isnumeric(stack) || ~isscalar(stack))
    error('read_solartron_csv:stack', 'STACK is not a numeric scalar.')
elseif ((stack < 1) || mod(stack, 1))
    error('read_solartron_csv:stack_int', 'STACK is not an integer greater than zero.')
end
if (~isnumeric(num_pt) || ~isscalar(num_pt))
    error('read_solartron_csv:num_pt', 'NUM_PT is not a numeric scalar.')
elseif ((num_pt < 1) || mod(num_pt, 1))
    error('read_solartron_csv:num_pt_int', 'NUM_PT is not an integer greater than zero.')
end
if (~isnumeric(num_file) || ~isscalar(num_file))
    error('read_solartron_csv:num_file', 'NUM_FILE is not a numeric scalar.')
elseif ((num_file < 1) || mod(num_file, 1))
    error('read_solartron_csv:num_file_int', 'NUM_FILE is not an integer greater than zero.')
end
if (nargout ~= 7)
    error('read_solartron_csv:nargout', 'Number of outputs must be 7.')
end

fid                         = fopen([name_path name_file '.csv'], 'r'); % load .csv data

% get start date/time into a useful format
fgetl(fid);                                                                             % skip line
start_str                   = fgetl(fid);                                               % keep second line
[~, start_str]              = strtok(start_str, ':');                                   % skip junk
start_str                   = start_str(3:end);                                         % remove junk
start_all                   = textscan(start_str, '%d/%d/%d %d:%d:%d %s');              % scan remaining string
start_vec                   = zeros(1, 6);                                              % initialize
[start_vec(2), start_vec(3), start_vec(1), start_vec(4), start_vec(5), start_vec(6), ~] ...
                            = deal(start_all{:});                                       % assign values to a start date vector (MATLAB datevec format)
if (strcmp(start_all{7}, 'PM') && (start_vec(4) < 12))                                  % all this fiddling about is why grown-up programs use 24-h time
    start_vec(4)            = start_vec(4) + 12;                                        % adjust hours for PM
elseif (strcmp(start_all{7}, 'AM') && (start_vec(4) == 12)) 
    start_vec(4)            = 0;                                                        % adjust hours because 00:00 is 12:00 AM in Solartron crazy-land date formatting 
end
start_num                   = datenum(start_vec);

% skip 2 more lines 
fgetl(fid); % skip line
fgetl(fid); % skip line

% initialize, read in and assign data
time_str                    = cell(num_file, num_pt, stack);
[time, freq, Z_abs, Z_phase, Z_real, Z_imag, Z_comp] ...
                            = deal(cell(num_file, 1));

% loop through stack number
for ii = 1:num_file
    [freq{ii}, Z_abs{ii}, Z_phase{ii}, Z_real{ii}, Z_imag{ii}, Z_comp{ii}] ...
                            = deal(zeros(num_pt, stack));
    for jj = 1:stack
        data_all            = textscan(fid, '%f%f%f%s%f%f%f%s%f%f%f%f', num_pt, 'delimiter', ','); % limited to num_pts, could be more records
        [~, ~, ~, time_str(ii, :, jj), freq{ii}(:, jj), ~, ~, ~, Z_abs{ii}(:, jj), Z_phase{ii}(:, jj), Z_real{ii}(:, jj), Z_imag{ii}(:, jj)] ...
                            = deal(data_all{:});                                        % assign variables, skip ones we don't care about (~)
        Z_phase{ii}(:, jj)  = deg2rad(Z_phase{ii}(:, jj));                              % degrees to radians
        Z_comp{ii}(:, jj)   = Z_real{ii}(:, jj) + (1i .* Z_imag{ii}(:, jj));            % complex impedance
    end
end
fclose(fid);

% convert time string into date numbers
for ii = 1:num_file
    time{ii}                = start_num .* ones(num_pt, stack); % absolute reference (date number)
    for jj = 1:stack
        for kk = 1:num_pt
            time_str_sep    = textscan(time_str{ii, kk, jj}, '%f%f%f', 'delimiter', ':'); % scan values from each time string
            [time_hr, time_min, time_sec] ...
                            = deal(time_str_sep{:}); % assign values from time string
            time{ii}(kk, jj)= time{ii}(kk, jj) + (time_hr / 24) + (time_min / (24 * 60)) + (time_sec / (24 * 60 * 60)); % add measurement time to start time to get absolute time
        end
    end
end
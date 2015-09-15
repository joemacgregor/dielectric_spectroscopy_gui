function [time, freq, Z_abs, Z_phase, Z_real, Z_imag, Z_comp] ...
                            = read_solartron_csv(name_path, name_file, stack, num_pts, num_files)
% READ_SOLARTRON_CSV
% 
% Read in Solartron .csv data file.
% 
% David Stillman, Joe MacGregor
% Last updated: 02/07/15

fid                         = fopen([name_path name_file '.csv'], 'r'); % load .csv data

% get start date/time into a useful format
fgetl(fid);                                                                             % skip line
start_str                   = fgetl(fid);                                               % keep second line
[~, start_str]              = strtok(start_str, ':');                                   % skip junk
start_str                   = start_str(3:end);                                         % remove junk
start_all                   = textscan(start_str, '%d/%d/%d %d:%d:%d %s');              % scan remaining string
start_vec                   = zeros(1, 6);                                              % initialize
[start_vec(2), start_vec(3), start_vec(1), start_vec(4), start_vec(5), start_vec(6), ~] ...
                            = deal(start_all{:});                                       % assign values to a start date vector (Matlab style)
if (strcmp(start_all{7}, 'PM') && (start_vec(4) < 12))                                  % all this fiddling about is why grown-up programs use 24-h time
    start_vec(4)            = start_vec(4) + 12;                                        % adjust hours for PM
elseif (strcmp(start_all{7}, 'AM') && (start_vec(4) == 12)) 
    start_vec(4)            = 0;                                                        % adjust hours because 00:00 is 12:00 AM in Solartron crazy-land date formatting 
end
start_num                   = datenum(start_vec);

% skip some more lines 
fgetl(fid); % skip line
fgetl(fid); % skip line

% read in and assign data
time_str                    = cell(num_files, num_pts, stack);
[time, freq, Z_abs, Z_phase, Z_real, Z_imag, Z_comp] ...
                            = deal(cell(num_files, 1));

% loop through stack number
for ii = 1:num_files
    [freq{ii}, Z_abs{ii}, Z_phase{ii}, Z_real{ii}, Z_imag{ii}, Z_comp{ii}] ...
                            = deal(zeros(num_pts, stack));
    for jj = 1:stack
        data_all            = textscan(fid, '%f%f%f%s%f%f%f%s%f%f%f%f', num_pts, 'delimiter', ','); % limited to num_pts, could be more records
        [~, ~, ~, time_str(ii, :, jj), freq{ii}(:, jj), ~, ~, ~, Z_abs{ii}(:, jj), Z_phase{ii}(:, jj), Z_real{ii}(:, jj), Z_imag{ii}(:, jj)] ...
                            = deal(data_all{:});                                        % assign variables, skip ones we don't care about (~)
        Z_phase{ii}(:, jj)  = Z_phase{ii}(:, jj) * (pi / 180);                          % degrees to radians
        Z_comp{ii}(:, jj)   = Z_real{ii}(:, jj) + (1i .* Z_imag{ii}(:, jj));            % complex impedance
    end
end
fclose(fid);

% convert time string into date numbers
for ii = 1:num_files
    time{ii}                = start_num .* ones(num_pts, stack); % absolute reference (date number)
    for jj = 1:stack
        for kk = 1:num_pts
            time_str_sep    = textscan(time_str{ii, kk, jj}, '%f%f%f', 'delimiter', ':'); % scan values from each time string
            [time_hr, time_min, time_sec] ...
                            = deal(time_str_sep{:}); % assign values from time string
            time{ii}(kk, jj)= time{ii}(kk, jj) + (time_hr / 24) + (time_min / (24 * 60)) + (time_sec / (24 * 60 * 60)); % add measurement time to start time to get absolute time
        end
    end
end
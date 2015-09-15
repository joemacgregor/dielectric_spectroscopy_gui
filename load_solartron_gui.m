function load_solartron_gui
% LOAD_SOLARTRON_GUI
% 
% Load Solartron measurement data and synchronize with temperature data.
% 
% David Stillman, Joe MacGregor
% Last updated: 09/15/15

%% initialize a bunch of stuff

path_os                     = 'cal/open_short/'; % open/short (os) path

% gap_guard                   = 0.001;    % gap between guard and guarded electrode (used in old calibration approach)
% rad_elec                    = 0.012;    % guarded electrode radius
threshold_std               = 0.04;     % threshold in standard deviation below which temperature match is satisfactory
pts_stable_min              = 10;       % minimum numbers pts (+1) over which standard deviation test is done

[name_file_data, path_data, name_file_temp, path_temp, name, name_os, path_fix, time_ignore_str, name_file_save, path_save] ...
                            = deal('');
[thick, area_elec, num_files, num_files_tot, num_loops, freq_max, freq_min, num_pts_decade, stack, num_pts_stack, num_pts, time, freq, Z_comp, temp_all, temp_sample, setpt, time_start_vec, ...
 time_start_num, time_temp, ind_temp_change, d_t, area_corr, stack_o, num_pts_o, freq_os, Z_comp_o, stack_s, num_pts_s, Z_real_s, Z_imag_s, Z_comp_s, freq_rad_os, C_os, resist_os, Z_corr, C_abs, ...
 C_phase, C_real, C_imag, C_loss_tan, permitt_real, permitt_imag, loss_tan, resist_abs, resist_phase, resist_real, resist_imag, freq_rad, setpt_fix, num_setpt, ind_setpt_start, ind_setpt_ord, ind_stable_start, ...
 ind_stable_end, ind_setpt_curr, std_temp, ind_setpt_stable, temp_mean, temp_std, temp_min, temp_max, temp_match, sweep_stable, temp_curr, ind_match, temp_match_unique, temp_match_best, colors, ...
 xticks, good_temp, tmp_load, pt, pis, fl, ph, lg, dia_sample, dia_elec, conduct, permitt_norm, dia_elec_ind, dia_sample_ind] ...
                            = deal(NaN);
time_ignore                 = [];

permitt_vacuum              = 8.854187818e-12; % permittivity of vacuum (F/m)

% load calibration curves
load cal/permitt_cal permitt_norm_cal dia_lo_cal dia_sample_cal thick_cal type_ratios permitt_norm_mean

if ispc % adjust for Mac/PC in path merging later on
    path_fix                = '\';
else
    path_fix                = '/';
end

subplot_start               = [0.05 0.27 0.43 0.19]; % start position of 1st subplot
plots                       = {'permitt_real{good_temp(jj)}' 'conduct{good_temp(jj)}' 'permitt_imag{good_temp(jj)}' '(resist_phase{good_temp(jj)} .* (180 / pi))'};
ylabels                     = {'\epsilon''' '\sigma (S/m)' '\epsilon''''' '\phi (deg)'}; % ylabels for each supblot
ax                          = NaN(1, 4); % axes handles for data subplots
pd                          = cell(1, 4); % data plots
freq_default                = [1e-3 1e7]; % freq range to display
ylim_default                = [1 1e4; 1e-12 1e-2; 1e-4 1e4; -90 0]; % ylims for each subplot

%% draw the GUI
loadgui                     = figure('toolbar', 'figure', 'name', 'LOAD_SOLARTRON_GUI', 'position', [1 1 1600 1200]);
tempax                      = subplot('position', [0.05 0.51 0.93 0.32]);
hold on;
set(gca, 'fontsize', 18, 'color', [0.9 0.9 0.9])
xlabel('Local time')
ylabel('Sample temperature (^oC)')
grid on
box on
for ii = 1:4 %#ok<*FXUP>
    subplot_curr            = subplot_start;
    if any(ii == [2 4])
        subplot_curr        = subplot_curr + [0.50 0 0 0];
    end
    if (ii > 2)
        subplot_curr        = subplot_curr + [0 -0.2 0 0];
    end
    ax(ii)                  = subplot('position', subplot_curr);
    hold on
    set(gca, 'fontsize', 16, 'xscale', 'log', 'color', [0.9 0.9 0.9])
    if (ii < 4) % all plots except phase should be y log
        set(gca, 'yscale', 'log')
    end
    axis([freq_default ylim_default(ii, :)])
    if (ii > 2)
        xlabel('Frequency (Hz)')
    else
        set(gca, 'xticklabel', {})
    end
    ylabel(ylabels{ii})
    grid on
    box on
end
linkaxes(ax, 'x') % links frequency axes together when zooming
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Load Solartron data', 'units', 'normalized', 'position', [0.02 0.95 0.16 0.04], 'callback', @load_solar, 'fontsize', 16)
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Load temperature data', 'units', 'normalized', 'position', [0.19 0.95 0.2 0.04], 'callback', @load_temp, 'fontsize', 16)
solar_box                   = annotation('textbox', [0.02 0.90 0.16 0.04], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
temp_box                    = annotation('textbox', [0.19 0.90 0.2 0.04], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
% annotation('textbox', [0.4 0.95 0.15 0.05], 'string', 'Gap guard (mm)', 'fontsize', 16, 'edgecolor', 'none');
% gap_guard_edit              = uicontrol(loadgui, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.92 0.05 0.04], 'string', sprintf('%3.1f', (1e3 * gap_guard)), 'fontsize', 16, 'foregroundcolor', 'b');
% annotation('textbox', [0.4 0.87 0.18 0.05], 'string', 'Electrode radius (mm)', 'fontsize', 16, 'edgecolor', 'none');
% rad_elec_edit              = uicontrol(loadgui, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.84 0.05 0.04], 'string', sprintf('%3.1f', (1e3 * rad_elec)), 'fontsize', 16, 'foregroundcolor', 'b');
annotation('textbox', [0.4 0.95 0.18 0.05], 'string', 'Temp. threshold (K)', 'fontsize', 16, 'edgecolor', 'none')
threshold_std_edit         = uicontrol(loadgui, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.92 0.05 0.04], 'string', sprintf('%6.4f', threshold_std), 'fontsize', 16, 'foregroundcolor', 'b');
annotation('textbox', [0.4 0.87 0.15 0.05], 'string', 'Min. pts stable', 'fontsize', 16, 'edgecolor', 'none')
pts_stable_min_edit        = uicontrol(loadgui, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.84 0.05 0.04], 'string', sprintf('%d', pts_stable_min), 'fontsize', 16, 'foregroundcolor', 'b');
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Merge', 'units', 'normalized', 'position', [0.02 0.85 0.06 0.04], 'callback', @do_merge, 'fontsize', 16, 'foregroundcolor', 'b')
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.09 0.85 0.06 0.04], 'callback', @do_save, 'fontsize', 16, 'foregroundcolor', 'g')
annotation('textbox', [0.15 0.85 0.06 0.04], 'string', 'Reset', 'fontsize', 16, 'edgecolor', 'none', 'color', 'r');
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Data', 'units', 'normalized', 'position', [0.21 0.85 0.06 0.04], 'callback', @do_reset_data, 'fontsize', 16, 'foregroundcolor', 'r')
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Temp.', 'units', 'normalized', 'position', [0.27 0.85 0.06 0.04], 'callback', @do_reset_temp, 'fontsize', 16, 'foregroundcolor', 'r')
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'All', 'units', 'normalized', 'position', [0.33 0.85 0.06 0.04], 'callback', @do_reset_all, 'fontsize', 16, 'foregroundcolor', 'r')
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Pick time to ignore', 'units', 'normalized', 'position', [0.56 0.95 0.15 0.04], 'callback', @do_ignore, 'fontsize', 16, 'foregroundcolor', 'b')
ignore_box                  = annotation('textbox', [0.56 0.90 0.15 0.04], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');
uicontrol(loadgui, 'style', 'pushbutton', 'string', 'Clear ignore time', 'units', 'normalized', 'position', [0.56 0.85 0.15 0.04], 'callback', @clear_ignore, 'fontsize', 16, 'foregroundcolor', 'r')
annotation('textbox', [0.72 0.95 0.2 0.04], 'string', 'Status', 'fontsize', 16, 'edgecolor', 'none')
status_box                  = annotation('textbox', [0.72 0.90 0.27 0.04], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none');

%% sub-functions

    function load_solar(source, eventdata) %#ok<*INUSD>

        do_reset_data % always reset data before loading a single file
        
        if isempty(path_data) % go through sequence if data path unavailable
            if ~isempty(path_temp)
                try
                    [name_file_data, path_data] ...
                            = uigetfile('*.txt', 'Load Solartron file:', path_temp); % dialog box
                catch %#ok<*CTCH>
                    set(status_box, 'string', 'Load data canceled')
                    [name_file_data, path_data] ...
                            = deal('');
                end
            else
                try
                    [name_file_data, path_data] ...
                            = uigetfile('*.txt', 'Load Solartron file:');
                catch
                    set(status_box, 'string', 'Load data canceled')
                    [name_file_data, path_data] ...
                            = deal('');
                end
            end
        else
            try
                [name_file_data, path_data] ...
                            = uigetfile('*.txt', 'Load Solartron file:', path_data); % dialog box if a file was previously loaded
            catch
                set(status_box, 'string', 'Load data canceled')
                    [name_file_data, path_data] ...
                            = deal('');
            end
        end

        if isempty(name_file_data) % must have hit cancel!
            [name_file_data, path_data] ...
                            = deal('');
            set(status_box, 'string', 'Load data canceled')
        end
        
        if ~isempty(name_file_data) % only do if there is a file to load (didn't cancel)
            
            name_file_data  = name_file_data(1:(end - 4)); % extract filename (without ".mat")
            set(solar_box, 'string', name_file_data) % show filename
            
            % load measurement metadata
            [name, name_os, thick, dia_sample, dia_elec, area_elec, num_files, num_files_tot, num_loops, freq_max, freq_min, num_pts_decade, stack, num_pts_stack, num_pts] ...
                            = read_meta(path_data, name_file_data, 'data');
            % load measurement data
            [time, freq, ~, ~, ~, ~, Z_comp] ...
                            = read_solartron_csv(path_data, name, stack(1), num_pts, num_files);
            
            % correct permittivities using calibration curves derived from Delrin and Teflon (see cal_work.m)
            dia_elec_ind    = find(dia_elec == dia_lo_cal); % narrow based on lower electrode diameter
            dia_sample_ind  = interp1(dia_sample_cal, 1:length(dia_sample_cal), dia_sample, 'nearest', 'extrap'); % narrow based on sample diameter
            if (~isempty(permitt_norm_cal{1}{dia_elec_ind, dia_sample_ind}) && ~isempty(permitt_norm_cal{2}{dia_elec_ind, dia_sample_ind}))
                permitt_norm = (type_ratios(1) * interp1(thick_cal{1}{dia_elec_ind, dia_sample_ind}, permitt_norm_cal{1}{dia_elec_ind, dia_sample_ind}, thick, 'linear', 'extrap')) + ...
                               (type_ratios(2) * interp1(thick_cal{2}{dia_elec_ind, dia_sample_ind}, permitt_norm_cal{2}{dia_elec_ind, dia_sample_ind}, thick, 'linear', 'extrap'));
            else
                permitt_norm = permitt_norm_mean; % settle for average ratio between the two if we don't have calibration values for both Delrin and Teflon
            end
            
            set(status_box, 'string', [num2str(num_files) ' data files loaded'])
        end
    end

    function load_temp(source, eventdata)
        
        do_reset_temp % reset temperature file first
        
        try % again try to load
            if ~isempty(path_data)
                [name_file_temp, path_temp] ...
                            = uigetfile({'*.mat'; '*.csv'}, 'Load temperature file:', path_data); % dialog box
            else
                [name_file_temp, path_temp] ...
                            = uigetfile({'*.mat'; '*.csv'}, 'Load temperature file:'); % dialog box  
            end
        catch
            set(status_box, 'string', 'Load temperature canceled')
            [name_file_temp, path_temp]  ...
                            = deal('');
        end
        
        if isempty(name_file_temp) % canceled!
            set(status_box, 'string', 'Load temperature canceled')
            [name_file_temp, path_temp] ...
                            = deal('');
        end
        
        if ~isempty(name_file_temp)
            
            set(temp_box, 'string', name_file_temp(1:(end - 4)))
            
            switch name_file_temp((end - 2):end) % read in temperature data
                
                case 'csv' % the old data file (much simpler to use .mat!)
                    
                    temp_all = dlmread([path_data name_file_temp], ',');
                    % clean up temperature data
                    temp_all = temp_all(setdiff(1:size(temp_all, 1), find(~temp_all(:, 1))), :); % removes pts where day of month == 0, i.e., glitchy pts
                    if (size(temp_all, 2) == 7) % correct for brief mistake in early version of code that didn't reproduce old .csv style (missing column decimal hours of current day)
                        disp('We can''t use the temperature .csv file for this run because of time ambiguity. Set temp_type to ''mat'' and run again.')
                        return
                    else
                        temp_sample ...
                            = temp_all(:, 5); % temperature of sample (B)
                        setpt ...
                            = temp_all(:, 8); % set point temperature
                    end
                    % convert temperature times to date numbers
                    time_start_vec ...
                            = datevec(time{1}(1));          % day of month of start of data acquisiton
                    if ((time_start_vec(3) == (temp_all(1, 1) + 1)) || ((time_start_vec(3) == 1) && (temp_all(1, 1) == 30)) || ((time_start_vec(3) == 1) && (temp_all(1, 1) == 31)))
                        time_start_num ...
                            = floor(time{1}(1)) - 1;        % subtract a day if it appears that temperature measurements started a day before data
                    else
                        time_start_num ...
                            = floor(time{1}(1));            % data and temperatures started on same day
                    end
                    time_temp ...
                            = time_start_num + (temp_all(:, 2) ./ 24);
                    ind_temp_change ...
                            = find(diff(temp_all(:, 1)));   % indices where date changes
                    for ii = 1:length(ind_temp_change)
                        time_temp((ind_temp_change(ii) + 1):end) ...
                            = time_temp((ind_temp_change(ii) + 1):end) + ii;
                    end
                    clear time_start_num time_start_vec ind_temp_change
                    
                case 'mat' % the new way
                    
                    tmp_load = load([path_temp name_file_temp]);
                    [time_temp, temp_sample, setpt] ...
                            = deal(tmp_load.time_num, tmp_load.temp_b, tmp_load.setpt);
                    clear tmp_load
                    [temp_sample, setpt] ...
                            = deal(temp_sample(logical(time_temp)), setpt(logical(time_temp))); % remove time glitches, which cause time_temp=0
                    time_temp ...
                            = time_temp(logical(time_temp)); % only keep times that are not 0
                    
            end
            
            set(status_box, 'string', 'Temperature file loaded')
        end
        
    end
        
    function do_merge(source, eventdata) % sort through temperature and sweep data, find temperature-stable sweeps
        
        % first remove old plots if they exist
        if ishandle(lg)
            delete(lg)
        end
        if any(ishandle(pd{1}))
            delete(pd{:})
        end
        if ishandle(pt)
            delete(pt)
        end
        if ishandle(pis(1))
            delete([pis fl ph])
        end

%         gap_guard           = 1e-3 * str2double(get(gap_guard_edit, 'string'));
%         rad_elec            = 1e-3 * str2double(get(rad_elec_edit, 'string'));
        threshold_std       = str2double(get(threshold_std_edit, 'string')); % std. dev. of temp. during sweep, below which is stable 
        pts_stable_min      = str2double(get(pts_stable_min_edit, 'string')); % numbers of pts that had stable temp. (not so important)
    
%         % correct for guard gap
%         d_t                 = ((4 * thick) / (pi * gap_guard)) * log(cosh((pi * gap_guard) / (4 * thick)));
%         area_corr           = (gap_guard / rad_elec) * (1 + (rad_elec / gap_guard) - d_t); % should be about 0.85
%         area_elec           = area_corr * area_elec;
        
        % load open/short data
        [stack_o, num_pts_o] = read_meta([path_os 'open' path_fix], name_os, 'os');
        [~, freq_os, ~, ~, ~, ~, Z_comp_o] ...
                            = read_solartron_csv([path_os 'open' path_fix], name_os, stack_o(1), num_pts_o, 1);
        [stack_s, num_pts_s] = read_meta([path_os 'short' path_fix], name_os, 'os');
        [~, ~, ~, ~, Z_real_s, Z_imag_s, Z_comp_s] ...
                            = read_solartron_csv([path_os 'short' path_fix], name_os, stack_s(1), num_pts_s, 1);
        [freq_os, Z_comp_o, Z_real_s, Z_imag_s, Z_comp_s] ...
                            = deal(freq_os{1}, Z_comp_o{1}, Z_real_s{1}, Z_imag_s{1}, Z_comp_s{1}); % get out of cells
        
        % open/short impedances
        freq_rad_os         = (2 * pi) .* freq_os;
        Z_imag_s            = Z_imag_s ./ freq_rad_os;
        if (stack_s(1) > 1)
            Z_real_s        = mean(Z_real_s, 2);
            Z_imag_s        = mean(Z_imag_s, 2);
        end
        C_os                = real((1 ./ ((Z_comp_o - Z_comp_s))) ./ (1i .* freq_rad_os)); % capacitance
        resist_os           = real(1 ./ ((1 ./ (Z_comp_o - Z_comp_s)) - (1i .* freq_rad_os .* C_os))); % resistance
        if (stack_s(1) > 1)
            C_os            = mean(C_os, 2);
            resist_os       = mean(resist_os, 2);
        end
        
        % trim open/short values
        [Z_real_s, Z_imag_s, resist_os, C_os] ...
                            = deal(Z_real_s(1:num_pts), Z_imag_s(1:num_pts), resist_os(1:num_pts), C_os(1:num_pts));
        
        [Z_corr, C_abs, C_phase, C_real, C_imag, C_loss_tan, permitt_real, permitt_imag, loss_tan, resist_abs, resist_phase, resist_real, resist_imag, conduct] ...
                            = deal(cell(num_files, 1));
        
        % extract dielectric measurements of interest from raw data
        for ii = 1:num_files
            % correct measurements for open/short values
            freq_rad        = (2 * pi) .* freq{ii}; % radial frequency
            Z_corr{ii}      = 1 ./ ((1 ./ (Z_comp{ii} - Z_real_s - (1i .* freq_rad .* Z_imag_s))) - (1 ./ resist_os) - (1i .* freq_rad .* C_os));
            % derived properties from complex impedance
            C_abs{ii}       = 1 ./ ((2 * pi) .* freq{ii} .* abs(Z_corr{ii}));
            C_phase{ii}     = (pi / 2) + angle(Z_corr{ii});
            C_real{ii}      = C_abs{ii} .* cos(C_phase{ii});
            C_imag{ii}      = C_abs{ii} .* sin(C_phase{ii});
            C_loss_tan{ii}  = C_imag{ii} ./ C_real{ii};
            permitt_real{ii}= C_real{ii} .* (thick / (area_elec * permitt_vacuum));
            permitt_imag{ii}= C_imag{ii} .* (thick / (area_elec * permitt_vacuum));
            permitt_real{ii}(permitt_real{ii} < 0) ...
                            = NaN;
            permitt_imag{ii}(permitt_imag{ii} < 0) ...
                            = NaN; % nip zaniness in the bud
            permitt_real{ii}= permitt_real{ii} ./ permitt_norm; % calibrate!
            permitt_imag{ii}= permitt_imag{ii} ./ permitt_norm;
            loss_tan{ii}    = permitt_imag{ii} ./ permitt_real{ii};
            resist_abs{ii}  = abs(Z_corr{ii}) .* (area_elec / thick);
            resist_abs{ii}(resist_abs{ii} < 0) ...
                            = NaN; % nip zaniness in the bud
            resist_phase{ii}= angle(Z_corr{ii});
            resist_real{ii} = resist_abs{ii} .* cos(resist_phase{ii});
            resist_imag{ii} = -resist_abs{ii} .* sin(resist_phase{ii});
            conduct{ii}     = freq_rad .* permitt_imag{ii} .* permitt_vacuum;
        end
        
        % trim temperature data beyond a certain time because they are bad
        if ~isempty(time_ignore)
            [temp_sample, setpt] ...
                            = deal(temp_sample(time_temp <= time_ignore), setpt(time_temp <= time_ignore));
            time_temp       = time_temp(time_temp <= time_ignore);
        end
        
        % assign temperatures to appropriate temperature set points
        setpt_fix           = unique(setpt(~diff(setpt))); % all the set points, sorted by unique from lowest to highest (not chronological order)
        num_setpt           = length(setpt_fix);
        ind_setpt_start     = zeros(num_setpt, 1);
        for ii = 1:length(setpt_fix)
            ind_setpt_start(ii) ...
                            = find((setpt == setpt_fix(ii)), 1);    % first index at each set point
        end
        [~, ind_setpt_ord]  = sort(ind_setpt_start);                % sort indices into chronological order
        setpt_fix           = setpt_fix(ind_setpt_ord);             % return set points to chronological order
        clear temp_all ind_setpt_ord
        
        % find indices where set point was fixed and where sample temperature was stable
        [ind_stable_start, ind_stable_end] ...
                            = deal(NaN(num_setpt, 1));
        for ii = 1:num_setpt
            ind_setpt_curr  = find(setpt == setpt_fix(ii)); % indices relevant to current set point temperature
            std_temp        = zeros((length(ind_setpt_curr) - pts_stable_min), 1);
            for jj = 1:(length(ind_setpt_curr) - pts_stable_min);
                std_temp(jj) ...
                            = std(temp_sample(ind_setpt_curr(jj):ind_setpt_curr(end))); % standard deviation of temperature with progressively fewer pts
            end
            if ~isempty(find((std_temp < threshold_std), 1)) % assign value if standard deviation was below threshold
                ind_stable_start(ii) ...
                            = ind_setpt_curr(find((std_temp < threshold_std), 1)); % index where standard deviation is first below threshold
                ind_stable_end(ii) ...
                            = ind_setpt_curr(end); % index above
            end
        end
        
        ind_setpt_stable    = find(~isnan(ind_stable_start));
        clear ind_setpt_curr std_temp jj
        
        % temperature statistics for each file/sweep, and determine if temperature was stable and match to individual set points
        [temp_mean, temp_std, temp_min, temp_max, temp_match] ...
                            = deal(NaN(num_files, 1));
        sweep_stable        = false(num_files, 1); % logical test initialization
        for ii = 1:num_files
            temp_curr       = temp_sample((time_temp >= time{ii}(1)) & (time_temp <= (time{ii}(end - 1) + ((time{ii}(end) - time{ii}(end - 1)) / 2))));
            if ~isempty(temp_curr) % only do if there are recorded temperatures during sweep!
                temp_mean(ii) ...
                            = mean(temp_curr);
                temp_std(ii) ...
                            = std(temp_curr);
                temp_min(ii) ...
                            = min(temp_curr);
                temp_max(ii) ...
                            = max(temp_curr);
                % determine if this sweep occurred during a temperature-stable period
                if any((time{ii}(1) >= time_temp(ind_stable_start(ind_setpt_stable))) & ((time{ii}(end - 1) + ((time{ii}(end) - time{ii}(end - 1)) / 2)) <= time_temp(ind_stable_end(ind_setpt_stable))))
                    sweep_stable(ii) ...
                            = true; % current file/sweep was temperature-stable
                    ind_match ...
                            = (time{ii}(1) >= time_temp(ind_stable_start(ind_setpt_stable))) & ((time{ii}(end - 1) + ((time{ii}(end) - time{ii}(end - 1)) / 2)) <= time_temp(ind_stable_end(ind_setpt_stable)));
                    temp_match(ii) ...
                            = setpt_fix(ind_setpt_stable(ind_match)); % set point temperature during this sweep
                end
            end
        end
        clear temp_curr ind_match
        
        % choose last of sweeps with each unique set point as best to use
        temp_match_unique   = unique(temp_match); % unique set points
        temp_match_best     = false(num_files, 1);
        for ii = 1:length(temp_match_unique)
            temp_match_best(find((temp_match_unique(ii) == temp_match), 1, 'last')) ...
                            = true;
        end
        clear ii temp_match_unique;
        
        % temperature matching verification
        if (num_files > 1)
            colors         = colormap(jet(double(num_files)));
        else
            colors         = [0 0 1];
        end
        xticks              = (time_temp(1) - mod(time_temp(1), (30 / (24 * 60)))):(30 / (24 * 60)):time_temp(end);
        axes(tempax) %#ok<*MAXES>
        pt                  = plot(time_temp, temp_sample, 'k', 'linewidth', 2); % sample temperature profile (black)
        pis                 = zeros(1, length(ind_setpt_stable));
        for ii = 1:length(ind_setpt_stable) % portions of sample temperature profile deemed stable (gray)
            pis(ii)         = plot(time_temp(ind_stable_start(ind_setpt_stable(ii)):ind_stable_end(ind_setpt_stable(ii))), ...
                                   temp_sample(ind_stable_start(ind_setpt_stable(ii)):ind_stable_end(ind_setpt_stable(ii))), 'color', [0.7 0.7 0.7], 'linewidth', 2);
        end
        [fl, ph]            = deal(zeros(1, num_files));
        for ii = 1:num_files % mean temperatures for each sweep (solid rainbow) and +/- box for standard deviation (default solid rainbow)
            fl(ii)          = fill([repmat(time{ii}(1), 1, 2) repmat((time{ii}(end - 1) + ((time{ii}(end) - time{ii}(end - 1)) / 2)), 1, 2) time{ii}(1)], ...
                                   [(temp_mean(ii) - temp_std(ii)) repmat((temp_mean(ii) + temp_std(ii)), 1, 2) repmat((temp_mean(ii) - temp_std(ii)), 1, 2)], ...
                                   'w', 'facecolor', 'none', 'edgecolor', colors(ii, :), 'linewidth', 2); % box for standard deviation
            if ~sweep_stable(ii) % temperature-unstable sweeps are dashed
                set(fl(ii), 'linestyle', '--')
            end
            ph(ii)          = plot([time{ii}(1) ((time{ii}(end - 1) + ((time{ii}(end) - time{ii}(end - 1)) / 2)))], repmat(temp_mean(ii), 1, 2), 'color', colors(ii, :), 'linewidth', 2); ...
        end
        axis([time_temp(1) time_temp(end) min(temp_sample) max(temp_sample)])
        set(gca, 'xtick', xticks)
        datetick('x', 15, 'keepticks', 'keeplimits')
        
        good_temp           = find(temp_match_best);
        if ~isempty(good_temp) % only do if there are any good temp
            for ii = 1:4
                axes(ax(ii)) %#ok<*LAXES>
                pd{ii}      = NaN(1, length(good_temp));
                for jj = 1:length(good_temp)
                    if (ii ~= 4) % all plots except phase should be y log
                        pd{ii}(jj) ...
                            = loglog(freq{good_temp(jj)}, eval(plots{ii}), 'linewidth', 2, 'color', colors(good_temp(jj), :));
                        set(gca, 'yscale', 'log')
                    else
                        pd{ii}(jj) ...
                            = semilogx(freq{good_temp(jj)}, eval(plots{ii}), 'linewidth', 2, 'color', colors(good_temp(jj), :));
                    end
                end
                xlim([freq{1}(end) freq{1}(1)])
                if (ii == 3)
                    lg      = legend(num2str(round(temp_mean(good_temp))));
                end
            end
            set(status_box, 'string', [num2str(length(good_temp)) ' stable sweep(s)'])
        else
            set(status_box, 'string', 'No stable sweeps!')
        end
        
    end
    
    function do_reset_all(source, eventdata) % reset everything
       do_reset_data
       do_reset_temp
%         gap_guard           = 0.001;
%         rad_elec            = 0.012;
        threshold_std       = 0.04;
        pts_stable_min      = 10;
        time_ignore         = [];
        time_ignore_str     = '';        
%         set(gap_guard_edit, 'string', sprintf('%3.1f', gap_guard));
%         set(rad_elec_edit, 'string', sprintf('%3.1f', rad_elec));
        set(threshold_std_edit, 'string', sprintf('%6.4f', threshold_std))
        set(pts_stable_min_edit, 'string', sprintf('%d', pts_stable_min))
        set(ignore_box, 'string', '')
        [d_t, area_corr, stack_o, num_pts_o, freq_os, Z_comp_o, stack_s, num_pts_s, Z_real_s, Z_imag_s, Z_comp_s, ...
            freq_rad_os, C_os, resist_os, Z_corr, C_abs, C_phase, C_real, C_imag, C_loss_tan, permitt_real, permitt_imag, loss_tan, resist_abs, resist_phase, resist_real, resist_imag, freq_rad, setpt_fix, num_setpt, ...
            ind_setpt_start, ind_setpt_ord, ind_stable_start, ind_stable_end, ind_setpt_curr, std_temp, ind_setpt_stable, temp_mean, temp_std, temp_min, temp_max, temp_match, sweep_stable, temp_curr, ind_match, ...
            temp_match_unique, temp_match_best, colors, xticks, good_temp, conduct] ...
                            = deal(NaN); %#ok<SETNU>
        set(status_box, 'string', 'All cleared')
    end

    function do_reset_data(source, eventdata) % reset data
        if ishandle(lg)
            delete(lg)
        end
        if any(ishandle(pd{1}))
            delete(pd{:})
        end
        if any(ishandle(pis(1)))
            delete([pis fl ph])
        end
        [name_file_data, name, name_os] ...
                            = deal('');
        set(solar_box, 'string', '');
        [thick, area_elec, num_files, num_files_tot, num_loops, freq_max, freq_min, num_pts_decade, stack, num_pts_stack, num_pts, time, freq, Z_comp, dia_elec, dia_sample, dia_elec_ind, dia_sample_ind, permitt_norm] ...
                            = deal(NaN); %#ok<SETNU>
        for ii = 1:4
            axes(ax(ii))
            axis([freq_default ylim_default(ii, :)]) % reset frequency axis to default
        end
        set(status_box, 'string', 'Data cleared.')
    end

    function do_reset_temp(source, eventdata) % reset temperature data
        if ishandle(lg)
            delete(lg)
        end
        if any(ishandle(pd{1}))
            delete(pd{:})
        end
        if ishandle(pt)
            delete(pt)
        end
        if ishandle(pis(1))
            delete([pis fl ph])
        end
        [name_file_temp, path_temp, name_file_save, path_save] ...
                            = deal('');
        set(temp_box, 'string', '')
        [temp_all, temp_sample, setpt, time_start_vec, time_start_num, time_temp, ind_temp_change, tmp_load] ...
                            = deal(NaN);
        axes(tempax)
        axis([0 1 0 1])
        set(status_box, 'string', 'Temp. data cleared.')
    end

    function do_ignore(source, eventdata)
        axes(tempax)
        [time_ignore, ~]    = ginput(1);
        time_ignore_str     = datestr(time_ignore);
        if (length(time_ignore) > 1)
            time_ignore     = 0;
            set(ignore_box, 'string', 'only click once!')
        end
        set(ignore_box, 'string', time_ignore_str(13:end))
    end

    function clear_ignore(source, eventdata) % reset ignore box
        time_ignore         = [];
        time_ignore_str     = '';
        set(ignore_box, 'string', '')
    end

    function do_save(source, eventdata) % save everything (even stuff that is not really needed, but it's hard to keep it out)
        [name_file_save, path_save] ...
                            = uiputfile('*.mat', 'Save merged file:', [path_data name_file_data '_proc.mat']);
        save([path_save name_file_save])
    end

end
function merge_solartron_gui
% MERGE_SOLARTRON_GUI Merge Solartron data recorded at different
% temperatures, based on temperature matching with each sweep and
% comparison between different sweeps at the same temperature.
% 
% See dielectric_spectroscopy_gui_man.pdf for operation.
%
% David Stillman, Joe MacGregor
% Last updated: 09/22/15

do_big                      = false; % if true, uses 5 columns instead of 2

% initialize loads of variables
[path_curr, file_save, leg_str, vars2trim_str, var_merge, curr_plot, name_file_cfg] ...
                            = deal('');
[num_load, ind_val, ind_merge, num_temp, num_merge, tmp1, tmp2, tmp3, tmp4] ...
                            = deal(0);
[name_file_tmp, path2load, file2load, file2load_ref, vars_all, num_set] ...
                            = deal({});
[curr_runs, curr_inds, temp_all] ...
                            = deal([]);
curr_box                    = 1;

[C_abs, C_imag, C_loss_tan, C_phase, C_real, Z_comp, Z_corr, freq, loss_tan, permitt_imag, permitt_real, resist_abs, resist_phase, resist_real, conduct, temp_match, temp_max, temp_mean, temp_min, temp_std, time] ...
                            = deal(0);
% variables to be trimmed according to each file's temp_match_best
vars2trim                   = {'C_abs' 'C_imag' 'C_loss_tan' 'C_phase' 'C_real' 'Z_comp' 'Z_corr' 'freq' 'loss_tan' 'conduct' 'permitt_imag' 'permitt_real' 'resist_abs' 'resist_phase' 'resist_real' 'temp_match' 'temp_max' 'temp_mean' 'temp_min' 'temp_std' 'time'};
num_trim                    = length(vars2trim);

% initialize GUI stuff
if do_big
    subplot_start           = [0.05 0.65 0.43 0.17]; % start position of 1st subplot    
else
    subplot_start           = [0.05 0.48 0.35 0.36]; % start position of 1st subplot
end
plots                       = {'permitt_real{jj}' 'conduct{jj}' 'permitt_imag{jj}' '(resist_phase{jj} .* (180 / pi))'};
ylabels                     = {'\epsilon''' '\sigma (S m^{-1})' '\epsilon''''' '\theta (\circ)'}; % ylabels for each supblot
plots_single                = {'vars_all{curr_runs(jj)}.permitt_real{curr_inds(jj)}' 'vars_all{curr_runs(jj)}.conduct{curr_inds(jj)}' 'vars_all{curr_runs(jj)}.permitt_imag{curr_inds(jj)}' '(vars_all{curr_runs(jj)}.resist_phase{curr_inds(jj)} .* (180 / pi))'};
[ax, ax2]                   = deal(NaN(1, 4));
freq_default                = [1e-3 1e7]; % freq range to display
ylim_default                = [1 1e4; 1e-12 1e-2; 1e-4 1e4; -90 0]; % ylims for each subplot
[colors, lg, lg2]           = deal(NaN);
pd                          = cell(1, 4);
if do_big
    num_box                 = 100;
else
    num_box                 = 40; % number of boxes for different temperatures (effectively limits number of different temperatures)
end
[temp_box, ind_box]         = deal(NaN(1, num_box));

%% draw the gui
merge_gui                   = figure('toolbar', 'figure', 'name', 'MERGE_SOLARTRON_GUI', 'position', [1 1 1600 1200]);
if do_big
    for ii = 1:4 %#ok<*FXUP>
        ax(ii)              = subplot('position', (subplot_start + [0 (-0.19 * (ii - 1)) 0 0]));
        hold on
        set(gca, 'fontsize', 18, 'xscale', 'log', 'color', [0.9 0.9 0.9])
        if (ii < 4) % all plots except phase should be y log
            set(gca, 'yscale', 'log')
        end
        axis([freq_default ylim_default(ii, :)])
        if (ii == 4)
            xlabel('Frequency (Hz)')
        else
            set(gca, 'xticklabel', {})
        end
        ylabel(ylabels{ii})
        grid on
        box on
    end
else
    for ii = 1:4
        subplot_curr        = subplot_start;
        if any(ii == [2 4])
            subplot_curr    = subplot_curr + [0.40 0 0 0];
        end
        if (ii > 2)
            subplot_curr    = subplot_curr + [0 -0.39 0 0];
        end
        ax(ii)              = subplot('position', subplot_curr);
        hold on
        set(gca, 'fontsize', 16, 'xscale', 'log', 'color', [0.9 0.9 0.9])
        if (ii < 4) % all plots except phase should be y log
            set(gca, 'yscale', 'log')
        end
        axis([freq_default ylim_default(ii, :)])
        if (ii > 2)
            xlabel('Frequency (Hz)');
        else
            set(gca, 'xticklabel', {});
        end
        ylabel(ylabels{ii})
        grid on
        box on
    end
end
linkaxes(ax, 'x')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Load data', 'units', 'normalized', 'position', [0.01 0.95 0.1 0.04], 'callback', @load_data, 'fontsize', 16, 'foregroundcolor', 'k')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Load config.', 'units', 'normalized', 'position', [0.01 0.89 0.1 0.04], 'callback', @load_cfg, 'fontsize', 16, 'foregroundcolor', 'k')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Plot single', 'units', 'normalized', 'position', [0.11 0.95 0.1 0.04], 'callback', @plot_single, 'fontsize', 16, 'foregroundcolor', 'm')
check_box                   = uicontrol(merge_gui, 'style', 'popupmenu', 'string', 'N/A', 'units', 'normalized', 'position', [0.21 0.96 0.07 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Plot merged', 'units', 'normalized', 'position', [0.11 0.89 0.1 0.04], 'callback', @plot_merge, 'fontsize', 16, 'foregroundcolor', 'm')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Pop plot', 'units', 'normalized', 'position', [0.21 0.89 0.07 0.04], 'callback', @plot_pop, 'fontsize', 16, 'foregroundcolor', 'm')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Save config.', 'units', 'normalized', 'position', [0.28 0.95 0.1 0.04], 'callback', @save_cfg, 'fontsize', 16, 'foregroundcolor', 'g')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Save', 'units', 'normalized', 'position', [0.28 0.89 0.1 0.04], 'callback', @do_save, 'fontsize', 16, 'foregroundcolor', 'g')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Merge', 'units', 'normalized', 'position', [0.38 0.95 0.1 0.04], 'callback', @do_merge, 'fontsize', 16, 'foregroundcolor', 'b')
uicontrol(merge_gui, 'style', 'pushbutton', 'string', 'Reset', 'units', 'normalized', 'position', [0.38 0.89 0.1 0.04], 'callback', @do_reset, 'fontsize', 16, 'foregroundcolor', 'r')
if do_big
    list_box                = uicontrol(merge_gui, 'style', 'popupmenu', 'string', 'N/A', 'units', 'normalized', 'position', [0.01 0.84 0.47 0.04], 'fontsize', 16, 'foregroundcolor', 'k');    
    for ii = 1:20 % 1st column of temps/checkboxes
        temp_box(ii)        = annotation('textbox', [0.50 (0.96 - ((ii - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.55 (0.96 - ((ii - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
    for ii = 21:40 % 2nd column
        temp_box(ii)        = annotation('textbox', [0.60 (0.96 - (((ii - 20) - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.65 (0.96 - (((ii - 20) - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
    for ii = 41:60 % 3rd column
        temp_box(ii)        = annotation('textbox', [0.70 (0.96 - (((ii - 40) - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.75 (0.96 - (((ii - 40) - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
    for ii = 61:80 % 4th column
        temp_box(ii)        = annotation('textbox', [0.80 (0.96 - (((ii - 60) - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.85 (0.96 - (((ii - 60) - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
    for ii = 81:100 % 5th column
        temp_box(ii)        = annotation('textbox', [0.90 (0.96 - (((ii - 80) - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.95 (0.96 - (((ii - 80) - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
else
    annotation('textbox', [0.48 0.95 0.18 0.05], 'string', 'File reference list', 'fontsize', 16, 'edgecolor', 'none')
    list_box                = uicontrol(merge_gui, 'style', 'popupmenu', 'string', 'N/A', 'units', 'normalized', 'position', [0.48 0.92 0.29 0.04], 'fontsize', 16, 'foregroundcolor', 'k');
    for ii = 1:(num_box / 2) % 1st column of temps/checkboxes
        temp_box(ii)        = annotation('textbox', [0.81 (0.96 - ((ii - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.86 (0.96 - ((ii - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
    for ii = ((num_box / 2) + 1):num_box % 2nd column
        temp_box(ii)        = annotation('textbox', [0.91 (0.96 - (((ii - (num_box / 2)) - 1) * 0.05)) 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
        ind_box(ii)         = uicontrol(merge_gui, 'style', 'popupmenu', 'string', ' ', 'units', 'normalized', 'position', [0.96 (0.96 - (((ii - (num_box / 2)) - 1) * 0.05)) 0.04 0.03], 'fontsize', 16, 'foregroundcolor', 'k');
    end
end

%% sub-functions
    
    function load_data(source, eventdata) % load data into workspace
        
        if ~isempty(path_curr)
            [name_file_tmp, path_curr] ...
                            = uigetfile('*.mat', 'Load data/temperature Matlab files:', 'multiselect', 'on', path_curr); % dialog box for selecting multiple files in same directory
        else
            [name_file_tmp, path_curr] ...
                            = uigetfile('*.mat', 'Load data/temperature Matlab files:', 'multiselect', 'on');
        end
        
        if ~iscell(name_file_tmp)
            if ~name_file_tmp
                [name_file_tmp, path_curr] ...
                            = deal('');
            end
        end
        
        if isempty(name_file_tmp) % only process if didn't cancel
            return
        end

        if ischar(name_file_tmp) % only selected one file
            name_file_tmp   = {name_file_tmp};
        end
        name_file_tmp       = fliplr(name_file_tmp); % flip order of files because negative is usually dropped in name
        file2load           = [file2load name_file_tmp];
        path2load           = [path2load repmat({path_curr}, 1, length(name_file_tmp))];
        num_load            = length(file2load);
        
        % cycle through each new file (no check for repeated files so watch out)
        for ii = 1:length(name_file_tmp)
            file2load_ref   = [file2load_ref [num2str(curr_box) ': ' name_file_tmp{ii}(1:(end - 4))]]; %#ok<*AGROW>
            tmp1            = load([path_curr name_file_tmp{ii}]);
            vars_all        = [vars_all tmp1]; % load each file separately and immediately rename
            if ~isempty(strfind(name_file_tmp{ii}, 'proc'))
                for jj = 1:num_trim
                    eval(['vars_all{curr_box}.' vars2trim{jj} ' = vars_all{curr_box}.' vars2trim{jj} '(vars_all{curr_box}.temp_match_best);'])
                end
            end
            temp_all        = [temp_all; vars_all{curr_box}.temp_match]; % add this file's match temperatures to set
            curr_box        = curr_box + 1;
        end
        temp_all            = unique(temp_all); % keep unique set points
        num_temp            = length(temp_all);
        
        % assign values to boxes and only include values that make sense for any given match temperature
        num_set             = cell(1, num_temp); % set of indices to files for each temperature
        for ii = 1:num_temp
            set(temp_box(ii), 'string', sprintf('%4.1f', temp_all(ii)))
            num_set{ii}     = 0;
            for jj = 1:num_load % test each merge temperature to see if it exists in files
                if any(temp_all(ii) == vars_all{jj}.temp_match)
                    num_set{ii} ...
                            = [num_set{ii} jj]; % only add index if there was a match
                end
            end
            set(ind_box(ii), 'string', num2cell(num_set{ii}), 'value', 2) % load set and 
            if (length(num_set{ii}) > 2)
                set(temp_box(ii), 'color', 'r') % make temperatures with more then one option red (as in red bad, choice must be made)
            end
        end
        
        set(list_box, 'string', file2load_ref) % update file list popup menu
        set(check_box, 'string', num2cell(temp_all)) % for single temperature plotting

    end
%%
    function do_merge(source, eventdata)
        [C_abs, C_imag, C_loss_tan, C_phase, C_real, Z_comp, Z_corr, freq, loss_tan, permitt_imag, permitt_real, resist_abs, resist_phase, resist_real, temp_match, temp_max, temp_mean, temp_min, temp_std, time, conduct] ...
                            = deal(0); % reset for merging
        % prepare merging
        ind_merge           = zeros(num_temp, 1);
        for ii = 1:num_temp
            ind_merge(ii)   = num_set{ii}(get(ind_box(ii), 'value')); % get index based on where it is in num_set
        end
        num_merge           = length(ind_merge(logical(ind_merge))); % number of temperatures to merge is equal to the number of non-zero file choices
        
        % merge runs
        for ii = 1:num_trim
            var_merge       = '';
            for jj = 1:num_temp % make long string that sorts out where each cell element from each run belongs
                if ind_merge(jj)
                    var_merge ...
                            = [var_merge 'vars_all{' num2str(ind_merge(jj)) '}.' vars2trim{ii} '(vars_all{' num2str(ind_merge(jj)) '}.temp_match == ' num2str(temp_all(jj)) '); '];
                end
            end
            eval([vars2trim{ii} ' = [' var_merge '];']) % assign variable without suffix for concatenated mess
        end
        plot_merge % plot results after each merge
    end
%%
    function plot_single(source, eventdata) %#ok<*INUSD> % individual temperature QC
        temp_check          = temp_all(get(check_box, 'value')); % temperature to compare/check
        if any(ishandle(pd{1})) % cleanup if necessary
            delete(pd{:})
            [curr_runs, curr_inds] ...
                            = deal([]);
        end
        if ishandle(lg)
            delete(lg)
            leg_str         = '';
        end
        leg_str             = {};
        for ii = 1:num_load
            tmp2            = find(vars_all{ii}.temp_match == temp_check);
            if ~isempty(tmp2)
                curr_runs   = [curr_runs repmat(ii, 1, length(tmp2))]; % if that temperature is present in that run, then include
                curr_inds   = [curr_inds tmp2];
            end
        end
        for ii = 1:length(curr_runs)
            leg_str{ii}     = sprintf('%s: %4.3f +/- %1.3f', file2load{curr_runs(ii)}(1:(end - 9)), vars_all{curr_runs(ii)}.temp_mean(curr_inds(ii)), vars_all{curr_runs(ii)}.temp_std(curr_inds(ii)));
        end
        if (length(curr_runs) > 1)
            colors          = colormap(jet(length(curr_runs)));
        else
            colors          = [0 0 1]; % just make it blue if there's only one color
        end
        for ii = 1:4
            axes(ax(ii)) %#ok<*LAXES>
            pd{ii}          = zeros(1, length(curr_runs));
            for jj = 1:length(curr_runs)
                if (ii ~= 4) % all plots except phase should be y log
                    pd{ii}(jj) ...
                            = loglog(vars_all{curr_runs(jj)}.freq{curr_inds(jj)}, eval(plots_single{ii}), 'linewidth', 2, 'color', colors(jj, :));
                else
                    pd{ii}(jj) ...
                            = semilogx(vars_all{curr_runs(jj)}.freq{curr_inds(jj)}, eval(plots_single{ii}), 'linewidth', 2, 'color', colors(jj, :));
                end
            end
            xlim([vars_all{jj}.freq{1}(end) vars_all{jj}.freq{1}(1)])
            if (ii == 3)
                lg          = legend(leg_str);
                set(lg, 'interpreter', 'none')
            end
        end
        curr_plot           = 'single';
    end
%%
    function plot_merge(source, eventdata)
        if any(ishandle(pd{1})) % cleanup
            delete(pd{:})
        end
        if ishandle(lg)
            delete(lg)
            leg_str         = '';
        end
        if (num_merge > 1);
            colors          = colormap(jet(num_merge));
        else
            colors          = [0 0 1];
        end
        leg_str             = cell(1, num_merge);
        for ii = 1:num_merge
            leg_str{ii}     = sprintf('%3.2f ', temp_mean(ii));
        end
        for ii = 1:4
            pd{ii}          = zeros(1, num_merge);
            axes(ax(ii))
            for jj = 1:num_merge
                if (ii ~= 4) % all plots except phase should be y log
                    pd{ii}(jj) ...
                            = loglog(freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(jj, :));
                else
                    pd{ii}(jj) ...
                            = semilogx(freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(jj, :));
                end
            end
            xlim([freq{1}(end) freq{1}(1)])
            if (ii == 3)
                lg          = legend(leg_str);
            end
        end
        curr_plot           = 'merge';
    end
%%
    function plot_pop(source, eventdata)
        switch curr_plot
            case 'single'
                figure('name', 'MERGE SINGLE')
                for ii = 1:4
                    ax2(ii) = subplot(2, 2, ii);
                    hold on
                    set(gca, 'fontsize', 16, 'xscale', 'log', 'color', [0.9 0.9 0.9])
                    for jj = 1:length(curr_runs)
                        if (ii ~= 4) % all plots except phase should be y log
                            loglog(vars_all{curr_runs(jj)}.freq{curr_inds(jj)}, eval(plots_single{ii}), 'linewidth', 2, 'color', colors(jj, :))
                            set(gca, 'yscale', 'log')
                        else
                            semilogx(vars_all{curr_runs(jj)}.freq{curr_inds(jj)}, eval(plots_single{ii}), 'linewidth', 2, 'color', colors(jj, :))
                        end
                    end
                    xlim([vars_all{jj}.freq{1}(end) vars_all{jj}.freq{1}(1)])
                    axis([freq_default ylim_default(ii, :)])
                    xlabel('Frequency (Hz)')
                    ylabel(ylabels{ii})
                    if (ii == 3)
                        lg2 = legend(leg_str);
                        set(lg2, 'interpreter', 'none')
                    end
                    grid on
                    box on
                end
                linkaxes(ax2, 'x')
            case 'merge'
                figure('name', 'MERGE ALL')
                for ii = 1:4
                    ax2(ii)  = subplot(2, 2, ii);
                    hold on
                    set(gca, 'fontsize', 16, 'xscale', 'log', 'color', [0.9 0.9 0.9])
                    for jj = 1:num_merge
                        if (ii ~= 4) % all plots except phase should be y log
                            loglog(freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(jj, :))
                            set(gca, 'yscale', 'log')
                        else
                            semilogx(freq{jj}, eval(plots{ii}), 'linewidth', 2, 'color', colors(jj, :))
                        end
                    end
                    xlim([freq{1}(end) freq{1}(1)])
                    xlabel('Frequency (Hz)')
                    ylabel(ylabels{ii})
                    if (ii == 3)
                        legend(leg_str)
                    end
                    grid on
                    box on
                end
                linkaxes(ax2, 'x')
        end
    end
%%
    function do_save(source, eventdata) % save trimmed set
        save_cfg % first save config, if config filename is normal, then no need to reenter name of file
        if ~strcmp(name_file_cfg((end - 13):end), '_merge_cfg.mat')
            file_save       = uiputfile('*.mat', 'Save merged data:', '_merge.mat');
        else
            file_save       = name_file_cfg(1:(end - 8));
        end
        for ii = 1:num_trim
            vars2trim_str   = [vars2trim_str vars2trim{ii} ' '];
        end
        eval(['save ' path_curr file_save ' ' vars2trim_str ';'])
    end
%%
    function do_reset(source, eventdata) % reset everything
        if ishandle(lg)
            delete(lg)
        end
        if any(ishandle(pd{1}))
            delete(pd{:})
        end
        [file_save, leg_str, vars2trim_str, var_merge, curr_plot, name_file_cfg] ...
                            = deal('');
        [num_load, ind_merge, num_temp, num_merge, tmp1, tmp2, tmp3, tmp4] ...
                            = deal(0);
        [name_file_tmp, path2load, file2load, file2load_ref, vars_all, num_set] ...
                            = deal({});
        [curr_runs, curr_inds, temp_all] ...
                            = deal([]);
        curr_box            = 1;
        [C_abs, C_imag, C_loss_tan, C_phase, C_real, Z_comp, Z_corr, freq, loss_tan, permitt_imag, permitt_real, resist_abs, resist_phase, resist_real, conduct, temp_match, temp_max, temp_mean, temp_min, temp_std, time] ...
                            = deal(0); %#ok<*SETNU>
        set(check_box, 'string', 'N/A', 'value', 1)
        set(list_box, 'string', 'N/A', 'value', 1)
        for ii = 1:num_box
            set(ind_box(ii), 'string', ' ', 'value', 1)
            set(temp_box(ii), 'string', ' ', 'color', 'k')
        end
    end
%%
    function load_cfg(source, eventdata) % load configuration
        
        do_reset % reset first
        
        if ~isempty(path_curr)
            [name_file_cfg, path_curr] ...
                            = uigetfile('*.mat', 'Load config. file (*_merge_cfg.mat):', path_curr); % dialog box for selecting config file
        else
            [name_file_cfg, path_curr] ...
                            = uigetfile('*.mat', 'Load config. file (*_merge_cfg.mat):');
        end
        
        if ~isempty(name_file_cfg)
            
            tmp3            = load([path_curr name_file_cfg]);
            [file2load, path2load, num_load, path_curr, temp_all, num_temp, ind_val, file_save, num_set] ...
                            = deal(tmp3.files2load, tmp3.paths2load, tmp3.num_load, tmp3.path_curr, tmp3.temp_all, tmp3.num_temp, tmp3.ind_val, tmp3.file_save, tmp3.num_set);
            
            for ii = 1:length(file2load)
                set(temp_box(ii), 'string', file2load{ii}(1:(end - 4)))
                file2load_ref ...
                            = [file2load_ref [num2str(curr_box) ': ' file2load{ii}(1:(end - 4))]]; %#ok<*AGROW>
                tmp         = load([path2load{ii} file2load{ii}]);
                vars_all    = [vars_all tmp]; % load each file separately and immediately rename
                if ~isempty(strfind(file2load{ii}, 'proc'))
                    for jj = 1:num_trim
                        eval(['vars_all{curr_box}.' vars2trim{jj} ' = vars_all{curr_box}.' vars2trim{jj} '(vars_all{curr_box}.temp_match_best);'])
                    end
                end
                curr_box    = curr_box + 1;
            end
            
            for ii = 1:num_temp
                set(temp_box(ii), 'string', sprintf('%4.1f', temp_all(ii)))
                set(ind_box(ii), 'string', num2cell(num_set{ii}), 'value', ind_val(ii))
                if (length(num_set{ii}) > 2)
                    set(temp_box(ii), 'color', 'r') % make temperatures with more then one option red (as in red requires attention; choice must be made)
                end
            end
            
            set(list_box, 'string', file2load_ref) % update file list popup menu
            set(check_box, 'string', num2cell(temp_all))
            
        end
        
    end
%%
    function save_cfg(source, eventdata) % save configuration
        if ~isempty(path_curr)
            [name_file_cfg, path_curr] ...
                            = uiputfile('*.mat', 'Save configuration:', [path_curr '_merge_cfg.mat']);
        else
            [name_file_cfg, path_curr] ...
                            = uiputfile('*.mat', 'Save configuration:', '_merge_cfg.mat');
        end
        ind_val             = zeros(num_temp, 1);
        for ii = 1:num_temp
            ind_val(ii)     = get(ind_box(ii), 'value');
        end
        save([path_curr name_file_cfg], 'file2load', 'path2load', 'num_load', 'path_curr', 'temp_all', 'num_temp', 'ind_val', 'file_save', 'num_set')
    end
%%
end
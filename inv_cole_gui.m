function inv_cole_gui
% INV_COLE_GUI Interactive Cole-Cole modeling of dielectric spectra.
% 
% INV_COLE_GUI fits the Cole-Cole parameters, DC conductivity and HF
% permittivity of Solartron data at a set of temperatures and for up to 4
% relaxations.
% 
% INV_COLE_GUI requires a license for the Optimization Toolbox.
% 
% See dielectric_spectroscopy_gui_man.pdf for operation.
% 
% David Stillman (SwRI-Boulder), Joe MacGregor (UTIG)
% Last updated: 09/22/15

if ~license('test', 'optimization_toolbox')
    error('inv_cole_gui:optimtoolbox', 'Optimization Toolbox must be licensed to invert dielectric spectra (uses LSQNONLIN).')
end

% adjust directory to search for reference arrhenius profiles
if ispc
    dir_arrhenius           = '';
else
    dir_arrhenius           = 'data/arrhenius/';
end

% initialize some fitting parameters
permitt_vacuum              = 8.8541878176e-12; % F/m
freq_repeat                 = [1 100];          % frequencies that get repeated in sweep, keep HF version
[freq_min, freq_max]        = deal(50e-2, 1e6); % max/min frequencies to analyze
num_relax                   = 2;                % number of relaxations
do_dc                       = false;            % default to not including DC conductivity
[do_weight, do_sep]         = deal(true);       % default to weighting frequencies (important for good inverse fits)
ind_remove                  = zeros(length(freq_repeat), 1); % indices with frequency range to remove because those frequencies are repeated
data_loaded                 = false;            % silly check for loaded data do then run forward model

% David's magic option set for good fits
options                     = optimset('display', 'off', 'maxiter', 5e4, 'maxfunevals', 8e2, 'largescale', 'on', 'algorithm', 'trust-region-reflective', 'tolx', 1e-25, 'tolfun', 1e-25, 'jacobian', 'off');

boltzmann                   = 1.3806503e-23; % Boltzmann's constant
conv_eV_J                   = 1.602176565e-19; % J/eV

% guess values to start with (for first temperature)
permitt_hf_guess            = 3;
conduct_dc_guess            = 1e-10;
%                             dielectric strength   relaxation freq.    Cole-Cole parameter (each row is a different relaxation)
relax_guess                 = [30                   5e2                 0.01;
                               4                    10                  0.05;
                               1                    0.1                 0.1;
                               1                    0.01                0.2];

% min/max values for sliders
permitt_hf_min              = 1; % dimensionless
permitt_hf_max              = 20;
conduct_dc_min              = 1e-12; % S/m
conduct_dc_max              = 1e-2;
permitt_diff_min            = 0; %#ok<*NASGU>
permitt_diff_max            = 500; % dimensionless
freq_relax_min              = 1e-3; % Hz
freq_relax_max              = 1e5;
alpha_cole_min              = 0; % dimensionless
alpha_cole_max              = 1;
tol                         = 0.95; % confidence interval tolerance for uncertainty calculations, e.g., 0.95 is 95%
frac_test                   = [0.25 0.75 2 5 10]; % range around model_final value about which to test (fraction of model_final value)
err_rel                     = [0.05 0.1]; % relative error of data [real imag]

% repeat mins/maxes for each parameter
[permitt_diff_min_1, permitt_diff_min_2, permitt_diff_min_3, permitt_diff_min_4] ...
                            = deal(permitt_diff_min); %#ok<ASGLU>
[permitt_diff_max_1, permitt_diff_max_2, permitt_diff_max_3, permitt_diff_max_4] ...
                            = deal(permitt_diff_max); %#ok<ASGLU>
[freq_relax_min_1, freq_relax_min_2, freq_relax_min_3, freq_relax_min_4] ...
                            = deal(freq_relax_min); %#ok<ASGLU>
[freq_relax_max_1, freq_relax_max_2, freq_relax_max_3, freq_relax_max_4] ...
                            = deal(freq_relax_max); %#ok<ASGLU>
[alpha_cole_min_1, alpha_cole_min_2, alpha_cole_min_3, alpha_cole_min_4] ...
                            = deal(alpha_cole_min); %#ok<ASGLU>
[alpha_cole_max_1, alpha_cole_max_2, alpha_cole_max_3, alpha_cole_max_4] ...
                            = deal(alpha_cole_max); %#ok<ASGLU>

% min/max values for model that bound the inversion
model_min_abs               = [permitt_hf_min; repmat([permitt_diff_min; freq_relax_min; alpha_cole_min], 4, 1)];
model_max_abs               = [permitt_hf_max; repmat([permitt_diff_max; freq_relax_max; alpha_cole_max], 4, 1)];

% convert relaxation frequency bounds into relaxation time constants, which are what the inversion uses to keep magnitudes similar
for ii = 1:4
    model_min_abs(3 + ((ii - 1) * 3)) ...
                            = 1 / (2 * pi * model_max_abs(3 + ((ii - 1) * 3)));
    model_max_abs(3 + ((ii - 1) * 3)) ...
                            = 1 / (2 * pi * model_min_abs(3 + ((ii - 1) * 3)));
end

% initialize a bunch of stuff (must be done prior to sub-functions)
load_all                    = struct;
[file_merge, file_save, path_curr, lr] ...
                            = deal('');
[ax, relax_check, period_relax_edit, permitt_diff_edit, freq_relax_edit, alpha_cole_edit, permitt_diff_box, freq_relax_box, alpha_cole_box, permitt_diff_slide, freq_relax_slide, alpha_cole_slide, p_model, p_freq_min, p_freq_max, activ_relax] ...
                            = deal(NaN(1, 4));
p_data                      = NaN(4, 2);
[curr_temp, curr_relax_var, freq, num_temp, permitt_real, permitt_imag, conduct, resist_phase, temp_mean, temp_std, temp_mean_alt, permitt_comp_model, permitt_hf, conduct_hf, conduct_dc, conduct_dc_std, freq_trim, permitt_real_trim, permitt_imag_trim, conduct_trim, weight_freq, model_guess, ...
 model_final, res_model, permitt_diff, freq_relax, alpha_cole, permitt_diff_std, freq_relax_std, alpha_cole_std, permitt_real_model, permitt_imag_model, conduct_model, phase_model, temp_min, temp_max, temp_vec, temp_inv_vec, temp_mean_inv, permitt_diff_cat, freq_relax_cat, alpha_cole_cat, ...
 conduct_hf_model, load_all, model_std, inv_done, permitt_hf_model, permitt_hf_model_std, model_min, model_max, relax_ord, model_final_old, num_inv, permitt_comp_model_full, permitt_real_model_full, conduct_model_full, permitt_imag_model_full, phase_model_full, curr_param, ...
 curr_range, param_tmp, model_tmp, res, res_def, chisq_diff, num_param, err, activ_dc, activ_dc_poly, curr_ind, ind_trim, permitt_comp_model_sep, permitt_real_model_sep, permitt_imag_model_sep, conduct_model_sep, length_freq, dia, curr_slide, tmp1, tmp2, period_relax, period_relax_std, ...
 period_relax_cat, period_relax_std_cat] ...
                            = deal(0);
[pdc, pdcf]                 = deal(NaN);
focus_check                 = true;
[permitt_diff_std_cat, freq_relax_std_cat, alpha_cole_std_cat, activ_relax_poly, period_relax_std_cat] ...
                            = deal(cell(4, 1));
[p_arr_fit, p_arr_err, p_arr_model] ...
                            = deal(NaN(4, 6));
p_model_sep                 = NaN(3, 6);
p_arr_leg                   = NaN(1, 6);
temp_str                    = {};

% intialize GUI parameters
freq_def                    = [1e-3 1e7]; % frequency range to display
ylim_def                    = [1 1e2; 1e-12 1e-2; 1e-3 1e3]; % ylims for each subplot
plots1                      = {'permitt_real{curr_temp}' 'conduct{curr_temp}' 'permitt_imag{curr_temp}' '(resist_phase{curr_temp} .* (180 / pi))'}; % data plot strings
plots2                      = {'permitt_real_model_full{curr_temp}' 'conduct_model_full{curr_temp}' 'permitt_imag_model_full{curr_temp}' '(phase_model_full{curr_temp} .* (180 / pi))'}; % model plot strings
plots3                      = {'permitt_real_model_sep{curr_temp}{jj}' 'conduct_model_sep{curr_temp}{jj}' 'permitt_imag_model_sep{curr_temp}{jj}'}; % model plot strings
ylabels                     = {'\epsilon''' '\sigma (S m^{-1})' '\epsilon''''' '\theta (\circ)'}; % ylabels for each supblot
relax_labels                = {'\Delta \epsilon' 'f_r' '\alpha'}; % labels for each slider
relax_labels_alt            = {'D-e-r' 'f-r' 'alpha'}; % labels for each slider
relax_simple                = {'permitt_diff' 'freq_relax' 'alpha_cole'}; % parameters within each relaxation
relax_str                   = {'%5.2f' '%2.2e' '%4.3f' '%1.3f'}; % sprintf precisions
subplot_start               = [0.05 0.47 0.35 0.35]; % start position of 1st subplot

% initialize odds and ends for the results plot
axf                         = NaN(1, 6); % axis handles
colors_fit                  = {'k' 'r' 'b' 'm'}; % easy colors for up to 4 different relaxations
ylabels_fit                 = {'{\Delta}{\epsilon}''' 'f_r (Hz)' '\alpha' '\epsilon''_{HF}' '\sigma_{DC}, \sigma_{HF} (S m^{-1})' 'Residual'}; % ylabels on fit plot
colors_sep                  = ['g' colors_fit 'c'];

% load existing relaxation frequency background data
auty_cole                   = load([dir_arrhenius 'Auty_Cole_1952']);
temp_inv_AC                 = 1e3 ./ (auty_cole.temp + 273.15); % inverse temperature
freq_relax_AC               = auty_cole.freq_relax;
kawada                      = load([dir_arrhenius 'Kawada_1978']);
temp_inv_K                  = 1e3 ./ (kawada.temp + 273.15);
freq_relax_K                = kawada.freq_relax;
sat                         = load([dir_arrhenius 'NaCl_Sat_Tau']);
temp_inv_sat                = 1e3 ./ sat.T;
freq_relax_sat              = 1 ./ ((2 * pi) .* sat.tau);

%% draw the GUI
set(0, 'DefaultFigureWindowStyle', 'docked')
inv_gui                     = figure('toolbar', 'figure', 'name', 'INV_COLE_GUI', 'position', [1 1 1600 1200], 'keypressfcn', @keypress);
for ii = 1:4
    subplot_curr            = subplot_start;
    if any(ii == [2 4])
        subplot_curr        = subplot_curr + [0.41 0 0 0];
    end
    if (ii > 2)
        subplot_curr        = subplot_curr + [0 -0.38 0 0];
    end
    ax(ii)                  = subplot('position', subplot_curr);
    hold on
    set(gca, 'fontsize', 20, 'xscale', 'log', 'color', [0.9 0.9 0.9])
    if (ii < 4) % all plots except phase should be y log
        set(gca, 'yscale', 'log')
        axis([freq_def ylim_def(ii, :)])
    else
        axis([freq_def -90 0])
    end
    if (ii < 4)
        p_freq_min(ii)      = loglog(repmat(freq_min, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
    else
        p_freq_min(ii)      = semilogx(repmat(freq_min, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
    end
    if (ii < 4)
        p_freq_max(ii)      = loglog(repmat(freq_max, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
    else
        p_freq_max(ii)      = semilogx(repmat(freq_max, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
    end
    if (ii > 2)
        xlabel('Frequency (Hz)')
    else
        set(gca, 'xticklabel', {})
    end
    ylabel(ylabels{ii})
    grid on
    box on
end
linkaxes(ax, 'x')
annotation('textbox', [0.01 0.95 0.03 0.04], 'string', 'Load', 'fontsize', 16, 'color', 'k', 'edgecolor', 'none')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'data', 'units', 'normalized', 'position', [0.06 0.95 0.055 0.04], 'callback', @load_merge, 'fontsize', 16, 'foregroundcolor', 'k')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'config.', 'units', 'normalized', 'position', [0.13 0.95 0.06 0.04], 'callback', @load_cfg, 'fontsize', 16, 'foregroundcolor', 'k')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'Clear plot', 'units', 'normalized', 'position', [0.01 0.84 0.08 0.04], 'callback', @nuke_plot, 'fontsize', 16, 'foregroundcolor', 'r')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'Erase', 'units', 'normalized', 'position', [0.74 0.83 0.07 0.04], 'callback', @nuke_inv, 'fontsize', 16, 'foregroundcolor', 'r')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'Invert', 'units', 'normalized', 'position', [0.74 0.89 0.07 0.04], 'callback', @do_inv, 'fontsize', 16, 'foregroundcolor', 'b')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'Next', 'units', 'normalized', 'position', [0.59 0.95 0.07 0.04], 'callback', @do_next, 'fontsize', 16, 'foregroundcolor', 'b')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'Save data', 'units', 'normalized', 'position', [0.67 0.95 0.07 0.04], 'callback', @save_inv, 'fontsize', 16, 'foregroundcolor', 'g')
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'config.', 'units', 'normalized', 'position', [0.74 0.95 0.07 0.04], 'callback', @save_cfg, 'fontsize', 16, 'foregroundcolor', 'g')
weight_check                = uicontrol(inv_gui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.10 0.84 0.08 0.04], 'string', 'weight', 'callback', @check_weight, 'fontsize', 16, 'value', double(do_weight));
sep_check                   = uicontrol(inv_gui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.15 0.84 0.08 0.04], 'string', 'disp. sep.', 'callback', @check_sep, 'fontsize', 16, 'value', double(do_sep));
file_box                    = annotation('textbox', [0.01 0.89 0.20 0.04], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'interpreter', 'none', 'linewidth', 1);
annotation('textbox', [0.22 0.94 0.05 0.05], 'string', 'f_{min}', 'fontsize', 16, 'edgecolor', 'none', 'color', 'b')
freq_min_edit               = annotation('textbox', [0.26 0.95 0.05 0.03], 'string', sprintf(relax_str{2}, freq_min), 'fontsize', 16, 'color', 'b', 'edgecolor', 'none', 'linewidth', 1);
freq_min_slide              = uicontrol(inv_gui, 'style', 'slider', 'units', 'normalized', 'position', [0.22 0.91 0.17 0.02], 'callback', @slide_freq_min, 'min', -3, 'max', 4, 'value', log10(freq_min), 'sliderstep', [0.01 0.1]);
annotation('textbox', [0.22 0.85 0.05 0.05], 'string', 'f_{max}', 'fontsize', 16, 'edgecolor', 'none', 'color', 'b')
freq_max_edit               = annotation('textbox', [0.26 0.87 0.05 0.03], 'string', sprintf(relax_str{2}, freq_max), 'fontsize', 16, 'color', 'b', 'edgecolor', 'none', 'linewidth', 1);
freq_max_slide              = uicontrol(inv_gui, 'style', 'slider', 'units', 'normalized', 'position', [0.22 0.83 0.17 0.02], 'callback', @slide_freq_max, 'min', 3, 'max', 6, 'value', log10(freq_max), 'sliderstep', [0.01 0.1]);
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'eHF', 'units', 'normalized', 'position', [0.40 0.94 0.04 0.04], 'callback', @lim_permitt_hf, 'fontsize', 16, 'foregroundcolor', 'b')
permitt_hf_edit             = annotation('textbox', [0.445 0.95 0.05 0.03], 'string', sprintf(relax_str{1}, permitt_hf_guess), 'fontsize', 16, 'color', 'b', 'edgecolor', 'none', 'linewidth', 1);
permitt_hf_box              = annotation('textbox', [0.50 0.95 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'color', 'r', 'linewidth', 1);
permitt_hf_slide            = uicontrol(inv_gui, 'style', 'slider', 'units', 'normalized', 'position', [0.40 0.91 0.17 0.02], 'callback', @slide_permitt_hf, 'min', permitt_hf_min, 'max', permitt_hf_max, 'value', permitt_hf_guess);
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'DC', 'units', 'normalized', 'position', [0.40 0.85 0.04 0.04], 'callback', @lim_conduct_dc, 'fontsize', 16, 'foregroundcolor', 'b')
conduct_dc_edit             = annotation('textbox', [0.445 0.87 0.05 0.03], 'string', sprintf(relax_str{2}, conduct_dc_guess), 'fontsize', 16, 'color', 'b', 'edgecolor', 'none', 'linewidth', 1);
conduct_dc_box              = annotation('textbox', [0.50 0.87 0.05 0.03], 'string', '', 'color', 'k', 'fontsize', 16, 'backgroundcolor', 'w', 'edgecolor', 'k', 'color', 'r', 'linewidth', 1);
conduct_dc_check            = uicontrol(inv_gui, 'style', 'checkbox', 'units', 'normalized', 'position', [0.56 0.86 0.02 0.03], 'callback', @check_dc, 'fontsize', 16, 'value', double(do_dc));
conduct_dc_slide            = uicontrol(inv_gui, 'style', 'slider', 'units', 'normalized', 'position', [0.40 0.83 0.17 0.02], 'callback', @slide_conduct_dc, 'min', log10(conduct_dc_min), 'max', log10(conduct_dc_max), 'value', log10(conduct_dc_guess));
temp_box                    = uicontrol(inv_gui, 'style', 'popupmenu', 'string', 'N/A', 'units', 'normalized', 'position', [0.59 0.89 0.07 0.05], 'callback', @plot_data, 'fontsize', 16);
fm_edit                     = annotation('textbox', [0.59 0.83 0.06 0.04], 'string', '', 'color', 'b', 'backgroundcolor', 'w', 'edgecolor', 'k', 'fontsize', 16, 'linewidth', 1);
inv_edit                    = annotation('textbox', [0.66 0.83 0.06 0.04], 'string', '', 'color', 'r', 'backgroundcolor', 'w', 'edgecolor', 'k', 'fontsize', 16, 'linewidth', 1);
inv_done_edit               = annotation('textbox', [0.67 0.89 0.04 0.04], 'string', '0 /', 'color', 'r', 'fontsize', 16, 'edgecolor', 'none', 'linewidth', 1);
inv_tot                     = annotation('textbox', [0.70 0.89 0.03 0.04], 'string', '', 'color', 'r', 'fontsize', 16, 'edgecolor', 'none', 'linewidth', 1);
relax_group                 = uibuttongroup('position', [0.82 0.91 0.17 0.08], 'selectionchangefcn', @relax_radio);
uicontrol(inv_gui, 'style', 'text', 'parent', relax_group, 'units', 'normalized', 'position', [0 0.6 0.9 0.3], 'string', '# relaxations', 'fontsize', 16)
for ii = 1:4
    relax_check(ii)         = uicontrol(inv_gui, 'style', 'radio', 'string', num2str(ii), 'units', 'normalized', 'position', [(0.05 + ((ii - 1) * 0.25)) 0 0.20 0.5], 'parent', relax_group, 'fontsize', 16, 'handlevisibility', 'off');
end
set(relax_group, 'selectedobject', relax_check(num_relax))
for ii = 1:4 % set up sliders
    for jj = 1:3
        eval(['uicontrol(inv_gui, ''style'', ''pushbutton'', ''units'', ''normalized'', ''position'', [0.82 ' num2str(0.84 - ((ii - 1) * 0.22) - ((jj - 1) * 0.06)) ' 0.05 0.04], ''callback'', @lim_' relax_simple{jj} '_' num2str(ii) ', ''fontsize'', 16, ''string'', [''' ...
              relax_labels_alt{jj} ''' ''(' num2str(ii) ')''], ''foregroundcolor'', ''' colors_fit{ii} ''')'])
        eval([relax_simple{jj} '_edit(' num2str(ii) ') = annotation(''textbox'', [0.88 ' num2str(0.84 - ((ii - 1) * 0.22) - ((jj - 1) * 0.06)) ' 0.05 0.04], ''string'', sprintf(''' relax_str{jj} ''', relax_guess(' num2str(ii) ',' num2str(jj) ...
              ')), ''fontsize'', 16, ''color'', ''b'', ''edgecolor'', ''none'', ''linewidth'', 1);'])
        eval([relax_simple{jj} '_box(' num2str(ii) ') = annotation(''textbox'', [0.93 ' num2str(0.85 - ((ii - 1) * 0.22) - ((jj - 1) * 0.06)) ' 0.06 0.03], ''string'', '''', ''fontsize'', 16, ''backgroundcolor'', ''w'', ''edgecolor'', ''k'', ''color'', ''r'', ''linewidth'', 1);'])
        eval([relax_simple{jj} '_slide(' num2str(ii) ') = uicontrol(inv_gui, ''style'', ''slider'', ''units'', ''normalized'', ''position'', [0.82 ' ...
            num2str(0.82 - ((ii - 1) * 0.22) - ((jj - 1) * 0.06)) ' 0.17 0.02], ''callback'', @slide_' relax_simple{jj} '_' num2str(ii) ', ''fontsize'', 16, ''min'', ' relax_simple{jj} '_min, ''max'', ' relax_simple{jj} '_max, ''value'', ' num2str(relax_guess(ii, jj)) ');'])
        if (jj == 2)
            set(eval(['' relax_simple{jj} '_edit(' num2str(ii) ')' '']), 'position', (eval(['' 'get(' relax_simple{jj} '_edit(' num2str(ii) '), ''position'')' '']) + [-0.01 0 0 0]))
            set(eval(['' relax_simple{jj} '_slide(' num2str(ii) ')' '']), 'min', log10(eval(['' relax_simple{jj} '_min' ''])), 'max', log10(eval(['' relax_simple{jj} '_max' ''])), 'value', log10(relax_guess(ii, jj)))
            period_relax_edit(ii) ...
                            = annotation('textbox', [0.88 (0.82 - ((ii - 1) * 0.22) - ((jj - 1) * 0.06)) 0.05 0.04], 'string', sprintf(relax_str{jj}, (1 / (2 * pi * relax_guess(ii, jj)))), 'fontsize', 14, 'color', 'b', 'edgecolor', 'none', 'linewidth', 1);
        end
        set(eval([relax_simple{jj} '_slide(' num2str(ii) ')']), 'sliderstep', ([0.1 0.1] .* get(eval([relax_simple{jj} '_slide(' num2str(ii) ')']), 'sliderstep')))
    end
end
uicontrol(inv_gui, 'style', 'pushbutton', 'string', 'test', 'units', 'normalized', 'position', [0.945 0.95 0.04 0.04], 'callback', @misctest, 'fontsize', 16, 'foregroundcolor', 'b')

pause(1)

%% arrhenius inversion results plot
arr_fig                     = figure('toolbar', 'figure', 'name', 'INV_RESULTS', 'position', [1 1 1600 1200], 'keypressfcn', @keypress2);
for ii = 1:6
    axf(ii)                 = subplot(2, 3, ii);
    hold on
    if (ii == 2)
        semilogy(temp_inv_AC, freq_relax_AC, 'c', 'linewidth', 2);
        semilogy(temp_inv_K, freq_relax_K, 'g', 'linewidth', 2);
        semilogy(temp_inv_sat, freq_relax_sat, 'color', [0.7 0.7 0.7], 'linewidth', 2);
    end
    if any(ii == [2 5])
        set(gca, 'yscale', 'log')
    end
    set(gca, 'fontsize', 20, 'color', [0.9 0.9 0.9])
    ylabel(ylabels_fit{ii})
    grid on
    box on
    xlim([min([temp_inv_AC(:); temp_inv_K(:); temp_inv_sat(:)]) max([temp_inv_AC(:); temp_inv_K(:); temp_inv_sat(:)])])
end
linkaxes(axf, 'x')

%%
figure(inv_gui)

%% load data

    function load_merge(source, eventdata) %#ok<*INUSD>
        
        figure(inv_gui)
        
        if ~focus_check
            focus_check     = true;
        end
        
        curr_temp           = 1;
        num_inv             = 0;
        inv_done            = false;
        nuke_inv
        
        % remove previously loaded data and forward model
        set(file_box, 'string', '')
        delete(p_data(ishandle(p_data(:))))
        delete(p_model(ishandle(p_model(:))))
        delete(p_model_sep(ishandle(p_model_sep(:))))
        
        % dialog box to get merged data file to load
        if ~isempty(path_curr)
            [file_merge, path_curr] ...
                            = uigetfile('*.mat', 'Load merged data:', path_curr);
        else
            [file_merge, path_curr] ...
                            = uigetfile('*.mat', 'Load merged data:');
        end
        if ~ischar(file_merge)
            [file_merge, path_curr] ...
                            = deal('');
        end
        
        if isempty(file_merge)
            return
        end
        
        % load file and display in filename box
        load_all            = load([path_curr file_merge]);
        set(file_box, 'string', file_merge(1:(end - 4)))
        
        % distribute loaded data to better variable names
        [freq, permitt_real, permitt_imag, conduct, resist_phase, temp_mean, temp_std] ...
                            = deal(load_all.freq, load_all.permitt_real, load_all.permitt_imag, load_all.conduct, load_all.resist_phase, load_all.temp_mean, load_all.temp_std); %#ok<*SETNU>
        
        num_temp            = length(temp_mean); % number of temperatures is determined from temp_mean
        
        temp_str            = cell(num_temp, 1);
        % put all temperatures in the pull-down menu
        for ii = 1:num_temp
            temp_str{ii}    = num2str(round(temp_mean(ii)));
        end
        set(temp_box, 'string', temp_str)
        
        % initialize a whole bunch of variables now that we have num_temp (note they start life as NaNs not zeros)
        [permitt_hf, conduct_hf, conduct_dc, conduct_hf_model, permitt_hf_model, res_def, chisq_diff, num_param, res_model, length_freq] ...
                            = deal(NaN(num_temp, 1));
        for ii = 1:num_temp
            length_freq(ii) = length(freq{ii});
        end
        [permitt_hf_model_std, conduct_dc_std] ...
                            = deal(NaN(num_temp, 2));
        [freq_trim, permitt_real_trim, permitt_imag_trim, conduct_trim, weight_freq, model_guess, permitt_real_model, permitt_imag_model, conduct_model, phase_model, model_final, permitt_diff, freq_relax, alpha_cole, permitt_real_model_full, permitt_imag_model_full, conduct_model_full, ...
         phase_model_full, model_std, permitt_diff_std, freq_relax_std, alpha_cole_std, model_min, model_max, permitt_real_model_sep, permitt_imag_model_sep, conduct_model_sep, period_relax, period_relax_std] ...
                            = deal(cell(num_temp, 1));
        [permitt_diff_cat, freq_relax_cat, alpha_cole_cat, period_relax_cat] ...
                            = deal(NaN(4, num_temp));
        
        for ii = 1:4
            [permitt_diff_std_cat{ii}, freq_relax_std_cat{ii}, alpha_cole_std_cat{ii}, period_relax_std_cat{ii}] ...
                            = deal(NaN(num_temp, 2));
        end
        num_relax           = 2 .* ones(1, num_temp); % default to two relaxations
        [do_dc, inv_done]   = deal(false(1, num_temp)); % default no DC conductivity
        
        temp_mean_alt       = (1 / (temp_mean(1) + 273.15)) - (1 ./ (temp_mean + 273.15)); % Arrhenius inverse temperature of temp_mean
        
        % odds and ends for fit plot
        temp_min            = min(temp_mean) - mod(min(temp_mean), 10) - 10;
        temp_max            = max(temp_mean) - mod(max(temp_mean), 10) + 10;
        temp_vec            = temp_min:10:temp_max;
        temp_inv_vec        = 1e3 ./ (temp_vec + 273.15);
        temp_mean_inv       = 1e3 ./ (temp_mean + 273.15);
        
        set(temp_box, 'value', 1)
        set(inv_tot, 'string', num2str(num_temp))
        
        data_loaded         = true;        
        plot_data % plot the first (should be lowest) temperature 
    end

%% load configuration file

    function load_cfg(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        if ~isempty(path_curr)
            [file_merge, path_curr] ...
                            = uigetfile('*.mat', 'Load configuration:', path_curr);
        else
            [file_merge, path_curr] ...
                            = uigetfile('*.mat', 'Load configuration:');
        end
        if ~file_merge
            [file_merge, path_curr] ...
                            = deal('');
        end
        if ~isempty(file_merge)
            load_all        = load([path_curr file_merge]);
            try
                [freq_min, freq_max, freq_trim, do_weight, weight_freq, num_relax, do_dc, relax_guess, model_guess, permitt_hf_guess, conduct_dc_guess, model_final, permitt_hf, permitt_hf_model, permitt_hf_model_std, conduct_hf, conduct_hf_model, conduct_dc, conduct_dc_std, ...
                 permitt_real_model, permitt_imag_model, conduct_model, phase_model, res_model, freq_relax, freq_relax_cat, freq_relax_std, freq_relax_std_cat, alpha_cole, alpha_cole_cat, alpha_cole_std, alpha_cole_std_cat, permitt_diff, permitt_diff_cat, permitt_diff_std, ...
                 permitt_diff_std_cat, permitt_real_trim, permitt_imag_trim, activ_relax, activ_relax_poly, activ_dc, activ_dc_poly, period_relax, period_relax_std, period_relax_cat, period_relax_std_cat] ...
                            = deal(load_all.freq_min, load_all.freq_max, load_all.freq_trim, load_all.do_weight, load_all.weight_freq, load_all.num_relax, load_all.do_dc, load_all.relax_guess, load_all.model_guess, load_all.permitt_hf_guess, load_all.conduct_dc_guess, load_all.model_final, ...
                                   load_all.permitt_hf, load_all.permitt_hf_model, load_all.permitt_hf_model_std, load_all.conduct_hf, load_all.conduct_hf_model, load_all.conduct_dc, load_all.conduct_dc_std, load_all.permitt_real_model, load_all.permitt_imag_model, load_all.conduct_model, ...
                                   load_all.phase_model, load_all.res_model, load_all.freq_relax, load_all.freq_relax_cat, load_all.freq_relax_std, load_all.freq_relax_std_cat, load_all.alpha_cole, load_all.alpha_cole_cat, load_all.alpha_cole_std, load_all.alpha_cole_std_cat, ...
                                   load_all.permitt_diff, load_all.permitt_diff_cat, load_all.permitt_diff_std, load_all.permitt_diff_std_cat, load_all.permitt_real_trim, load_all.permitt_imag_trim, load_all.activ_relax, load_all.activ_relax_poly, load_all.activ_dc, load_all.activ_dc_poly, ...
                                   load_all.period_relax, load_all.period_relax_std, load_all.period_relax_cat, load_all.period_relax_std_cat);
            catch %#ok<CTCH>
                [freq_min, freq_max, freq_trim, do_weight, weight_freq, num_relax, do_dc, relax_guess, model_guess, permitt_hf_guess, conduct_dc_guess, model_final, permitt_hf, permitt_hf_model, permitt_hf_model_std, conduct_hf, conduct_hf_model, conduct_dc, conduct_dc_std, ...
                 permitt_real_model, permitt_imag_model, conduct_model, phase_model, res_model, freq_relax, freq_relax_cat, freq_relax_std, freq_relax_std_cat, alpha_cole, alpha_cole_cat, alpha_cole_std, alpha_cole_std_cat, permitt_diff, permitt_diff_cat, permitt_diff_std, ...
                 permitt_diff_std_cat, permitt_real_trim, permitt_imag_trim, activ_relax, activ_relax_poly, activ_dc, activ_dc_poly] ...
                            = deal(load_all.freq_min, load_all.freq_max, load_all.freq_trim, load_all.do_weight, load_all.weight_freq, load_all.num_relax, load_all.do_dc, load_all.relax_guess, load_all.model_guess, load_all.permitt_hf_guess, load_all.conduct_dc_guess, load_all.model_final, ...
                                   load_all.permitt_hf, load_all.permitt_hf_model, load_all.permitt_hf_model_std, load_all.conduct_hf, load_all.conduct_hf_model, load_all.conduct_dc, load_all.conduct_dc_std, load_all.permitt_real_model, load_all.permitt_imag_model, load_all.conduct_model, ...
                                   load_all.phase_model, load_all.res_model, load_all.freq_relax, load_all.freq_relax_cat, load_all.freq_relax_std, load_all.freq_relax_std_cat, load_all.alpha_cole, load_all.alpha_cole_cat, load_all.alpha_cole_std, load_all.alpha_cole_std_cat, ...
                                   load_all.permitt_diff, load_all.permitt_diff_cat, load_all.permitt_diff_std, load_all.permitt_diff_std_cat, load_all.permitt_real_trim, load_all.permitt_imag_trim, load_all.activ_relax, load_all.activ_relax_poly, load_all.activ_dc, load_all.activ_dc_poly);
                [period_relax, period_relax_std] ...
                            = deal(cell(num_temp, 1));
                for ii = 1:num_temp
                    [period_relax{ii}, period_relax_std{ii}] ...
                            = deal((1 ./ ((2 * pi) .* freq_relax{ii})), (1 ./ ((2 * pi) .* freq_relax_std{ii})));
                end
                period_relax_cat ...
                            = 1 ./ ((2 * pi) .* freq_relax_cat);
                for ii = 1:4
                    period_relax_std_cat{ii} ...
                            = 1 ./ ((2 * pi) .* freq_relax_std_cat{ii});
                end
            end
            set(weight_check, 'value', double(do_weight))
            set(freq_min_edit, 'string', sprintf(relax_str{2}, freq_min))
            set(freq_min_slide, 'value', log10(freq_min))
            set(freq_max_edit, 'string', sprintf(relax_str{2}, freq_max))
            set(freq_max_slide, 'value', log10(freq_max))
            set(permitt_hf_slide, 'value', permitt_hf_guess)
            set(permitt_hf_edit, 'string', sprintf(relax_str{1}, permitt_hf_guess))
            set(conduct_dc_slide, 'value', log10(conduct_dc_guess))
            set(conduct_dc_edit, 'string', sprintf(relax_str{2}, conduct_dc_guess))
            set(conduct_dc_check, 'value', double(do_dc(1)))
            set(relax_group, 'selectedobject', relax_check(num_relax(1)))
            for ii = 1:3
                for jj = 1:4
                    if (ii == 2)
                        if (get(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'min') > log10(relax_guess(jj, ii)))
                            set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'min', log10(relax_guess(jj, ii)))
                        elseif (get(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'max') < log10(relax_guess(jj, ii)))
                            set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'max', log10(relax_guess(jj, ii)))
                        end                        
                        set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'value', log10(relax_guess(jj, ii)))
                    else
                        if (get(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'min') > relax_guess(jj, ii))
                            set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'min', relax_guess(jj, ii))
                        elseif (get(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'max') < relax_guess(jj, ii))
                            set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'max', relax_guess(jj, ii))
                        end
                        set(eval([relax_simple{ii} '_slide(' num2str(jj) ')']), 'value', relax_guess(jj, ii))
                    end
                    set(eval([relax_simple{ii} '_edit(' num2str(jj) ')']), 'string', sprintf(relax_str{ii}, relax_guess(jj, ii)))
                end
            end
            for ii = 1:4
                set(period_relax_edit(ii), 'string', num2str(1 / (2 * pi * relax_guess(ii, 2))))
            end
            temp_str        = cell(num_temp, 1);
            inv_done        = false(1, num_temp);
            for ii = 1:num_temp
                temp_str{ii} ...
                            = num2str(round(temp_mean(ii)));
                if ~isempty(model_final{ii})
                    temp_str{ii} ...
                            = [temp_str{ii} 'X'];
                    inv_done(ii) ...
                            = true;
                end
            end
            set(temp_box, 'string', temp_str)
        end
        curr_temp           = 1;
        num_inv             = length(find(~isnan(res_model)));
        set(inv_done_edit, 'string', [num2str(num_inv) ' /'])
        inv_report
        plot_fits
        do_fm
        
    end

%% plot a single temperature's data

    function plot_data(source, eventdata)
        
        figure(inv_gui)
        
        if ~focus_check
            focus_check     = true;
        end
        
        curr_temp           = get(temp_box, 'value'); % current temperature is determined from pull-down menu (very important variable)
        
        if any(ishandle(p_data(:))) % get rid of old plotted data
            delete(p_data(ishandle(p_data(:))))
        end
        
        % plot the data for current temperature
        for ii = 1:4 %#ok<*FXUP>
            axes(ax(ii)) %#ok<*LAXES>
            if (ii < 4)
                p_data(ii, 1) ...
                            = loglog(freq{curr_temp}, eval(plots1{ii}), 'r', 'linewidth', 2);
                p_data(ii, 2) ...
                            = loglog(freq{curr_temp}, eval(plots1{ii}), 'ko', 'markerfacecolor', 'r', 'markersize', 8, 'linewidth', 1);
                set(ax(ii), 'ylim', [min(eval(plots1{ii})) max(eval(plots1{ii}))])
            else
                p_data(ii, 1) ...
                            = semilogx(freq{curr_temp}, eval(plots1{ii}), 'r', 'linewidth', 2);
                p_data(ii, 2) ...
                            = semilogx(freq{curr_temp}, eval(plots1{ii}), 'ko', 'markerfacecolor', 'r', 'markersize', 8, 'linewidth', 1);
            end
        end
        xlim([min(freq{curr_temp}) max(freq{curr_temp})]) % restrict x axis to current data's limits
        
        set(conduct_dc_check, 'value', do_dc(curr_temp)) % to DC or not to DC
        
        % figure out if we've already inverted for this temperature, and adjust based on outcome
        if ~isempty(model_final{curr_temp})
            set(freq_min_slide, 'value', log10(min(freq_trim{curr_temp})))
            set(freq_max_slide, 'value', log10(max(freq_trim{curr_temp})))
            inv_report
        else
            set([inv_edit, permitt_hf_box conduct_dc_box], 'string', '')
            set([permitt_diff_box freq_relax_box alpha_cole_box], 'string', '')
        end
        slide_freq_min
        slide_freq_max
        set(relax_group, 'selectedobject', relax_check(num_relax(curr_temp)))
        do_fm
    end

%% forward model

    function do_fm(source, eventdata)
        
        figure(inv_gui)
        
        if any(ishandle(p_model(:))) % get rid of old plotted forward model
            delete(p_model(ishandle(p_model(:))))
        end
        if any(ishandle(p_model_sep(:)))
            delete(p_model_sep(ishandle(p_model_sep(:))))
        end
        
        % initalize model guess
        if do_dc(curr_temp) % throw DC conductivity at end of model if including
            model_guess{curr_temp} ...
                            = NaN((2 + (3 * num_relax(curr_temp))), 1);
            model_guess{curr_temp}(end) ...
                            = -log10(conduct_dc_guess); % convert DC conductivity to log of DC resitivity (small and not negative), makes inversion work better
        else
            model_guess{curr_temp} ...
                            = NaN((1 + (3 * num_relax(curr_temp))), 1);
        end
        
        % assign HF permittivity to first model bin
        model_guess{curr_temp}(1) ...
                            = permitt_hf_guess;
        
        % load in each relaxation
        for ii = 1:num_relax(curr_temp) 
            [model_guess{curr_temp}(2 + ((ii - 1) * 3)), model_guess{curr_temp}(3 + ((ii - 1) * 3)), model_guess{curr_temp}(4 + ((ii - 1) * 3))] ...
                            = deal(relax_guess(ii, 1), (1 / (2 * pi * relax_guess(ii, 2))), relax_guess(ii, 3));
        end
        
        % trim frequencies as per freq_min and freq_max
        ind_trim            = find((freq{curr_temp} >= freq_min) & (freq{curr_temp} <= freq_max)); 
        [freq_trim{curr_temp}, permitt_real_trim{curr_temp}, permitt_imag_trim{curr_temp}, conduct_trim{curr_temp}] ...
                            = deal(freq{curr_temp}(ind_trim), permitt_real{curr_temp}(ind_trim), permitt_imag{curr_temp}(ind_trim), conduct{curr_temp}(ind_trim));
        ind_remove          = zeros(length(freq_repeat), 1);
        for ii = 1:length(freq_repeat)
            if any(freq_trim{curr_temp} == freq_repeat(ii));
                ind_remove(ii) ...
                            = find((freq_trim{curr_temp} == freq_repeat(ii)), 1, 'last'); % trim repeated frequencies
            end
        end
        ind_remove          = ind_remove(logical(ind_remove));
        ind_remove          = unique([ind_remove; find(isnan(permitt_real_trim{curr_temp}) | isnan(permitt_imag_trim{curr_temp}))]);
        if ~isempty(ind_remove)
            [freq_trim{curr_temp}, permitt_real_trim{curr_temp}, permitt_imag_trim{curr_temp}, conduct_trim{curr_temp}] ...
                            = deal(freq_trim{curr_temp}(setdiff(1:length(freq_trim{curr_temp}), ind_remove)), permitt_real_trim{curr_temp}(setdiff(1:length(freq_trim{curr_temp}), ind_remove)), permitt_imag_trim{curr_temp}(setdiff(1:length(freq_trim{curr_temp}), ind_remove)), ...
                                   conduct_trim{curr_temp}(setdiff(1:length(freq_trim{curr_temp}), ind_remove)));
        end
        
        % run forward model and break out data into plotted parts
        permitt_comp_model  = fm_cole(model_guess{curr_temp}, freq_trim{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum);
        [permitt_real_model{curr_temp}, permitt_imag_model{curr_temp}] ...
                            = deal(permitt_comp_model(1:(length(freq_trim{curr_temp}))), permitt_comp_model((1 + length(freq_trim{curr_temp})):end)); % break out modeled complex permittivity into real and imaginary parts
        conduct_model{curr_temp} ...
                            = (2 * pi * permitt_vacuum) .* freq_trim{curr_temp} .* permitt_imag_model{curr_temp}; % modeled conductivity derived from modeled permittivity
        phase_model{curr_temp} ...
                            = -atan(permitt_real_model{curr_temp} ./ permitt_imag_model{curr_temp});
        
        % redo forward to show data outside of trimmed frequency range
        permitt_comp_model_full ...
                            = fm_cole(model_guess{curr_temp}, freq{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum);
        [permitt_real_model_full{curr_temp}, permitt_imag_model_full{curr_temp}] ...
                            = deal(permitt_comp_model_full(1:length_freq(curr_temp)), permitt_comp_model_full((1 + length_freq(curr_temp)):end)); % break out modeled complex permittivity into real and imaginary parts
        conduct_model_full{curr_temp} ...
                            = (2 * pi * permitt_vacuum) .* freq{curr_temp} .* permitt_imag_model_full{curr_temp}; % modeled complex conductivity derived from modeled permittivity
        phase_model_full{curr_temp} ...
                            = -atan(permitt_real_model_full{curr_temp} ./ permitt_imag_model_full{curr_temp});
        conduct_hf_model(curr_temp) ...
                            = conduct_model{curr_temp}(1); % HF conductivity at max frequency (needs to be ~1 MHz)
        
        % separate each component of complex permittivity
        permitt_comp_model_sep ...
                            = fm_cole_sep(model_guess{curr_temp}, freq{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum);
        [permitt_real_model_sep{curr_temp}, permitt_imag_model_sep{curr_temp}, conduct_model_sep{curr_temp}] ...
                            = deal(cell(1, length(permitt_comp_model_sep)));
        permitt_real_model_sep{curr_temp}{1} ...
                            = permitt_comp_model_sep{1}(1:(length_freq(curr_temp)));
        permitt_imag_model_sep{curr_temp}{1} ...
                            = NaN(length_freq(curr_temp), 1);
        for ii = 1:num_relax(curr_temp)
            [permitt_real_model_sep{curr_temp}{ii + 1}, permitt_imag_model_sep{curr_temp}{ii + 1}] ...
                            = deal(permitt_comp_model_sep{ii + 1}(1:length_freq(curr_temp)), permitt_comp_model_sep{ii + 1}((1 + length_freq(curr_temp)):end));
            permitt_real_model_sep{curr_temp}{ii + 1} ...
                            = permitt_real_model_sep{curr_temp}{ii + 1} + permitt_real_model_sep{curr_temp}{ii};
        end
        if do_dc(curr_temp)
            permitt_real_model_sep{curr_temp}{end} ...
                            = NaN(length_freq(curr_temp), 1);
            permitt_imag_model_sep{curr_temp}{end} ...
                            = permitt_comp_model_sep{end}((1 + length_freq(curr_temp)):end);
        end
        for ii = 1:length(permitt_comp_model_sep)
            conduct_model_sep{curr_temp}{ii} ...
                            = (2 * pi * permitt_vacuum) .* freq{curr_temp} .* permitt_imag_model_sep{curr_temp}{ii};
        end
        
        % plot forward model
        if do_sep
            disp_sep
        end
        for ii = 1:4
            axes(ax(ii))
            if (ii < 4)
                p_model(ii) = loglog(freq{curr_temp}, eval(plots2{ii}), 'b--', 'linewidth', 2);
            else
                p_model(ii) = semilogx(freq{curr_temp}, eval(plots2{ii}), 'b--', 'linewidth', 2);
            end
        end
        
        % do weights if necessary so that normalized residual can be computed
        if do_weight
            weight
        else
            weight_freq{curr_temp} ...
                            = ones((2 * length(freq_trim{curr_temp})), 1);
        end
        
        set(fm_edit, 'string', sprintf(relax_str{2}, sum(res_cole(model_guess{curr_temp}, freq_trim{curr_temp}, permitt_real_trim{curr_temp}, permitt_imag_trim{curr_temp}, weight_freq{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum) .^ 2)))
        
    end

%% inverse fitting
    
    function do_inv(source, eventdata)
        
        figure(inv_gui)
        
        if ~focus_check
            focus_check     = true;
        end
        
        delete(p_model(ishandle(p_model(:))))
        delete(p_model_sep(ishandle(p_model_sep(:))))
        
        % get data's HF permittivity and conductivity from first frequency (should be highest)
        [permitt_hf(curr_temp), conduct_hf(curr_temp)] ...
                            = deal(permitt_real_trim{curr_temp}(1), conduct_trim{curr_temp}(1));
        
        % empty inverse fit boxes
        set(permitt_hf_box, 'string', '')
        if ~isempty(get(conduct_dc_box, 'string'))
            set(conduct_dc_box, 'string', '')
        end
        set([permitt_diff_box freq_relax_box alpha_cole_box], 'string', '')
        
        % assign model min/max based on current number of relaxations
        model_min{curr_temp}= model_min_abs(1:(1 + (3 * num_relax(curr_temp))));
        model_max{curr_temp}= model_max_abs(1:(1 + (3 * num_relax(curr_temp))));
        if do_dc(curr_temp)
            model_min{curr_temp} ...
                            = [model_min{curr_temp}; -log10(conduct_dc_max)];
            model_max{curr_temp} ...
                            = [model_max{curr_temp}; -log10(conduct_dc_min)];
        end
        
        % weight (or not) the frequencies, typically this should be done
        if do_weight
            weight
        else
            weight_freq{curr_temp} ...
                            = ones((2 * length(freq_trim{curr_temp})), 1);
        end
        
        % run the inversion
        [model_final{curr_temp}, res_model(curr_temp)] ...
                            = lsqnonlin(@res_cole, model_guess{curr_temp}, model_min{curr_temp}, model_max{curr_temp}, options, freq_trim{curr_temp}, permitt_real_trim{curr_temp}, permitt_imag_trim{curr_temp}, weight_freq{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum);
        
        % add 1 to completed inversions displayed if this one hasn't already been done
        if ~inv_done(curr_temp)
            num_inv         = num_inv + 1;
            set(inv_done_edit, 'string', [num2str(num_inv) ' /'])
            inv_done(curr_temp) ...
                            = true;
            temp_str{curr_temp} ...
                            = [temp_str{curr_temp} 'X'];
            set(temp_box, 'string', temp_str)
        end
        
        % assign inverse model results
        permitt_hf_model(curr_temp) ...
                            = model_final{curr_temp}(1);
        if (permitt_hf_model(curr_temp) < permitt_hf_min)
            permitt_hf_model(curr_temp) ...
                            = permitt_hf_min;
        elseif (permitt_hf_model(curr_temp) > permitt_hf_max)
            permitt_hf_model(curr_temp) ...
                            = permitt_hf_max;
        end
        if do_dc(curr_temp)
            conduct_dc(curr_temp) ...
                            = 10 ^ -model_final{curr_temp}(end); % convert back to DC conductivity
            if (conduct_dc(curr_temp) < conduct_dc_min);
                conduct_dc(curr_temp) ...
                            = conduct_dc_min;
            elseif (conduct_dc(curr_temp) > conduct_dc_max);
                conduct_dc(curr_temp) ...
                            = conduct_dc_max;
            end
        else
            conduct_dc(curr_temp) ...
                            = NaN;
        end
        [permitt_diff{curr_temp}, freq_relax{curr_temp}, alpha_cole{curr_temp}, period_relax{curr_temp}] ...
                            = deal(NaN(num_relax(curr_temp), 1));
        for ii = 1:num_relax(curr_temp)
            [permitt_diff{curr_temp}(ii), freq_relax{curr_temp}(ii), alpha_cole{curr_temp}(ii), period_relax{curr_temp}(ii)] ...
                            = deal(model_final{curr_temp}(2 + ((ii - 1) * 3)), (1 / (2 * pi * model_final{curr_temp}(3 + ((ii - 1) * 3)))), model_final{curr_temp}(4 + ((ii - 1) * 3)), model_final{curr_temp}(3 + ((ii - 1) * 3)));
        end
        
        % reorder so first relaxation is highest frequency
        [~, relax_ord]      = sort(freq_relax{curr_temp});
        relax_ord           = flipud(relax_ord);
        if ~all(relax_ord' == 1:num_relax(curr_temp))
            [permitt_diff{curr_temp}, freq_relax{curr_temp}, alpha_cole{curr_temp}, period_relax{curr_temp}] ...
                            = deal(permitt_diff{curr_temp}(relax_ord), freq_relax{curr_temp}(relax_ord), alpha_cole{curr_temp}(relax_ord), period_relax{curr_temp}(relax_ord));
            model_final_old = model_final{curr_temp};
            for ii = 1:num_relax(curr_temp);
                [model_final{curr_temp}(2 + ((ii - 1) * 3)), model_final{curr_temp}(3 + ((ii - 1) * 3)), model_final{curr_temp}(4 + ((ii - 1) * 3))] ...
                            = deal(model_final_old(2 + ((relax_ord(ii) - 1) * 3)), model_final_old(3 + ((relax_ord(ii) - 1) * 3)), model_final_old(4 + ((relax_ord(ii) - 1) * 3)));
            end
        end
        
        % assign model uncertainties
        num_param(curr_temp)= length(model_final{curr_temp}); % number of model parameters for current temperature
        model_std{curr_temp}= NaN(num_param(curr_temp), 2);
        chisq_diff(curr_temp) ...
                            = delta_chisq(tol, ((2 * length(freq_trim{curr_temp})) - num_param(curr_temp))); % chi-squared difference given number of degrees of freedom
        err                 = [(err_rel(1) .* permitt_real_trim{curr_temp}); (err_rel(2) .* permitt_imag_trim{curr_temp})]; % relative data error
        res_def(curr_temp)  = chisq([permitt_real_trim{curr_temp}; permitt_imag_trim{curr_temp}], fm_cole(model_final{curr_temp}, freq_trim{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum), err);
        
        % test an increasingly large range about best-fit parameter
        for curr_param = 1:num_param(curr_temp)
            curr_range      = 1;
            while (any(isnan(model_std{curr_temp}(curr_param, :))) && (curr_range <= length(frac_test)))
                param_test
                curr_range  = curr_range + 1;
            end
        end
        
        % concatenate values for easier referencing
        for ii = 1:num_relax(curr_temp)
            [permitt_diff_cat(ii, curr_temp), freq_relax_cat(ii, curr_temp), alpha_cole_cat(ii, curr_temp), period_relax_cat(ii, curr_temp)] ...
                            = deal(permitt_diff{curr_temp}(ii), freq_relax{curr_temp}(ii), alpha_cole{curr_temp}(ii), period_relax{curr_temp}(ii));
        end
        if (num_relax(curr_temp) < 4) % do some NaNing if number of relaxations has changed
            [permitt_diff_cat((num_relax(curr_temp) + 1):4, curr_temp), freq_relax_cat((num_relax(curr_temp) + 1):4, curr_temp), alpha_cole_cat((num_relax(curr_temp) + 1):4, curr_temp), period_relax_cat((num_relax(curr_temp) + 1):4, curr_temp)] ...
                            = deal(NaN(length((num_relax(curr_temp) + 1):4), 1));
        end
        
        % extract standard deviations from model_std
        permitt_hf_model_std(curr_temp, :) ...
                            = model_std{curr_temp}(1, :);
        if do_dc(curr_temp)
            conduct_dc_std(curr_temp, :) ...
                            = fliplr(10 .^ -model_std{curr_temp}(end, :));
        else
            conduct_dc_std(curr_temp, :) ...
                            = NaN(1, 2);
        end
        [permitt_diff_std{curr_temp}, freq_relax_std{curr_temp}, alpha_cole_std{curr_temp}, period_relax_std{curr_temp}] ...
                            = deal(NaN(num_relax(curr_temp), 2));
        for ii = 1:num_relax(curr_temp)
            [permitt_diff_std{curr_temp}(ii, :), freq_relax_std{curr_temp}(ii, :), alpha_cole_std{curr_temp}(ii, :), period_relax_std{curr_temp}(ii, :)] ...
                            = deal(model_std{curr_temp}(2 + ((ii - 1) * 3), :), model_std{curr_temp}(3 + ((ii - 1) * 3), :), model_std{curr_temp}(4 + ((ii - 1) * 3), :), model_std{curr_temp}(3 + ((ii - 1) * 3), :));
            freq_relax_std{curr_temp}(ii, :) ...
                            = fliplr(1 ./ ((2 * pi) .* freq_relax_std{curr_temp}(ii, :)));
            [permitt_diff_std_cat{ii}(curr_temp, :), freq_relax_std_cat{ii}(curr_temp, :), alpha_cole_std_cat{ii}(curr_temp, :), period_relax_std_cat{ii}(curr_temp, :)] ...
                            = deal(permitt_diff_std{curr_temp}(ii, :), freq_relax_std{curr_temp}(ii, :), alpha_cole_std{curr_temp}(ii, :), period_relax_std{curr_temp}(ii, :));
        end
        
        inv_report % report inverse results to GUI screen
        
        do_fm % display final results
        do_activ
        plot_fits % plot fit parameters in separate sub-function
        if ~isempty(file_save)
            save_cfg_var
            save_var
        end
    end

%% calculate simple (all temperatures) activation energies of each relaxation frequency and DC conductivity

    function do_activ(source, eventdata)
        for ii = 1:max(num_relax)
            curr_ind    = find(~isnan(freq_relax_cat(ii, :)));
            if (length(curr_ind) > 1)
                activ_relax_poly{ii} ...
                        = polyfit(temp_mean_inv(curr_ind), log(freq_relax_cat(ii, curr_ind))', 1);
                activ_relax(ii) ...
                        = (-1e3 * boltzmann * activ_relax_poly{ii}(1)) / conv_eV_J;
            else
                activ_relax(ii) ...
                        = NaN;
            end
        end
        if (length(find(~isnan(conduct_dc))) > 1)
            activ_dc_poly ...
                        = polyfit(temp_mean_inv(~isnan(conduct_dc)), log(conduct_dc(~isnan(conduct_dc))), 1);
            activ_dc    = (-1e3 * boltzmann * activ_dc_poly(1)) / conv_eV_J;
        else
            [activ_dc_poly, activ_dc] ...
                        = deal(NaN);
        end
    end

%% report inverse results

    function inv_report(source, eventdata)
        
        figure(inv_gui)
        
        set(inv_edit, 'string', sprintf(relax_str{2}, res_model(curr_temp))) % display residual
        
        % put all fitted results in boxes
        set(permitt_hf_box, 'string', sprintf(relax_str{1}, permitt_hf_model(curr_temp)))
        set(permitt_hf_edit, 'string', sprintf(relax_str{1}, permitt_hf_model(curr_temp)))
        if ((permitt_hf_model(curr_temp) < permitt_hf_min) || (permitt_hf_model(curr_temp) > permitt_hf_max))
            set(permitt_hf_slide, 'min', (permitt_hf_model(curr_temp) - (0.5 * permitt_hf_model(curr_temp))), 'max', (permitt_hf_model(curr_temp) + (0.5 * permitt_hf_model(curr_temp))), 'sliderstep', [0.001 0.01])
        end
        set(permitt_hf_slide, 'value', permitt_hf_model(curr_temp))
        permitt_hf_guess    = permitt_hf_model(curr_temp);
        
        if do_dc(curr_temp)
            set(conduct_dc_box, 'string', sprintf(relax_str{2}, conduct_dc(curr_temp)))
            set(conduct_dc_edit, 'string', sprintf(relax_str{2}, conduct_dc(curr_temp)))
            if ((conduct_dc(curr_temp) < conduct_dc_min) || (conduct_dc(curr_temp) > conduct_dc_max))
                set(conduct_dc_slide, 'min', (log10(conduct_dc(curr_temp)) - (0.5 * abs(log10(conduct_dc(curr_temp))))), 'max', (log10(conduct_dc(curr_temp)) + (0.5 * abs(log10(conduct_dc(curr_temp))))))
            end
            set(conduct_dc_slide, 'value', log10(conduct_dc(curr_temp)))
            conduct_dc_guess= conduct_dc(curr_temp);
        end
        
        for ii = 1:num_relax(curr_temp)
            set(permitt_diff_box(ii), 'string', sprintf(relax_str{1}, permitt_diff{curr_temp}(ii)))
            set(permitt_diff_edit(ii), 'string', sprintf(relax_str{1}, permitt_diff{curr_temp}(ii)))
            if ((permitt_diff{curr_temp}(ii) < eval(['permitt_diff_min_' num2str(ii)])) || (permitt_diff{curr_temp}(ii) > eval(['permitt_diff_max_' num2str(ii)])))
                set(permitt_diff_slide(ii), 'min', (permitt_diff{curr_temp}(ii) - (0.5 * permitt_diff{curr_temp}(ii))), 'max', (permitt_diff{curr_temp}(ii) + (0.5 * permitt_diff{curr_temp}(ii))), 'sliderstep', [0.001 0.01])
            end
            set(permitt_diff_slide(ii), 'value', permitt_diff{curr_temp}(ii))
            relax_guess(ii, 1) ...
                            = permitt_diff{curr_temp}(ii);
            set(freq_relax_box(ii), 'string', sprintf(relax_str{2}, freq_relax{curr_temp}(ii)))
            set(freq_relax_edit(ii), 'string', sprintf(relax_str{2}, freq_relax{curr_temp}(ii)))
            if ((freq_relax{curr_temp}(ii) < eval(['freq_relax_min_' num2str(ii)])) || (freq_relax{curr_temp}(ii) > eval(['freq_relax_max_' num2str(ii)])))
                set(freq_relax_slide(ii), 'min', (log10(freq_relax{curr_temp}(ii)) - (0.5 * abs(log10(freq_relax{curr_temp}(ii))))), 'max', ...
                    (log10(freq_relax{curr_temp}(ii)) + (0.5 * abs(log10(freq_relax{curr_temp}(ii))))))
            end
            set(freq_relax_slide(ii), 'value', log10(freq_relax{curr_temp}(ii)))
            relax_guess(ii, 2) ...
                            = freq_relax{curr_temp}(ii);
            set(alpha_cole_box(ii), 'string', sprintf(relax_str{3}, alpha_cole{curr_temp}(ii)))
            set(alpha_cole_edit(ii), 'string', sprintf(relax_str{3}, alpha_cole{curr_temp}(ii)))
            if ((alpha_cole{curr_temp}(ii) < eval(['alpha_cole_min_' num2str(ii)])) || (alpha_cole{curr_temp}(ii) > eval(['alpha_cole_max_' num2str(ii)])))
                set(alpha_cole_slide(ii), 'min', (alpha_cole{curr_temp}(ii) - (0.5 * alpha_cole{curr_temp}(ii))), 'max', (alpha_cole{curr_temp}(ii) + (0.5 * alpha_cole{curr_temp}(ii))), 'sliderstep', [0.001 0.01])
            end
            set(alpha_cole_slide(ii), 'value', alpha_cole{curr_temp}(ii))
            relax_guess(ii, 3) ...
                            = alpha_cole{curr_temp}(ii);
            set(period_relax_edit(ii), 'string', num2str(1 / (2 * pi * relax_guess(ii, 2))))
        end
        
    end

%% plot parameters fits in separate figure

    function plot_fits(source, eventdata)
        figure(arr_fig)
        delete(p_arr_fit(ishandle(p_arr_fit(:))))
        delete(p_arr_err(ishandle(p_arr_err(:))))
        delete(p_arr_model(ishandle(p_arr_model(:))))
        delete(p_arr_leg(ishandle(p_arr_leg(:))))
        if ishandle(pdcf)
            delete(pdcf)
        end
        if ishandle(pdc)
            delete(pdc)
        end
        for ii = 1:6
            axes(axf(ii))
            lr              = {};
            switch ii
                case 1
                    for jj = 1:max(num_relax)
                        curr_ind ...
                            = find(~isnan(permitt_diff_cat(jj, :)));
                        p_arr_fit(jj, ii) ...
                            = plot(temp_mean_inv(curr_ind), permitt_diff_cat(jj, curr_ind), 's-', 'linewidth', 2, 'color', colors_fit{jj}, 'markerfacecolor', colors_fit{jj});
                        p_arr_err(jj, ii) ...
                            = errorbar(temp_mean_inv(curr_ind), permitt_diff_cat(jj, curr_ind), (permitt_diff_cat(jj, curr_ind)' - permitt_diff_std_cat{jj}(curr_ind, 1)), (permitt_diff_std_cat{jj}(curr_ind, 2) - permitt_diff_cat(jj, curr_ind)'), 'color', colors_fit{jj});
                    end
                case 2
                    for jj = 1:max(num_relax)
                        curr_ind ...
                            = find(~isnan(freq_relax_cat(jj, :)));
                        if (length(curr_ind) > 1)
                            p_arr_model(jj, ii) ...
                            = semilogy(temp_mean_inv(curr_ind), exp(polyval(activ_relax_poly{jj}, temp_mean_inv(curr_ind))), '--', 'linewidth', 2, 'color', colors_fit{jj});
                            lr  ...
                            = [lr; [sprintf(relax_str{4}, activ_relax(jj)) ' eV']]; %#ok<AGROW>
                        
                        end
                        p_arr_fit(jj, ii) ...
                            = semilogy(temp_mean_inv(curr_ind), freq_relax_cat(jj, curr_ind), 's', 'color', colors_fit{jj}, 'markerfacecolor', colors_fit{jj});
                        p_arr_err(jj, ii) ...
                            = errorbar(temp_mean_inv(curr_ind), freq_relax_cat(jj, curr_ind), (freq_relax_cat(jj, curr_ind)' - freq_relax_std_cat{jj}(curr_ind, 1)), (freq_relax_std_cat{jj}(curr_ind, 2) - freq_relax_cat(jj, curr_ind)'), 'color', colors_fit{jj}, 'linestyle', 'none');
                    end
                case 3
                    for jj = 1:max(num_relax)
                        curr_ind ...
                            = find(~isnan(alpha_cole_cat(jj, :)));
                        p_arr_fit(jj, ii) ...
                            = plot(temp_mean_inv(curr_ind), alpha_cole_cat(jj, curr_ind), 's-', 'linewidth', 2, 'color', colors_fit{jj}, 'markerfacecolor', colors_fit{jj});
                        p_arr_err(jj, ii) ...
                            = errorbar(temp_mean_inv(curr_ind), alpha_cole_cat(jj, curr_ind), (alpha_cole_cat(jj, curr_ind)' - alpha_cole_std_cat{jj}(curr_ind, 1)), (alpha_cole_std_cat{jj}(curr_ind, 2) - alpha_cole_cat(jj, curr_ind)'), 'color', colors_fit{jj});
                    end
                case 4
                    curr_ind ...
                            = find(~isnan(permitt_hf_model));
                    p_arr_fit(1, ii) ...
                            = plot(temp_mean_inv(curr_ind), permitt_hf(curr_ind), 'ks-', 'markerfacecolor', 'k');
                    p_arr_model(1, ii) ...
                            = plot(temp_mean_inv(curr_ind), permitt_hf_model(curr_ind), 'ro--', 'linewidth', 2, 'markerfacecolor', 'r');
                    p_arr_err(1, ii) ...
                            = errorbar(temp_mean_inv(curr_ind), permitt_hf_model(curr_ind), (permitt_hf_model(curr_ind) - permitt_hf_model_std(curr_ind, 1)), (permitt_hf_model_std(curr_ind, 2) - permitt_hf_model(curr_ind)), 'r');
                case 5
                    curr_ind ...
                            = find(~isnan(conduct_hf_model));
                    p_arr_fit(1, ii) ...
                            = semilogy(temp_mean_inv(curr_ind), conduct_hf(curr_ind), 'ks', 'linewidth', 2, 'markerfacecolor', 'k');
                    p_arr_model(1, ii) ...
                            = semilogy(temp_mean_inv(curr_ind), conduct_hf_model(curr_ind), 'ro--', 'linewidth', 2, 'markerfacecolor', 'r');
                    if any(do_dc) % only show dc conductivity if it got modeled
                        curr_ind ...
                            = find(~isnan(conduct_dc));
                        if (length(curr_ind) > 1)
                            pdcf ...
                            = semilogy(temp_mean_inv(curr_ind), exp(polyval(activ_dc_poly, temp_mean_inv(curr_ind))), 'r--', 'linewidth', 2);
                            lr ...
                            = [sprintf(relax_str{4}, activ_dc) ' eV'];
                        end
                        pdc = semilogy(temp_mean_inv(curr_ind), conduct_dc(curr_ind), 'ks', 'linewidth', 2, 'markerfacecolor', 'k');
                        p_arr_err(1, ii) ...
                            = errorbar(temp_mean_inv(curr_ind), conduct_dc(curr_ind), (conduct_dc(curr_ind) - conduct_dc_std(curr_ind, 1)), (conduct_dc_std(curr_ind, 2) - conduct_dc(curr_ind)), 'k');
                    else
                        curr_ind ...
                            = 0;
                    end
                case 6
                    p_arr_fit(1, ii) ...
                            = plot(temp_mean_inv(~isnan(res_model)), res_model(~isnan(res_model)), 'ks-', 'markerfacecolor', 'k');
            end
            if (num_temp > 1)
                xlim([min(temp_mean_inv) max(temp_mean_inv)])
            else
                axis tight
            end
            switch ii
                case 2
                    if any(ishandle(p_arr_model(:, ii)))
                        p_arr_leg(ii) ...
                            = legend(p_arr_model(ishandle(p_arr_model(:))), lr, 'location', 'southwest');
                    end
                case 4
                    p_arr_leg(ii)  ...
                            = legend([p_arr_fit(1, ii) p_arr_model(1, ii)], 'observed', 'modeled');
                case 5
                    if (length(curr_ind) > 1)
                        p_arr_leg(ii) ...
                            = legend(pdcf, lr);
                    end
            end
        end
    end

%% save configuration file

    function save_cfg(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        if ~isempty(path_curr)
            [file_save, path_curr] ...
                            = uiputfile('*.mat', 'Save configuration:', [path_curr file_merge(1:(end - 9)) '_inv_cfg.mat']);
        else
            [file_save, path_curr] ...
                            = uiputfile('*.mat', 'Save configuration:', [file_merge(1:(end - 9)) '_inv_cfg.mat']);
        end
        if ~file_save
            [file_save, path_curr] ...
                            = deal('');
        end
        if ~isempty(file_save)
            save_cfg_var
        end
    end

%% save configuration variables

    function save_cfg_var(source, eventdata)
        save([path_curr file_save], 'freq_min', 'freq_max', 'freq_trim', 'do_weight', 'weight_freq', 'num_relax', 'do_dc', 'relax_guess', 'model_guess', 'permitt_hf_guess', 'conduct_dc_guess', 'model_final', 'permitt_hf', 'permitt_hf_model', 'permitt_hf_model_std', 'conduct_hf', ...
                                    'conduct_hf_model', 'conduct_dc', 'conduct_dc_std', 'permitt_real_model', 'permitt_imag_model', 'conduct_model', 'phase_model', 'res_model', 'freq_relax', 'freq_relax_cat', 'freq_relax_std', 'freq_relax_std_cat', 'alpha_cole', 'alpha_cole_cat', ...
                                    'alpha_cole_std', 'alpha_cole_std_cat', 'permitt_diff', 'permitt_diff_cat', 'permitt_diff_std', 'permitt_diff_std_cat', 'permitt_real_trim', 'permitt_imag_trim', 'activ_relax', 'activ_relax_poly', 'activ_dc', 'activ_dc_poly', 'period_relax', ...
                                    'period_relax_std', 'period_relax_cat', 'period_relax_std_cat')
    end

%% save inverse results

    function save_inv(source, eventdata) % save model results
        save_cfg
        if ~isempty(path_curr)
            [file_save, path_curr] ...
                            = uiputfile('*.mat', 'Save inverse modeling:', [path_curr file_merge(1:(end - 4)) '_inv.mat']);
        else
            [file_save, path_curr] ...
                            = uiputfile('*.mat', 'Save inverse modeling:', [file_merge(1:(end - 4)) '_inv.mat']);
        end
        if ~file_save
            [file_save, path_curr] ...
                            = deal('');
        end
        if ~isempty(file_save)
            save_var
        end
    end
    
%% save inverse variables

    function save_var(source, eventdata)
        save([path_curr file_save], 'num_temp', 'temp_mean', 'temp_std', 'freq_trim', 'do_weight', 'weight_freq', 'num_relax', 'do_dc', 'model_final', 'permitt_hf', 'permitt_hf_model', 'permitt_hf_model_std', 'conduct_hf', 'conduct_hf_model', 'conduct_dc', 'conduct_dc_std', ...
                                    'permitt_real_model', 'permitt_imag_model', 'conduct_model', 'phase_model', 'res_model', 'freq_relax', 'freq_relax_cat', 'freq_relax_std', 'freq_relax_std_cat', 'alpha_cole', 'alpha_cole_cat', 'alpha_cole_std', 'alpha_cole_std_cat', 'permitt_diff', ...
                                    'permitt_diff_cat', 'permitt_diff_std', 'permitt_diff_std_cat', 'permitt_real_trim', 'permitt_imag_trim', 'activ_relax', 'activ_relax_poly', 'activ_dc', 'activ_dc_poly', 'period_relax', 'period_relax_std', 'period_relax_cat', 'period_relax_std_cat')
    end

%% sliders, sub-functions etc

    function slide_freq_min(source, eventdata) % update min/max frequency for fitting
        if ~focus_check
            focus_check     = true;
        end
        freq_min            = 10 ^ get(freq_min_slide, 'value');
        set(freq_min_edit, 'string', sprintf(relax_str{2}, freq_min))
        if any(ishandle(p_freq_min))
            delete(p_freq_min(ishandle(p_freq_min)))
        end
        for ii = 1:4
            axes(ax(ii))
            if (ii < 4)
                p_freq_min(ii) ...
                            = loglog(repmat(freq_min, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
            else
                p_freq_min(ii) ...
                            = semilogx(repmat(freq_min, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
            end
        end
        tmp1                = dbstack;
        if (data_loaded && (length(tmp1) == 1)) % only want to redo forward model if actually adjusting frequency
            do_fm
        end
    end
%%
    function slide_freq_max(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        freq_max            = 10 ^ (get(freq_max_slide, 'value'));
        set(freq_max_edit, 'string', sprintf(relax_str{2}, freq_max));
        if any(ishandle(p_freq_max))
            delete(p_freq_max(ishandle(p_freq_max)));
        end
        for ii = 1:4
            axes(ax(ii))
            if (ii < 4)
                p_freq_max(ii) ...
                            = loglog(repmat(freq_max, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
            else
                p_freq_max(ii) ...
                            = semilogx(repmat(freq_max, 1, 2), get(gca, 'ylim'), 'k--', 'linewidth', 2);
            end
        end
        tmp1                = dbstack;
        if (data_loaded && (length(tmp1) == 1));
            do_fm
        end
    end
%%
    function check_dc(source, eventdata) % update checkbox for DC conductivity
        if ~focus_check
            focus_check     = true;
        end
        do_dc(curr_temp)    = logical(get(conduct_dc_check, 'value'));
        if data_loaded
            do_fm
        end
        do_dc(isnan(permitt_hf_model)) ...
                            = do_dc(curr_temp);
    end
%%
    function check_weight(source, eventdata) % update checkbox for frequency weighting
        if ~focus_check
            focus_check     = true;
        end
        do_weight           = logical(get(weight_check, 'value')); 
    end
%%
    function check_sep(source, eventdata) % update checkbox for inversion separation
        if ~focus_check
            focus_check     = true;
        end
        do_sep              = logical(get(sep_check, 'value'));
        if ~do_sep
            if any(ishandle(p_model_sep(:)))
                delete(p_model_sep(ishandle(p_model_sep(:))))
            end
        else
            disp_sep
        end
    end
%%
    function relax_radio(~, eventdata) % update number of relaxations
        if ~focus_check
            focus_check     = true;
        end
        num_relax(curr_temp)= str2double(get(eventdata.NewValue, 'string'));
        if data_loaded
            do_fm
        end
        num_relax(isnan(permitt_hf_model)) ...
                            = num_relax(curr_temp);
    end
%%
    function slide_permitt_hf(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        permitt_hf_guess    = get(permitt_hf_slide, 'value');
        set(permitt_hf_edit, 'string', sprintf(relax_str{1}, permitt_hf_guess))
        if data_loaded
            do_fm
        end
    end
%%
    function slide_conduct_dc(source, eventdata)
        conduct_dc_guess    = 10 ^ get(conduct_dc_slide, 'value');
        set(conduct_dc_edit, 'string', sprintf(relax_str{2}, conduct_dc_guess))
        if (data_loaded && (do_dc(curr_temp)))
            do_fm
        end
    end
%%
    function slide_permitt_diff_1(source, eventdata) %#ok<*DEFNU> % update permitt_diff sliders
        tmp1                = 1;
        slide_permitt_diff
    end

    function slide_permitt_diff_2(source, eventdata)
        tmp1                = 2;
        slide_permitt_diff
    end

    function slide_permitt_diff_3(source, eventdata)
        tmp1                = 3;
        slide_permitt_diff
    end

    function slide_permitt_diff_4(source, eventdata)
        tmp1                = 4;
        slide_permitt_diff
    end

    function slide_permitt_diff(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        relax_guess(tmp1, 1)= get(permitt_diff_slide(tmp1), 'value');
        set(permitt_diff_edit(tmp1), 'string', sprintf(relax_str{1}, relax_guess(tmp1, 1)))
        if (data_loaded && (num_relax(curr_temp) >= tmp1))
            do_fm
        end
    end
%%
    function slide_freq_relax_1(source, eventdata)
        tmp1                = 1;
        slide_freq_relax
    end

    function slide_freq_relax_2(source, eventdata)
        tmp1                = 2;
        slide_freq_relax
    end

    function slide_freq_relax_3(source, eventdata)
        tmp1                = 3;
        slide_freq_relax
    end

    function slide_freq_relax_4(source, eventdata)
        tmp1                = 4;
        slide_freq_relax
    end

    function slide_freq_relax(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        relax_guess(tmp1, 2)= 10 ^ get(freq_relax_slide(tmp1), 'value');
        set(freq_relax_edit(tmp1), 'string', sprintf(relax_str{2}, relax_guess(tmp1, 2)))
        set(period_relax_edit(tmp1), 'string', num2str(1 / (2 * pi * relax_guess(tmp1, 2))))
        if (data_loaded && (num_relax(curr_temp) >= tmp1))
            do_fm
        end
    end
%%
    function slide_alpha_cole_1(source, eventdata)
        tmp1                = 1;
        slide_alpha_cole
    end

    function slide_alpha_cole_2(source, eventdata)
        tmp1                = 2;
        slide_alpha_cole
    end

    function slide_alpha_cole_3(source, eventdata)
        tmp1                = 3;
        slide_alpha_cole
    end

    function slide_alpha_cole_4(source, eventdata)
        tmp1                = 4;
        slide_alpha_cole
    end

    function slide_alpha_cole(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        relax_guess(tmp1, 3)= get(alpha_cole_slide(tmp1), 'value');
        set(alpha_cole_edit(tmp1), 'string', sprintf(relax_str{3}, relax_guess(tmp1, 3)))
        if (data_loaded && (num_relax(curr_temp) >= tmp1))
            do_fm
        end
    end
%%
    function weight(source, eventdata)
        weight_freq{curr_temp} ...
                            = [(max(permitt_real_trim{curr_temp}) ./ permitt_real_trim{curr_temp}) (max(permitt_imag_trim{curr_temp}) ./ permitt_imag_trim{curr_temp})];
        weight_freq{curr_temp}((freq_trim{curr_temp} >= freq_repeat(1)) & (freq_trim{curr_temp} < freq_repeat(2)), :) ...
                            = weight_freq{curr_temp}(((freq_trim{curr_temp} >= freq_repeat(1)) & (freq_trim{curr_temp} < freq_repeat(2))), :) ./ 2;
        weight_freq{curr_temp}(freq_trim{curr_temp} >= freq_repeat(2), :) ...
                            = weight_freq{curr_temp}((freq_trim{curr_temp} >= freq_repeat(2)), :) ./ 4;
        weight_freq{curr_temp} ...
                            = reshape(weight_freq{curr_temp}, numel(weight_freq{curr_temp}), 1);
    end
%%    
    function param_test(source, eventdata)
        param_tmp           = linspace((model_final{curr_temp}(curr_param) - (frac_test(curr_range) * model_final{curr_temp}(curr_param))), (model_final{curr_temp}(curr_param) + (frac_test(curr_range) * model_final{curr_temp}(curr_param))));
        model_tmp           = model_final{curr_temp};
        res                 = NaN(1, 100);
        for ii = 1:100
            model_tmp(curr_param) ...
                            = param_tmp(ii);
            res(ii)         = chisq([permitt_real_trim{curr_temp}; permitt_imag_trim{curr_temp}], fm_cole(model_tmp, freq_trim{curr_temp}, num_relax(curr_temp), do_dc(curr_temp), permitt_vacuum), err); ...
                % chi-squared of original model
        end
        try %#ok<TRYNC>
            if isnan(model_std{curr_temp}(curr_param, 1))
                model_std{curr_temp}(curr_param, 1) ...
                            = model_final{curr_temp}(curr_param) + interp1((res(1:50) - res_def(curr_temp)), (param_tmp(1:50) - model_final{curr_temp}(curr_param)), chisq_diff(curr_temp));
            end
            if isnan(model_std{curr_temp}(curr_param, 2))
                model_std{curr_temp}(curr_param, 2) ...
                            = model_final{curr_temp}(curr_param) + interp1((res(51:end) - res_def(curr_temp)), (param_tmp(51:end) - model_final{curr_temp}(curr_param)), chisq_diff(curr_temp));
            end
        end
    end
%%
    function nuke_plot(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        for ii = 1:4
            axes(ax(ii))
            delete(get(ax(ii), 'children'))
            if (ii < 4)
                p_data(ii, 1) ...
                            = loglog(freq{curr_temp}, eval(plots1{ii}), 'r', 'linewidth', 2);
                p_data(ii, 2) ...
                            = loglog(freq{curr_temp}, eval(plots1{ii}), 'ko', 'markerfacecolor', 'r', 'markersize', 8);
                set(ax(ii), 'ylim', [min(eval(plots1{ii})) max(eval(plots1{ii}))])
            else
                p_data(ii, 1) ...
                            = semilogx(freq{curr_temp}, eval(plots1{ii}), 'r', 'linewidth', 2);
                p_data(ii, 2) ...
                            = semilogx(freq{curr_temp}, eval(plots1{ii}), 'ko', 'markerfacecolor', 'r', 'markersize', 8);
            end
        end
        slide_freq_min
        slide_freq_max
        if data_loaded
            do_fm
        end
    end
%%
    function nuke_inv(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        if inv_done(curr_temp)
            [model_guess{curr_temp}, model_final{curr_temp}, res_model(curr_temp), permitt_hf_model(curr_temp), conduct_dc(curr_temp), permitt_diff{curr_temp}, freq_relax{curr_temp}, alpha_cole{curr_temp}, num_param(curr_temp), model_std{curr_temp}, chisq_diff(curr_temp), res_def(curr_temp), ...
             permitt_diff_std{curr_temp}, permitt_hf(curr_temp), conduct_hf(curr_temp), freq_relax_std{curr_temp}, alpha_cole_std{curr_temp}] ...
                            = deal(NaN);
            conduct_dc_std(curr_temp, :) ...
                            = NaN(1, 2);
            [permitt_diff_cat(:, curr_temp), freq_relax_cat(:, curr_temp), alpha_cole_cat(:, curr_temp)] ...
                            = deal(NaN(4, 1));
            for ii = 1:4
                [permitt_diff_std_cat{ii}(curr_temp, :), freq_relax_std_cat{ii}(curr_temp, :), alpha_cole_std_cat{ii}(curr_temp, :)] ...
                            = deal(NaN(1, 2));
            end
            num_inv         = num_inv - 1;
            inv_done(curr_temp) ...
                            = false;
            set(inv_done_edit, 'string', [num2str(num_inv) ' /'])
            temp_str{curr_temp} ...
                            = temp_str{curr_temp}(1:(end - 1));
            set(temp_box, 'string', temp_str)
        end
        set([inv_edit permitt_hf_box], 'string', '')
        if do_dc(curr_temp)
            set(conduct_dc_box, 'string', '')
        end
        set([permitt_diff_box freq_relax_box alpha_cole_box], 'string', '')
        if num_inv
            do_activ
            plot_fits
        end
    end
%%
    function do_next(source, eventdata)
        if ~focus_check
            focus_check     = true;
        end
        if (curr_temp < num_temp)
            curr_temp       = curr_temp + 1;
            set(temp_box, 'value', curr_temp);
            plot_data
            if ~inv_done(curr_temp)
                do_inv
            end
        end
    end
%%
    function disp_sep(source, eventdata)
        % plot separated forward model
        for ii = 1:3
            axes(ax(ii))
            for jj = 1:length(permitt_comp_model_sep)
                if ((jj == length(permitt_comp_model_sep)) && do_dc(curr_temp) && (ii ~= 1))
                    p_model_sep(ii, jj) ...
                            = loglog(freq{curr_temp}, eval(plots3{ii}), 'color', colors_sep{end}, 'linewidth', 1);
                else
                    p_model_sep(ii, jj) ...
                            = loglog(freq{curr_temp}, eval(plots3{ii}), 'color', colors_sep{jj}, 'linewidth', 1);
                end
            end
        end
    end
%%
    function lim_permitt_diff_1(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(1, 1);
        lim_gen
    end

    function lim_permitt_diff_2(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(1, 2);
        lim_gen
    end

    function lim_permitt_diff_3(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(1, 3);
        lim_gen
    end

    function lim_permitt_diff_4(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(1, 4);
        lim_gen
    end

    function lim_freq_relax_1(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(2, 1);
        lim_gen
    end

    function lim_freq_relax_2(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(2, 2);
        lim_gen
    end

    function lim_freq_relax_3(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(2, 3);
        lim_gen
    end

    function lim_freq_relax_4(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(2, 4);
        lim_gen
    end

    function lim_alpha_cole_1(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(3, 1);
        lim_gen
    end

    function lim_alpha_cole_2(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(3, 2);
        lim_gen
    end

    function lim_alpha_cole_3(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(3, 3);
        lim_gen
    end

    function lim_alpha_cole_4(source, eventdata)
        [curr_relax_var, curr_slide] ...
                            = deal(3, 4);
        lim_gen
    end

    function lim_gen(source, eventdata)
        dia                 = dialog;
        annotation('textbox', [0.4 0.75 0.1 0.1], 'string', [relax_labels{curr_relax_var} '(' num2str(curr_slide) ') min'], 'fontsize', 16, 'edgecolor', 'none')
        tmp1                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.6 0.2 0.1], 'fontsize', 16, 'string', num2str(eval([relax_simple{curr_relax_var} '_min_' num2str(curr_slide)])), 'callback', eval(['@adj_' relax_simple{curr_relax_var} '_min']));
        annotation('textbox', [0.4 0.35 0.1 0.1], 'string', [relax_labels{curr_relax_var} '(' num2str(curr_slide) ') max'], 'fontsize', 16, 'edgecolor', 'none')
        tmp2                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.2 0.2 0.1], 'fontsize', 16, 'string', num2str(eval([relax_simple{curr_relax_var} '_max_' num2str(curr_slide)])), 'callback', eval(['@adj_' relax_simple{curr_relax_var} '_max']));
    end
%%
    function lim_permitt_hf(source, eventdata)
        dia                 = dialog;
        annotation('textbox', [0.4 0.75 0.1 0.1], 'string', '\epsilon''_{HF} min', 'fontsize', 16, 'edgecolor', 'none')
        tmp1                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.6 0.2 0.1], 'fontsize', 16, 'string', num2str(permitt_hf_min), 'callback', @adj_permitt_hf_min);
        annotation('textbox', [0.4 0.35 0.1 0.1], 'string', '\epsilon''_{HF} max', 'fontsize', 16, 'edgecolor', 'none')
        tmp2                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.2 0.2 0.1], 'fontsize', 16, 'string', num2str(permitt_hf_max), 'callback', @adj_permitt_hf_max);
    end
%%
    function lim_conduct_dc(source, eventdata)
        dia                 = dialog;
        annotation('textbox', [0.4 0.75 0.1 0.1], 'string', '\sigma_{DC} min', 'fontsize', 16, 'edgecolor', 'none')
        tmp1                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.6 0.2 0.1], 'fontsize', 16, 'string', num2str(conduct_dc_min), 'callback', @adj_conduct_dc_min);
        annotation('textbox', [0.4 0.35 0.1 0.1], 'string', '\sigma_{DC} max', 'fontsize', 16, 'edgecolor', 'none')
        tmp2                = uicontrol(dia, 'style', 'edit', 'units', 'normalized', 'position', [0.4 0.2 0.2 0.1], 'fontsize', 16, 'string', num2str(conduct_dc_max), 'callback', @adj_conduct_dc_max);
    end
%%
    function adj_permitt_diff_min(source, eventdata)
        eval(['permitt_diff_min_' num2str(curr_slide) ' = ' get(tmp1, 'string') ';']);
        model_min_abs(2 + ((curr_slide - 1) * 3)) ...
                            = eval(['permitt_diff_min_' num2str(curr_slide)]);
        set(permitt_diff_slide(curr_slide), 'min', eval(['permitt_diff_min_' num2str(curr_slide)]))
    end
%%
    function adj_freq_relax_min(source, eventdata)
        eval(['freq_relax_min_' num2str(curr_slide) ' = ' get(tmp1, 'string') ';']);
        model_min_abs(3 + ((curr_slide - 1) * 3)) ...
                            = eval(['freq_relax_min_' num2str(curr_slide)]);
        set(freq_relax_slide(curr_slide), 'min', log10(eval(['freq_relax_min_' num2str(curr_slide)])))
    end
%%
    function adj_alpha_cole_min(source, eventdata)
        eval(['alpha_cole_min_' num2str(curr_slide) ' = ' get(tmp1, 'string') ';']);
        model_min_abs(4 + ((curr_slide - 1) * 3)) ...
                            = eval(['alpha_cole_min_' num2str(curr_slide)]);
        set(alpha_cole_slide(curr_slide), 'min', eval(['alpha_cole_min_' num2str(curr_slide)]))
    end
%%
    function adj_permitt_hf_min(source, eventdata)
        [permitt_hf_min, model_min_abs(1)] ...
                            = deal(str2double(get(tmp1, 'string')));
        set(permitt_hf_slide, 'min', permitt_hf_min)
    end
%%
    function adj_conduct_dc_min(source, eventdata)
        conduct_dc_min      = str2double(get(tmp1, 'string'));
        set(conduct_dc_slide, 'min', log10(conduct_dc_min))
    end
%%
    function adj_permitt_diff_max(source, eventdata)
        eval(['permitt_diff_max_' num2str(curr_slide) ' = ' get(tmp2, 'string') ';']);
        model_max_abs(2 + ((curr_slide - 1) * 3)) ...
                            = eval(['permitt_diff_max_' num2str(curr_slide)]);
        set(permitt_diff_slide(curr_slide), 'max', eval(['permitt_diff_max_' num2str(curr_slide)]))
    end
%%
    function adj_freq_relax_max(source, eventdata)
        eval(['freq_relax_max_' num2str(curr_slide) ' = ' get(tmp2, 'string') ';']);
        model_max_abs(3 + ((curr_slide - 1) * 3)) ...
                            = eval(['freq_relax_max_' num2str(curr_slide)]);
        set(freq_relax_slide(curr_slide), 'max', log10(eval(['freq_relax_max_' num2str(curr_slide)])))
    end
%%
    function adj_alpha_cole_max(source, eventdata)
        eval(['alpha_cole_max_' num2str(curr_slide) ' = ' get(tmp2, 'string') ';']);
        model_max_abs(4 + ((curr_slide - 1) * 3)) ...
                            = eval(['alpha_cole_max_' num2str(curr_slide)]);
        set(alpha_cole_slide(curr_slide), 'max', eval(['alpha_cole_min_' num2str(curr_slide)]))
    end
%%
    function adj_permitt_hf_max(source, eventdata)
        [permitt_hf_max, model_max_abs(1)] ...
                            = deal(str2double(get(tmp2, 'string')));
        set(permitt_hf_slide, 'max', permitt_hf_max)
    end
%%
    function adj_conduct_dc_max(source, eventdata)
        conduct_dc_max      = str2double(get(tmp2, 'string'));
        set(conduct_dc_slide, 'max', log10(conduct_dc_max))
    end

%% Keyboard shortcuts for various functions

    function keypress(~, eventdata)
        switch eventdata.Key
            case '1'
                set(relax_group, 'selectedobject', relax_check(1))
                num_relax(curr_temp) ...
                            = 1;
                if data_loaded
                    do_fm
                end
                num_relax(isnan(permitt_hf_model)) ...
                            = num_relax(curr_temp);
            case '2'
                set(relax_group, 'selectedobject', relax_check(2))
                num_relax(curr_temp) ...
                            = 2;
                if data_loaded
                    do_fm
                end
                num_relax(isnan(permitt_hf_model)) ...
                            = num_relax(curr_temp);
            case '3'
                set(relax_group, 'selectedobject', relax_check(3))
                num_relax(curr_temp) ...
                            = 3;
                if data_loaded
                    do_fm
                end
                num_relax(isnan(permitt_hf_model)) ...
                            = num_relax(curr_temp);
            case '4'
                set(relax_group, 'selectedobject', relax_check(4))
                num_relax(curr_temp) ...
                            = 4;
                if data_loaded
                    do_fm
                end
                num_relax(isnan(permitt_hf_model)) ...
                            = num_relax(curr_temp);
            case 'c'
                load_cfg
            case 'd'
                if get(conduct_dc_check, 'value')
                    set(conduct_dc_check, 'value', 0)
                else
                    set(conduct_dc_check, 'value', 1)
                end
                check_dc
            case 'e'
                nuke_inv
            case 'i'
                do_inv
            case 'l'
                load_merge
            case 'n'
                do_next
            case 'o'
                save_cfg
            case 'p'
                if get(sep_check, 'value')
                    set(sep_check, 'value', 0)
                else
                    set(sep_check, 'value', 1)
                end
                check_sep
            case 'r'
                nuke_plot
            case 's'
                save_inv
            case 'w'
                if get(weight_check, 'value')
                    set(weight_check, 'value', 0)
                else
                    set(weight_check, 'value', 1)
                end
                check_weight
            case 'uparrow'
                if (curr_temp > 1)
                    curr_temp ...
                            = curr_temp - 1;
                    set(temp_box, 'value', curr_temp)
                    plot_data
                end
            case 'downarrow'
                if (curr_temp < num_temp)
                    curr_temp ...
                            = curr_temp + 1;
                    set(temp_box, 'value', curr_temp)
                    plot_data
                end
            case 'space'
                figure(arr_fig)
        end
    end

    function keypress2(~, eventdata)
        switch eventdata.Key
            case 'space'
                figure(inv_gui)
        end
    end

%% Test function

    function misctest(source, eventdata)
        
        disp('Test done.')
    end
end
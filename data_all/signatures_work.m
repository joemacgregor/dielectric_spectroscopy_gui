% SIGNATURES_WORK Analysis and figures for dielectric signatures paper.
% 
% Joe MacGregor
% Last updated: 10/12/15

clear

plotting                    = false;

permitt_vacuum              = 8.8541878176e-12; % permittivity of the vacuum, F/m
boltzmann                   = 1.3806503e-23; % Boltzmann constant
conv_eV_J                   = 1.60217646e-19; % J/eV


name_core                   = dir('data/signatures_fig1/');
name_core                   = {name_core(4:end).name};

num_core                    = length(name_core);
num_core2                   = num_core + 2;

core                        = cell(1, num_core);
for ii = 1:num_core
    core{ii}                = load(['data/signatures_fig1/' name_core{ii}]);
end

% fix one of the Mullins samples' depths
core{3}.depth(1)            = 8.57;

core{num_core + 1}          = load('data/Glen_Paren_1975');
core{num_core + 2}          = load('data/Reynolds_1985_Fig3b.mat');
core{num_core + 2}.permitt_real ...
                            = NaN(length(core{num_core + 2}.permitt_imag), 1);
core{num_core + 1}.conduct  = (2 * pi * permitt_vacuum) .* core{num_core + 1}.freq .* core{num_core + 1}.permitt_imag;
core{num_core + 2}.conduct  = (2 * pi * permitt_vacuum) .* core{num_core + 2}.freq .* core{num_core + 2}.permitt_imag;

[core{num_core + 1}.freq, core{num_core + 1}.permitt_real, core{num_core + 1}.permitt_imag, core{num_core + 1}.conduct] ...
                            = deal(core{num_core + 1}.freq(3:(end - 4)), core{num_core + 1}.permitt_real(3:(end - 4)), core{num_core + 1}.permitt_imag(3:(end - 4)), core{num_core + 1}.conduct(3:(end - 4)));
ind_40                      = find(core{num_core + 2}.temp == -40.9);
[core{num_core + 2}.freq, core{num_core + 2}.temp, core{num_core + 2}.permitt_real, core{num_core + 2}.permitt_imag, core{num_core + 2}.conduct] ...
                            = deal(core{num_core + 2}.freq(ind_40), core{num_core + 2}.temp(ind_40), core{num_core + 2}.permitt_real(ind_40), core{num_core + 2}.permitt_imag(ind_40), core{num_core + 2}.conduct(ind_40));
                        
% core names to load
name_core_long              = {'Upper Fremont Gl.' 'GISP2' 'Mullins Gl.' 'Newall Gl.' 'Siple Dome' 'South Pole' 'Taylor Dome' 'Vostok' 'WAIS Divide'};
name_core_long2             = [name_core_long 'Byrd' 'Palmer Land'];
name_core_short             = {'UFG' 'G' ' MG' 'NG' 'SD' 'SP' 'TD' 'V' 'WD' 'B' 'PL'};

color_core_alt              = [0    0    0;
                               0    0    1;
                               0    0.5  1;
                               0    1    0;
                               0.5  1    0;
                               0.9  0.9  0;
                               1    0.5  0;
                               1    0    0;
                               1    0    0.5;
                               1    0    1;
                               0.5  0    1];

color_core                  = [0   0    0;
                               1   0    0;
                               1   0.5  0;
                               1   0.75 0;
                               1   0    1;
                               0   0.75 0;
                               0   1    0.75;
                               0.5 1    1
                               0   0    1];
color_core2                 = [color_core;
                               0.5 0.5  0.5;
                               1   1    0];

color_depth                 = [0    0    0;
                               1    0    0;
                               0    0    1;
                               0    0.5  0;
                               1    0    1;
                               0    1    1;
                               1    0.5  0;
                               0.5  0    1;
                               1    1    0;
                               0    1    0.5];

% trimming data based on frequency
freq_trim                    = {{-Inf 4e5; -Inf 4e5; -Inf 1e5; -Inf 1e5};
                                {3e-1 6e5; 3e-1 6e5; 3e-1 6e5; 3e-1 6e5; 3e-1 6e5; 3e-1 6e5};
                                {-Inf 4e5; -Inf 4e5; -Inf Inf; -Inf Inf; 3e-1 4e5; -Inf 4e5; -Inf 4e5; -Inf 4e5; 2e1 4e5; -Inf 4e5};
                                {-Inf 2e5; -Inf 2e5; -Inf 2e5; 2e0 2e5; -Inf 2e5};
                                {3e-1 1e5; 1e-1 1e6; 3e-1 4e5; 3e-1 1e5; 3e-1 1e5};
                                {2e0 5e5; 1e1 5e5; 1e0 5e5; -Inf 5e5};
                                {3e-1 6e5; -Inf 6e5; 3e-1 6e5};
                                {-Inf 1e6; 3e-1 4e5; 3e-1 2e5; -Inf 1e6; -Inf 8e5; -Inf 1e6; -Inf 1e5; 5e-2 1e5};
                                {-Inf 1e5; 3e-1 1e5; -Inf 1e5; 3e-1 1e5}};
for ii = 1:num_core
    for jj = 1:core{ii}.num_files
        ind_curr            = find((core{ii}.freq{jj} >= freq_trim{ii}{jj, 1}) & (core{ii}.freq{jj} <= freq_trim{ii}{jj, 2}));
        [core{ii}.permitt_real{jj}, core{ii}.permitt_imag{jj}, core{ii}.conduct{jj}, core{ii}.freq{jj}] ...
                            = deal(core{ii}.permitt_real{jj}(ind_curr), core{ii}.permitt_imag{jj}(ind_curr), core{ii}.conduct{jj}(ind_curr), core{ii}.freq{jj}(ind_curr));
    end
end

% plot parameters
letters                     = 'a':'z';
letters_cell                = cell(1, 27);
for ii = 1:26
    letters_cell{ii}        = letters(ii);
end
letters_cell{27}            = 'aa';
letters_cell{28}            = 'bb';
letters_cell{29}            = 'cc';
letters_cell{30}            = 'dd';
plots                       = {'permitt_real' 'permitt_imag', 'conduct'};
ylabels                     = {{'Real part'; 'of permittivity'} {'Imaginary part'; 'of permittivity'} 'Conductivity (S m^{-1})'};
data_range                  = [1e0 1e3; 1e-1 1e3; 1e-10 1e-4];
num_ytick                   = [4 5 7];
freq_range                  = [1e-2 1e6];
subplot_start               = [0.14 0.56 0.555 0.41];
y_letter                    = [550 4e2 3e-5];
yticks                      = {{'' '10^1' '10^2' '10^3'} {'' '10^0' '10^1' '10^2'} {'10^{-10}' '' '10^{-8}' '' '10^{-6}' '' ''}};
ind_disp                    = [1 3 4 2 5 9 7 6 8];
ind_disp_alt                = [1 11 3 5 4 10 9 2 7 6 8];
[~, ind_disp2]              = sort(ind_disp_alt);
ind_disp2                   = 12 - ind_disp2;
color_core_alt              = color_core_alt(ind_disp2, :);
ind_text                    = [];

% highlight samples (core number / depth number)
ind_highlight               = flipud([1 4; 3 6; 9 4; 7 2; 8 8]);

% load Bob's two relaxation analysis
fid                         = fopen('data/Rlxn_Correl_030215.txt', 'r');
tmp                         = textscan(fid, '%s%f%f%f%f%f%f%f%f%f%f%f%s', 'headerlines', 1);
fclose(fid);

% relaxation structure
rlx                         = struct;

[rlx.core, rlx.depth, rlx.Cl_1, rlx.Cl_2, rlx.conduct_hf, rlx.conduct_1, rlx.conduct_2, rlx.H_tot, rlx.Cl_tot, rlx.NH4_tot, rlx.temp, rlx.firn_depth, rlx.notes] ...
                            = deal(tmp{1}, tmp{2}, tmp{3}, tmp{4}, tmp{5}, tmp{6}, tmp{7}, tmp{8}, tmp{9}, tmp{10}, tmp{11}, tmp{12}, tmp{13});


% get rid of some bad points
ind_keep                    = find(~strcmp(rlx.notes, 'Remove') & ~strcmp(rlx.core, 'Kawada'));
[rlx.core, rlx.depth, rlx.Cl_1, rlx.Cl_2, rlx.conduct_hf, rlx.conduct_1, rlx.conduct_2, rlx.H_tot, rlx.Cl_tot, rlx.NH4_tot, rlx.temp, rlx.firn_depth, rlx.notes] ...
                            = deal(rlx.core(ind_keep), rlx.depth(ind_keep), rlx.Cl_1(ind_keep), rlx.Cl_2(ind_keep), rlx.conduct_hf(ind_keep), rlx.conduct_1(ind_keep), rlx.conduct_2(ind_keep), rlx.H_tot(ind_keep), rlx.Cl_tot(ind_keep), rlx.NH4_tot(ind_keep), rlx.temp(ind_keep), ...
                                   rlx.firn_depth(ind_keep), rlx.notes(ind_keep));
num_rlx                     = length(rlx.core);

name_core_convert           = {'Vostok' 'SouthPole' 'Taylor' 'GISP2' 'WAIS' 'Newall' 'SipleDome' 'Mullins' 'Fremont' 'Byrd' 'PalmerLand'};
ind_core_convert            = [8 6 7 2 9 4 5 3 1 10 11];
rlx.ind_core                = NaN(num_rlx, 1);
for ii = 1:num_rlx
    rlx.ind_core(ii)        = ind_core_convert(strcmp(rlx.core{ii}, name_core_convert));
end
%%
% sort from cold to warm with increasing depth
[~, ind_sort]               = sortrows([ind_disp2(rlx.ind_core)' rlx.depth]);
[rlx.core, rlx.depth, rlx.Cl_1, rlx.Cl_2, rlx.conduct_hf, rlx.conduct_1, rlx.conduct_2, rlx.H_tot, rlx.Cl_tot, rlx.NH4_tot, rlx.temp, rlx.firn_depth, rlx.notes, rlx.ind_core] ...
                            = deal(rlx.core(ind_sort), rlx.depth(ind_sort), rlx.Cl_1(ind_sort), rlx.Cl_2(ind_sort), rlx.conduct_hf(ind_sort), rlx.conduct_1(ind_sort), rlx.conduct_2(ind_sort), rlx.H_tot(ind_sort), rlx.Cl_tot(ind_sort), rlx.NH4_tot(ind_sort), rlx.temp(ind_sort), ...
                                   rlx.firn_depth(ind_sort), rlx.notes(ind_sort), rlx.ind_core(ind_sort));
[rlx.core(37:44), rlx.depth(37:44), rlx.Cl_1(37:44), rlx.Cl_2(37:44), rlx.conduct_hf(37:44), rlx.conduct_1(37:44), rlx.conduct_2(37:44), rlx.H_tot(37:44), rlx.Cl_tot(37:44), rlx.NH4_tot(37:44), rlx.temp(37:44), rlx.firn_depth(37:44), rlx.notes(37:44), rlx.ind_core(37:44)] ...
                            = deal(rlx.core([40 38 37 39 41:44]), rlx.depth([40 38 37 39 41:44]), rlx.Cl_1([40 38 37 39 41:44]), rlx.Cl_2([40 38 37 39 41:44]), rlx.conduct_hf([40 38 37 39 41:44]), rlx.conduct_1([40 38 37 39 41:44]), rlx.conduct_2([40 38 37 39 41:44]), ...
                                   rlx.H_tot([40 38 37 39 41:44]), rlx.Cl_tot([40 38 37 39 41:44]), rlx.NH4_tot([40 38 37 39 41:44]), rlx.temp([40 38 37 39 41:44]), rlx.firn_depth([40 38 37 39 41:44]), rlx.notes([40 38 37 39 41:44]), rlx.ind_core([40 38 37 39 41:44]));

% fix Mullins 15.34 m
[rlx.Cl_1(41), rlx.Cl_2(41)]= deal(129.04, 33.35);
                               
[rlx.H_1, rlx.H_2]          = deal((rlx.Cl_1 .* 0.48947), (rlx.Cl_2 .* 0.48947));
[rlx.conduct_1(rlx.conduct_1 <= 0), rlx.conduct_2(rlx.conduct_2 <= 0), rlx.H_1(rlx.H_1 <= 0), rlx.H_2(rlx.H_2 <= 0), rlx.Cl_1(rlx.Cl_1 <= 0), rlx.Cl_2(rlx.Cl_2 <= 0)] ...
                            = deal(NaN);

[~, ind_sort]               = sort(rlx.temp);
[rlx.core, rlx.depth, rlx.Cl_1, rlx.Cl_2, rlx.conduct_hf, rlx.conduct_1, rlx.conduct_2, rlx.H_tot, rlx.Cl_tot, rlx.NH4_tot, rlx.temp, rlx.firn_depth, rlx.notes, rlx.ind_core] ...
                            = deal(rlx.core(ind_sort), rlx.depth(ind_sort), rlx.Cl_1(ind_sort), rlx.Cl_2(ind_sort), rlx.conduct_hf(ind_sort), rlx.conduct_1(ind_sort), rlx.conduct_2(ind_sort), rlx.H_tot(ind_sort), rlx.Cl_tot(ind_sort), rlx.NH4_tot(ind_sort), rlx.temp(ind_sort), ...
                                   rlx.firn_depth(ind_sort), rlx.notes(ind_sort), rlx.ind_core(ind_sort));

                               
% temp_ord                    = [8 6 7 2 9 10 4 5 3 11 1]; % order of increasing surface temperature
% 
% rlx_ord                     = struct;
% jj                          = 0;
% for ii = 1:length(temp_ord)
%     rlx_ord(jj + (1:length(find(rlx.ind_core == temp_ord(ii))))) ...
%                             = rlx(logical(rlx.ind_core == temp_ord(ii)));
%     jj                      = jj + length(find(rlx.ind_core == temp_ord(ii)));
% end
%%
% rlx                         = rlx_ord;

% indices of various things
ind_firn                    = find(strcmp(rlx.notes, 'firn'));
ind_relax_1                 = find(isnan(rlx.Cl_2));
ind_relax_2                 = find(~isnan(rlx.Cl_2));
ind_relax_3                 = find(isnan(rlx.conduct_2));
ind_relax_4                 = find(~isnan(rlx.conduct_2));
ind_acc                     = 4:8; % accreted ice

name_core_fig2              = dir('data/signatures_fig2/');
name_core_fig2              = {name_core_fig2(4:end).name};

num_core_fig2               = length(name_core_fig2);

core_fig2                   = cell(1, num_core_fig2);
ind_40                      = NaN(1, num_core_fig2);
for ii = 1:num_core_fig2
    core_fig2{ii}           = load(['data/signatures_fig2/' name_core_fig2{ii}]);
    ind_40(ii)              = interp1(core_fig2{ii}.temp_mean, 1:core_fig2{ii}.num_temp, -40, 'nearest', 'extrap');
end

ind_highlight_match         = [4 3 5 2 1];

relax_description           = {'accreted from subglacial lake' 'meteoric polar ice sheet' 'meteoric polar ice sheet' 'ancient meteoric polar glacier' 'temperate glacier'}; % figure 2

% figure 3
temp_vec                    = -90:10:0;
temp_inv_vec                = 1e3 ./ (temp_vec + 273.15);

[activ_energy, activ_poly, vars_all, temp_inv, freq_relax, freq_relax_min, freq_relax_max] ...
                            = deal(cell(1, 2));

kawada                      = load('data/arrhenius/Kawada_1978');
vars_all{1}.temp_mean       = kawada.temp;
temp_inv{1}                 = 1e3 ./ (kawada.temp + 273.15);
freq_relax{1}               = kawada.freq_relax';

files                       = {'NaCl_0_85M_041408_merge_inv1'};
num_files                   = length(files);

num_relax                   = [1 1];
temp_break                  = [-51 NaN NaN;
                               NaN NaN NaN];

% load standardized files
for ii = 2
    vars_all{ii}            = load(['data/two_relax_fig/' files{ii - 1}]); % all the variables in the loaded data file
    temp_inv{ii}            = 1e3 ./ (vars_all{ii}.temp_mean + 273.15); % inverted temperature
    freq_relax{ii}          = vars_all{ii}.freq_relax_cat(1:num_relax(ii), :);
    [freq_relax_min{ii}, freq_relax_max{ii}] ...
                            = deal(cell(1, num_relax(ii)));
    for jj = 1:num_relax(ii)
        [freq_relax_min{ii}{jj}, freq_relax_max{ii}{jj}] ...
                            = deal(vars_all{ii}.freq_relax_std_cat{jj}(:, 1), vars_all{ii}.freq_relax_std_cat{jj}(:, 2));
    end
end


for ii = 1:2
    [activ_poly{ii}, activ_energy{ii}] ...
                            = deal(cell(2, num_relax(ii)));
    for jj = 1:num_relax(ii)
        tmp                 = log(freq_relax{ii}(jj, ~isnan(freq_relax{ii}(jj, :))));
        if isrow(tmp)
            tmp             = tmp';
        end
        if isnan(temp_break(ii, jj))
            activ_poly{ii}{1, jj} ...
                            = polyfit(temp_inv{ii}(~isnan(freq_relax{ii}(jj, :))), tmp, 1); % polynomials for only 1 activation energy
            activ_energy{ii}{jj} ...
                            = (-1e3 * boltzmann / conv_eV_J) .* activ_poly{ii}{1, jj}(1); % in eV
        else
            ind1            = find(vars_all{ii}.temp_mean(~isnan(freq_relax{ii}(jj, :))) < temp_break(ii, jj)); % low temp indices
            ind2            = find(vars_all{ii}.temp_mean(~isnan(freq_relax{ii}(jj, :))) > temp_break(ii, jj)); % high temp indices
            activ_poly{ii}{1, jj} ...
                            = polyfit(temp_inv{ii}(ind1), tmp(ind1), 1);
            activ_poly{ii}{2, jj} ...
                            = polyfit(temp_inv{ii}(ind2), tmp(ind2), 1);
            activ_energy{ii}{1, jj} ...
                            = (-1e3 * boltzmann / conv_eV_J) .* activ_poly{ii}{1, jj}(1);
            activ_energy{ii}{2, jj} ...
                            = (-1e3 * boltzmann / conv_eV_J) .* activ_poly{ii}{2, jj}(1);
        end
    end
end

%                              name in legend                           marker  color            temp range (C)                  activ energy x/y
plot_data                   = {'Pure'                                   '>'     [0.5 0.5 0.5]    [-130    (225 - 273.15)]        1   1   1;
                               'Pure'                                   '>'     [0.5 0.5 0.5]    [(225 - 273.15) 0]              1   1   2;
                               'Cl^--saturated'                         'd'     'm'              [vars_all{2}.temp_mean(1) -10]  2   1   1};

name_core_fig3              = dir('data/signatures_fig3/');
name_core_fig3              = {name_core_fig3(4:end).name};

num_core_fig3               = length(name_core_fig3);

core_fig3                   = cell(1, num_core_fig3);
for ii = 1:num_core_fig2
    core_fig3{ii}           = load(['data/signatures_fig3/' name_core_fig3{ii}]);
end
ind_core_fig3               = [1 3 7 8 9];

%%
if plotting

%% docked

    set(0, 'DefaultFigureWindowStyle', 'docked')
    
%% floating

    set(0, 'DefaultFigureWindowStyle', 'default')

%% ALL ON ONE FIGURE (OOF)
    figure('position', [200 200 600 1000])
    for ii = 1:3
        subplot('position', [0.15 (0.68 - (0.3 * (ii - 1))) 0.80 0.3])
        axis([freq_range data_range(ii, :)])
        hold on
        for jj = 1:num_core
            for kk = 1:core{jj}.num_files
                loglog(core{jj}.freq{kk}, eval(['core{jj}. ' plots{ii} '{kk}']), 'linewidth', 3, 'color', color_core(jj, :))
            end
        end
        set(gca, 'fontsize', 20, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'yscale', 'log', 'ytick', logspace(log10(data_range(ii, 1)), log10(data_range(ii, 2)), num_ytick(ii)))
        text(0.02, y_letter(ii), letters(ii), 'fontsize', 22, 'color', 'k', 'fontweight', 'bold')
        ylabel(ylabels{ii})
        box on
        grid on
        switch ii
            case 1
                set(gca, 'xticklabel', {})
                text(3e4, 550, '-40{\circ}C', 'fontsize', 22, 'color', 'k', 'fontweight', 'bold')
            case 2
                set(gca, 'xticklabel', {})
            case 3
                xlabel('Frequency (Hz)')
        end
    end
    
%% FIGURE 1 A-O SPECTRA DUMP

    figure('position', [50 50 1400 850])
    ll                      = 1;
    ax                      = NaN(3, 5);
    for ii = fliplr(ind_disp(5:end))
        for jj = 1:3
            ax(jj, ll)      = subplot('position', [(0.08 + (0.18 * (ll - 1))) (0.665 - ((0.18 * (1400 / 850)) * (jj - 1))) 0.18 (0.18 * (1400 / 850))]);
            axis([freq_range data_range(jj, :)])
            hold on
            name_depth      = cell(1, core{ii}.num_files);
            for kk = 1:core{ii}.num_files
                loglog(core{ii}.freq{kk}, eval(['core{ii}. ' plots{jj} '{kk}']), 'linewidth', 3, 'color', color_depth(kk, :))
                name_depth{kk} ...
                            = sprintf('%5.1f', core{ii}.depth(kk));
            end
            set(gca, 'fontsize', 22, 'xscale', 'log', 'xtick', logspace(-1, 5, 7), 'xticklabel', {'' '10^0' '' '10^2' '' '10^4' ''}, 'yscale', 'log', 'ytick', logspace(log10(data_range(jj, 1)), log10(data_range(jj, 2)), num_ytick(jj)), 'yticklabel', yticks{jj})
            text(0.02, y_letter(jj), letters_cell{jj + (3 * (ll - 1))}, 'fontsize', 24, 'color', 'k', 'fontweight', 'bold')%, 'edgecolor', 'k', 'backgroundcolor', 'w')
            if (ll == 1)
                ylabel(ylabels{jj})
            else
                set(gca, 'yticklabel', {})
            end
            box on
            grid on
            switch jj
                case 1
                    title(name_core_long{ii}, 'fontsize', 22, 'fontweight', 'bold', 'color', 'k')
                    set(gca, 'xticklabel', {})
                    legend(name_depth, 'fontsize', 20, 'location', 'southwest')
                case 2
                    set(gca, 'xticklabel', {})
                case 3
                    if (ii == 6)
                        annotation('arrow', [0.427 0.427], [0.24 0.29], 'color', 'k', 'linewidth', 2, 'headstyle', 'plain')
                        text(4e4, 1.5e-7, '\sigma_{HF}', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k')
                    end
                    if (ll == 3)
                        xlabel('Frequency (Hz)')
                    end
            end
        end
        linkaxes(ax(:, ll), 'x')
        ll                  = ll + 1;
    end
    linkaxes(ax(1, :), 'y')
    linkaxes(ax(2, :), 'y')
    linkaxes(ax(3, :), 'y')
    
%% FIGURE 1 P-AA SPECTRA DUMP

    figure('position', [50 50 1400 850])
    ax                      = NaN(3, 4);
    ll                      = 1;
    for ii = fliplr(ind_disp(1:4))
        for jj = 1:3
            ax(jj, ll)      = subplot('position', [(0.08 + (0.18 * (ll - 1))) (0.665 - ((0.18 * (1400 / 850)) * (jj - 1))) 0.18 (0.18 * (1400 / 850))]);
            axis([freq_range data_range(jj, :)])
            hold on
            name_depth      = cell(1, core{ii}.num_files);
            for kk = 1:core{ii}.num_files
                loglog(core{ii}.freq{kk}, eval(['core{ii}. ' plots{jj} '{kk}']), 'linewidth', 3, 'color', color_depth(kk, :))
                name_depth{kk} ...
                            = sprintf('%5.1f', core{ii}.depth(kk));
            end
            set(gca, 'fontsize', 22, 'xscale', 'log', 'xtick', logspace(-1, 5, 7), 'xticklabel', {'' '10^0' '' '10^2' '' '10^4' ''}, 'yscale', 'log', 'ytick', logspace(log10(data_range(jj, 1)), log10(data_range(jj, 2)), num_ytick(jj)), 'yticklabel', yticks{jj})
            text(0.02, y_letter(jj), letters_cell{jj + (3 * (ll + 5 - 1))}, 'fontsize', 24, 'color', 'k', 'fontweight', 'bold')%, 'edgecolor', 'k', 'backgroundcolor', 'w')
            if (ll == 1)
                ylabel(ylabels{jj})
            else
                set(gca, 'yticklabel', {})
            end
            box on
            grid on
            switch jj
                case 1
                    title(name_core_long{ii}, 'fontsize', 22, 'fontweight', 'bold', 'color', 'k')
                    set(gca, 'xticklabel', {})
                    lg      = legend(name_depth, 'location', 'southwest');
                    set(lg, 'fontsize', 20)
                case 2
                    set(gca, 'xticklabel', {})
                case 3
                    if (ll == 3)
                        xlabel('Frequency (Hz)')
                    end
            end
        end
        linkaxes(ax(:, ll), 'x')
        ll                  = ll + 1;
    end
    mm                      = 1;
    for ii = (num_core + 1):(num_core + 2)
        for jj = 1:3
            ax(jj, ll)          = subplot('position', [(0.08 + (0.18 * (ll - 1))) (0.665 - ((0.18 * (1400 / 850)) * (jj - 1))) 0.18 (0.18 * (1400 / 850))]);
            axis([freq_range data_range(jj, :)])
            hold on
            loglog(core{ii}.freq, eval(['core{ii}. ' plots{jj}]), 'linewidth', 3, 'color', color_depth(mm, :))
            set(gca, 'fontsize', 22, 'xscale', 'log', 'xtick', logspace(-1, 5, 7), 'xticklabel', {'' '10^0' '' '10^2' '' '10^4' ''}, 'yscale', 'log', 'ytick', logspace(log10(data_range(jj, 1)), log10(data_range(jj, 2)), num_ytick(jj)), 'yticklabel', yticks{jj})
            set(gca, 'yticklabel', {})
            text(0.02, y_letter(jj), letters_cell{jj + (3 * (ll + 5 - 1))}, 'fontsize', 24, 'color', 'k', 'fontweight', 'bold')%, 'edgecolor', 'k', 'backgroundcolor', 'w')
            box on
            grid on
            switch jj
                case 1
                    title('Other published spectra', 'fontsize', 22, 'fontweight', 'bold', 'color', 'k')
                    set(gca, 'xticklabel', {})
                    lg      = legend('Byrd', 'Palmer Land');
                    set(lg, 'fontsize', 20, 'location', 'southwest')
                case 2
                    set(gca, 'xticklabel', {})
            end
        end
        mm                   = mm + 1;
    end
    linkaxes(ax(1, :), 'y')
    linkaxes(ax(2, :), 'y')
    linkaxes(ax(3, :), 'y')
    
%% FIGURE 2 FITTED SPECTRA OF HIGHLIGHTED SAMPLES

    name_tmp                = cell(1, size(ind_highlight, 1));
    pl                      = NaN(1, size(ind_highlight, 1));
    color_relax             = {'r' 'b' 'm'};
    figure('position', [2800 200 1800 1000])
    for ii = 1:size(ind_highlight, 1)
        ind_curr            = ind_highlight_match(ii);
        ind_40_curr         = ind_40(ind_curr);
        num_relax_curr      = core_fig2{ind_curr}.num_relax(ind_40_curr);
        freq_curr           = core_fig2{ind_curr}.freq_trim{ind_40_curr};
        if core_fig2{ind_curr}.do_dc(ind_40_curr)
            model           = NaN((2 + (3 * num_relax_curr)), 1);
            model(end)      = -log10(core_fig2{ind_curr}.conduct_dc(ind_40_curr));
        else
            model           = NaN((1 + (3 * num_relax_curr)), 1);
        end
        try
            model(1)        = core_fig2{ind_curr}.permitt_hf_model(ind_40_curr);
        catch
            model(1)        = core_fig2{ind_curr}.permitt_hf(ind_40_curr);
        end
        for jj = 1:num_relax_curr
            [model(2 + ((jj - 1) * 3)), model(3 + ((jj - 1) * 3)), model(4 + ((jj - 1) * 3))] ...
                            = deal(core_fig2{ind_curr}.permitt_diff_cat(jj, ind_40_curr), (1 / (2 * pi * core_fig2{ind_curr}.freq_relax_cat(jj, ind_40_curr))), core_fig2{ind_curr}.alpha_cole_cat(jj, ind_40_curr));
        end
        permitt_complex_full= fm_cole(model, freq_curr, num_relax_curr, core_fig2{ind_curr}.do_dc(ind_40_curr), permitt_vacuum);
        [permitt_real_full, permitt_imag_full] ...
                            = deal(permitt_complex_full(1:length(freq_curr)), permitt_complex_full((length(freq_curr) + 1):end));
        conduct_full        = (2 * pi * permitt_vacuum) .* freq_curr .* permitt_imag_full;
        permitt_complex     = fm_cole_sep(model, freq_curr, num_relax_curr, core_fig2{ind_curr}.do_dc(ind_40_curr), permitt_vacuum);
        [permitt_real, permitt_imag, conduct] ...
                            = deal(cell(1, num_relax_curr));
        tmp                 = permitt_complex{1}(1);
        for jj = 1:num_relax_curr
            [permitt_real{jj}, permitt_imag{jj}] ...
                            = deal((permitt_complex{1 + jj}(1:length(freq_curr)) + tmp), permitt_complex{1 + jj}((length(freq_curr) + 1):end));
            conduct{jj}     = (2 * pi * permitt_vacuum) .* freq_curr .* permitt_imag{jj};
            tmp             = tmp + model(double(core_fig2{ind_curr}.do_dc(ind_40_curr)) + 1 + ((jj - 1) * 3));
        end
        for jj = 1:3
            subplot('position', [(0.06 + ((ii - 1) * 0.185)) (0.645 - (0.29 * (jj - 1))) 0.185 0.29])
            hold on
            set(gca, 'fontsize', 22, 'xscale', 'log', 'xtick', logspace(-2, 6, 9), 'yscale', 'log')
            axis([freq_range data_range(jj, :)])
            set(gca, 'ytick', logspace(log10(data_range(jj, 1)), log10(data_range(jj, 2)), num_ytick(jj)))
            plg(1)          = loglog(core{ind_highlight(ii, 1)}.freq{ind_highlight(ii, 2)}, eval(['core{ind_highlight(ii, 1)}. ' plots{jj} '{ind_highlight(ii, 2)}']), 'ko', 'linewidth', 1, 'markersize', 8, 'markerfacecolor', [0.75 0.75 0.75]);
            plg(2)          = loglog(freq_curr, eval([plots{jj} '_full']), 'k', 'linewidth', 3);
            for kk = 1:num_relax_curr
                if ((jj == 1) && (kk > 1))
                    plg(kk + 2) ...
                            = loglog(freq_curr(find((diff([permitt_real{kk} permitt_complex_full(1:length(freq_curr))], 1, 2) > -0.1), 1):end), ...
                                     eval([plots{jj} '{kk}((find((diff([permitt_real{kk} permitt_complex_full(1:length(freq_curr))], 1, 2) > -0.1), 1)):end)']), '--', 'color', color_relax{kk}, 'linewidth', 3);
                else
                    plg(kk + 2) ...
                            = loglog(freq_curr, eval([plots{jj} '{kk}']), '--', 'color', color_relax{kk}, 'linewidth', 3);
                end
                if ((kk == 3) || (core_fig2{ind_curr}.freq_relax_cat(kk, ind_40_curr) < 1e1))
                    set(plg(kk + 2), 'color', [0.75 0.75 0.75])
                end
            end
            if core_fig2{ind_curr}.do_dc(ind_40_curr)
                switch jj
                    case 2
                        loglog(freq_curr, (core_fig2{ind_curr}.conduct_dc(ind_40_curr(ones(length(freq_curr), 1))) ./ ((2 * pi * permitt_vacuum) .* freq_curr)), '--', 'color', [0.75 0.75 0.75], 'linewidth', 3)
                    case 3
                        loglog(freq_curr, core_fig2{ind_curr}.conduct_dc(ind_40_curr(ones(length(freq_curr), 1))), '--', 'color', [0.75 0.75 0.75], 'linewidth', 3)
                end
            end
            box on
            grid on
            text(0.02, y_letter(jj), letters(jj + ((ii - 1) * 3)), 'fontsize', 24, 'color', 'k', 'fontweight', 'bold')
            if (ii > 1)
                set(gca, 'yticklabel', {}, 'xtick', logspace(-2, 6, 9), 'xticklabel', {'' '' '10^0' '' '10^2' '' '10^4' '' '10^6'})
            else
                set(gca, 'xtick', logspace(-2, 6, 9), 'xticklabel', {'10^{-2}' '' '10^0' '' '10^2' '' '10^4' '' '10^6'})                
                ylabel(ylabels{jj})
            end
            switch jj
                case 1
                    set(gca, 'xticklabel', {})
                    title({[name_core_long{ind_highlight(ii, 1)} ' / ' sprintf('%4.1f', core{ind_highlight(ii, 1)}.depth(ind_highlight(ii, 2))) ' m']; relax_description{ii}}, 'fontweight', 'bold')
                case 2
                    set(gca, 'xticklabel', {})
                    if (ii == 1)
                        set(gca, 'ytick', logspace(-1, 3, 5), 'yticklabel', {'10^{-1}' '10^0' '10^1' '10^2' ''})
                    end
                    if (ii == size(ind_highlight, 1))
                        plg(4) = plot(NaN, NaN, 'b--', 'linewidth', 3);
                        lg  = legend(plg, 'data', 'full Cole-Cole fit', 'fast ice relaxation', 'slow ice relaxation', 'other', 'location', 'northeast');
                        set(lg, 'fontsize', 20)
                    end
                case 3
                    switch ii
                        case 1
                            set(gca, 'ytick', logspace(-10, -4, 7), 'yticklabel', {'10^{-10}' '10^{-9}' '10^{-8}' '10^{-7}' '10^{-6}' '10^{-5}' ''})
                        case 3
                            xlabel('Frequency (Hz)')
                    end
            end
        end
    end
    
%% FIGURE 3: TEMPERATURE DEPENDENCE OF HF CONDUCTIVITY AND RELAXATION FREQUENCIES

    color_relax             = {'k' 'r' 'b'};
    marker_relax            = {'o' '^' 'v'};
    figure('position', [2600 200 800 800])
    subplot('position', [0.10 0.07 0.8 0.85])
    axis([(1e3 / 273.15) (1e3 / (273.15 - 90)) 2e0 3e4])
    hold on
    prfo                    = NaN(1, size(plot_data, 1));
    for ii = 1:size(plot_data, 1)
        prfo(ii)            = plot((1e3 ./ (plot_data{ii, 4} + 273.15)), exp(polyval(activ_poly{plot_data{ii, 5}}{plot_data{ii, 7}, plot_data{ii, 6}}, (1e3 ./ (plot_data{ii, 4} + 273.15)))), '--', 'linewidth', 3, 'color', plot_data{ii, 3});
    end
    prf                     = NaN(1, 3);
    prf(1)                  = semilogy(NaN, NaN, 'ko', 'markersize', 10, 'linewidth', 1, 'markerfacecolor', [0.75 0.75 0.75]);
    prf(2)                  = semilogy(NaN, NaN, 'k^', 'markersize', 10, 'linewidth', 1, 'markerfacecolor', [0.75 0.75 0.75]);
    prf(3)                  = semilogy(NaN, NaN, 'kv', 'markersize', 10, 'linewidth', 1, 'markerfacecolor', [0.75 0.75 0.75]);
    for ii = fliplr(ind_highlight_match)
        if (all(core_fig3{ii}.num_relax == 1) || (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1))
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), (core_fig3{ii}.freq_relax_cat(1, :) - core_fig3{ii}.freq_relax_std_cat{1}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{1}(:, 2)' - core_fig3{ii}.freq_relax_cat(1, :)), ...
                     'color', color_core_alt(ind_core_fig3(ii), :), 'linestyle', 'none')
        else
            for jj = 1:core_fig3{ii}.num_temp
                semilogy(repmat((1e3 ./ (273.15 + core_fig3{ii}.temp_mean(jj))), 1, 2), core_fig3{ii}.freq_relax_cat(1:2, jj), '--', 'color', color_core_alt(ind_core_fig3(ii), :), 'linewidth', 2)
            end
            semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(2, :), 'k', 'marker', marker_relax{3}, 'markerfacecolor', color_core_alt(ind_core_fig3(ii), :), 'markersize', 10, 'linewidth', 1, 'linestyle', 'none')            
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), (core_fig3{ii}.freq_relax_cat(1, :) - core_fig3{ii}.freq_relax_std_cat{1}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{1}(:, 2)' - core_fig3{ii}.freq_relax_cat(1, :)), color_relax{2}, ...
                     'color', color_core_alt(ind_core_fig3(ii), :), 'linestyle', 'none')
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(2, :), (core_fig3{ii}.freq_relax_cat(2, :) - core_fig3{ii}.freq_relax_std_cat{2}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{2}(:, 2)' - core_fig3{ii}.freq_relax_cat(2, :)), color_relax{3}, ...
                     'color', color_core_alt(ind_core_fig3(ii), :), 'linestyle', 'none')
        end
%         if (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1)
%             for jj = 2:core_fig3{ii}.num_relax
%                 semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(jj, :), 'k', 'marker', 'p', 'markerfacecolor', [0.75 0.75 0.75], 'markersize', 10, 'linewidth', 1, 'linestyle', 'none');
%             end
%         end
    end
    for ii = fliplr(ind_highlight_match)
        if (all(core_fig3{ii}.num_relax == 1) || (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1))
            semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), 'k', 'marker', marker_relax{1}, 'markerfacecolor', color_core_alt(ind_core_fig3(ii), :), 'markersize', 10, 'linewidth', 1, 'linestyle', 'none')
        else
            semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), 'k', 'marker', marker_relax{2}, 'markerfacecolor', color_core_alt(ind_core_fig3(ii), :), 'markersize', 10, 'linewidth', 1, 'linestyle', 'none')
            semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(2, :), 'k', 'marker', marker_relax{3}, 'markerfacecolor', color_core_alt(ind_core_fig3(ii), :), 'markersize', 10, 'linewidth', 1, 'linestyle', 'none')
        end
%         if (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1)
%             for jj = 2:core_fig3{ii}.num_relax
%                 semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(jj, :), 'k', 'marker', 'p', 'markerfacecolor', [0.75 0.75 0.75], 'markersize', 10, 'linewidth', 1, 'linestyle', 'none');
%             end
%         end
    end    
    tmp_log                 = logspace(1.85, 1.25, 5);
    for ii = 1:num_core_fig3
        text(3.7, tmp_log(ii), [name_core_long{ind_core_fig3(ind_highlight_match(ii))} ' / ' sprintf('%4.1f', core{ind_highlight(ii, 1)}.depth(ind_highlight(ii, 2))) ' m'], 'color', color_core_alt(ind_core_fig3(ind_highlight_match(ii)), :), 'fontsize', 20, 'fontweight', 'bold')
    end
    set(gca, 'fontsize', 20, 'yscale', 'log', 'yticklabel', {'10^0' '10^1' '10^2' '10^3' '10^4' ''})
    xlabel('1000 / Temperature (K)')
    ylabel('Relaxation frequency (Hz)')
    lg                      = legend([prfo([2 3]) prf], {plot_data{[2 3], 1} 'Single ice relaxation' 'Fast ice relaxation' 'Slow ice relaxation'}, 'location', 'southwest');
    set(lg, 'fontsize', 20)
    box off
    ax                      = gca;
    axes('position', get(gca, 'position'), 'color', 'none', 'fontsize', 20, 'xaxislocation', 'top', 'yaxislocation', 'right', 'xlim', get(gca, 'xlim'), 'yscale', 'log', 'ydir', 'reverse')
    xlabel('Temperature (\circC)')
    ylabel('Relaxation time (s)', 'rotation', 270, 'position', [5.66 5e-4 1.00005])
    ylim(fliplr(1 ./ ((2 * pi) .* get(ax, 'ylim'))))
    set(gca, 'xtick', fliplr(temp_inv_vec), 'xticklabel', fliplr(temp_vec))
    
%% FIGURE 4: TWO-RELAXATION COMPARISON

    figure('position', [2600 200 1600 825], 'color', 'w')
    subplot('position', [0.055 0.52 0.89 0.44])
    hold on
    axis([0 (length(rlx.depth) + 1) 1e-7 4e-5])
    tmp                     = rlx.core{1};
    for ii = 1:num_rlx
        if ~strcmp(rlx.core{ii}, tmp)
            plot(repmat((ii - 0.5), 1, 2), get(gca, 'ylim'), 'k', 'linewidth', 1)
        end
        tmp                 = rlx.core{ii};
    end
    for ii = ind_relax_3'
        p_tmp               = plot(ii, rlx.conduct_1(ii), 'ko', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
        if any(ii == ind_acc)
            set(p_tmp, 'marker', 'd')
        end
        if any(ii == ind_relax_2)
            set(p_tmp, 'marker', '^')
        end
    end
    for ii = ind_relax_4'
        plot([ii ii], [rlx.conduct_1(ii), rlx.conduct_2(ii)], 'k--', 'linewidth', 3)
        p_tmp               = plot(ii, rlx.conduct_1(ii), 'k^', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
        p_tmp               = plot(ii, rlx.conduct_2(ii), 'kv', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
    end
    pdu2                    = NaN(1, num_core2);
    for ii = 1:num_core2
        pdu2(ii)            = plot(NaN, NaN, 'ko', 'markerfacecolor', color_core_alt(ii, :), 'markersize', 12, 'linewidth', 1.5);
%         text(6e-5, name_core_short{ii}, 'fontweight', 'bold', 'fontsize', 20, 'color', color_core_alt(ii, :))
    end
    set(gca, 'fontsize', 20, 'xtick', [], 'yscale', 'log', 'ytick', logspace(-7, -4, 4), 'ygrid', 'on')
    text(0.5, 2.5e-5, 'a', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold')
    ylabel('HF conductivity at -40\circC (S m^{-1})')
    tmp                     = {pdu2(fliplr(ind_disp_alt)) name_core_long2(fliplr(ind_disp_alt)) [3 9.5 13 17.5 21 25.5 27.5 32 38.5 43.75 46.5] [0 1 1 0 0 0 0 0 0 1 1]};
    for ii = 1:length(tmp{1})
        switch tmp{4}(ii)
            case 0
                text(tmp{3}(ii), 7e-8, tmp{2}(ii), 'fontsize', 20, 'color', get(tmp{1}(ii), 'markerfacecolor'), 'fontweight', 'bold')
            case 1
                [tmp2, tmp3]= strtok(tmp{2}{ii});
                text(tmp{3}(ii), 6e-8, {tmp2; tmp3(2:end)}, 'fontsize', 20, 'color', get(tmp{1}(ii), 'markerfacecolor'), 'fontweight', 'bold')
        end
    end
%     lg2                     = legend(pdu2(fliplr(ind_disp_alt)), name_core_long2(fliplr(ind_disp_alt)), 'location', 'south');
%     set(lg2, 'fontsize', 18, 'orientation', 'horizontal', 'position', [0.0050 0.0063 0.9857 0.0294])
    box on
    annotation('arrow', [0.06 0.20], [0.54 0.54], 'color', 'k', 'linewidth', 2, 'headstyle', 'plain', 'linestyle', '--')
    text(2, 1.7e-7, 'Increasing depth', 'color', 'k', 'fontsize', 20)
    annotation('arrow', [0.24 0.40], [0.98 0.98], 'color', 'k', 'linewidth', 2, 'headstyle', 'plain', 'linestyle', '--')
    text(20, 5.5e-5, 'Increasing mean annual surface temperature', 'color', 'k', 'fontsize', 20)    
    subplot('position', [0.055 0.015 0.89 0.44])
    hold on
    axis([0 (num_rlx + 1) 6e-8 1e-4])
    fill([0 0 (num_rlx + 1) (num_rlx + 1)], [6e-8 6e-7 6e-7 6e-8], [0.9 0.9 0.9], 'edgecolor', 'none')
    tmp                     = rlx.core{1};
    for ii = 1:num_rlx
        if ~strcmp(rlx.core{ii}, tmp)
            plot(repmat((ii - 0.5), 1, 2), get(gca, 'ylim'), 'k', 'linewidth', 1)
        end
        tmp                 = rlx.core{ii};
    end
    for ii = ind_relax_1'
        p_tmp               = plot(ii, (1e-6 .* rlx.H_1(ii)), 'ko', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
        if any(ii == ind_acc)
            set(p_tmp, 'marker', 'd')
        end
    end
    for ii = ind_relax_2'
        plot([ii ii], (1e-6 .* [rlx.H_1(ii), rlx.H_2(ii)]), 'k--', 'linewidth', 3)
        p_tmp               = plot(ii, (1e-6 .* rlx.H_1(ii)), 'k^', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
        p_tmp               = plot(ii, (1e-6 .* rlx.H_2(ii)), 'kv', 'markerfacecolor', color_core_alt(rlx.ind_core(ii), :), 'markersize', 12, 'linewidth', 1.5);
        if any(ii == ind_firn)
            set(p_tmp, 'color', color_core_alt(rlx.ind_core(ii), :), 'markerfacecolor', 'w')
        end
    end
    pdu                     = NaN(1, 5);
    pdu(1)                  = plot(NaN, NaN, 'k^', 'markersize', 12, 'markerfacecolor', [0.75 0.75 0.75], 'linewidth', 1.5);
    pdu(2)                  = plot(NaN, NaN, 'kv', 'markersize', 12, 'markerfacecolor', [0.75 0.75 0.75], 'linewidth', 1.5);
    pdu(3)                  = plot(NaN, NaN, 'ko', 'markersize', 12, 'markerfacecolor', [0.75 0.75 0.75], 'linewidth', 1.5);
    pdu(4)                  = plot(NaN, NaN, 'ko', 'markersize', 12, 'linewidth', 1.5);
    pdu(5)                  = plot(NaN, NaN, 'kd', 'markersize', 12, 'linewidth', 1.5);
    set(gca, 'fontsize', 20, 'xtick', [], 'yscale', 'log', 'ytick', logspace(-7, -4, 4), 'ygrid', 'on', 'layer', 'top')
    text(0.5, 6e-5, 'b', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold')
    text(1, 1.3e-7, '"Pure" relaxation', 'color', [0.4 0.4 0.4], 'fontsize', 20, 'fontweight', 'bold')
    ylabel('Apparent lattice [H^+] (M)')
    lg1                     = legend(pdu, 'Fast relaxation', 'Slow relaxation', 'Single relaxation', 'Firn', 'Accreted', 'location', 'southeast');
    set(lg1, 'fontsize', 18, 'orientation', 'horizontal', 'position', [0.485 0.53 0.4999 0.0294])
    axes('position', get(gca, 'position'), 'color', 'none', 'fontsize', 20, 'xaxislocation', 'top', 'yaxislocation', 'right', 'xtick', [], 'xscale', 'log', 'yscale', 'log')
    axis([0 (num_rlx + 1) ([6e-8 1e-4] ./ 0.48947)])
    ylabel('Apparent lattice [Cl^-] (M)', 'rotation', 270)
    
%% **** OLD FIGURE 3: TEMPERATURE DEPENDENCE OF HF CONDUCTIVITY AND RELAXATION FREQUENCIES ****** OLD

    color_relax             = {'k' 'r' 'b'};
    marker_relax            = {'o' '^' 'v'};
    figure('position', [2600 200 800 1200])
    subplot('position', [0.11 0.55 0.8 0.4])
    axis([(1e3 / 273.15) (1e3 / (273.15 - 90)) 1e-8 1e-4])
    hold on
    phf                     = NaN(1, 3);
    for ii = 1:num_core_fig3
        conduct_hf          = NaN(core_fig3{ii}.num_temp, 2);
        for jj = 1:core_fig3{ii}.num_temp
            num_relax_curr  = core_fig3{ii}.num_relax(jj);
            freq_curr       = core_fig3{ii}.freq_trim{jj};
            if core_fig3{ii}.do_dc(jj)
                model       = NaN((2 + (3 * num_relax_curr)), 1);
                model(end)  = -log10(core_fig3{ii}.conduct_dc(jj));
            else
                model       = NaN((1 + (3 * num_relax_curr)), 1);
            end
            try
                model(1)    = core_fig3{ii}.permitt_hf_model(jj);
            catch
                model(1)    = core_fig3{ii}.permitt_hf(jj);
            end
            for kk = 1:num_relax_curr
                [model(2 + ((kk - 1) * 3)), model(3 + ((kk - 1) * 3)), model(4 + ((kk - 1) * 3))] ...
                            = deal(core_fig3{ii}.permitt_diff_cat(kk, jj), (1 / (2 * pi * core_fig3{ii}.freq_relax_cat(kk, jj))), core_fig3{ii}.alpha_cole_cat(kk, jj));
            end
            if any(isnan(model))
                continue
            end
            permitt_complex = fm_cole_sep(model, freq_curr, num_relax_curr, core_fig3{ii}.do_dc(jj), permitt_vacuum);
            for kk = 1:min([2 num_relax_curr])
                conduct_hf(jj, kk) ...
                            = interp1(flipud(freq_curr), ((2 * pi * permitt_vacuum) .* flipud(freq_curr) .* flipud(permitt_complex{1 + kk}((length(freq_curr) + 1):end))), 3e5, 'nearest', 'extrap');
            end
        end
        if (all(core_fig3{ii}.num_relax == 1) || (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1))
            phf(1)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), conduct_hf(:, 1), 'color', [0.75 0.75 0.75], 'marker', marker_relax{1}, 'markerfacecolor', color_relax{1}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
        else
            phf(2)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), conduct_hf(:, 1), 'k', 'marker', marker_relax{2}, 'markerfacecolor', color_relax{2}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
            phf(3)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), conduct_hf(:, 2), 'k', 'marker', marker_relax{3}, 'markerfacecolor', color_relax{3}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
        end
    end
    set(gca, 'fontsize', 20, 'yscale', 'log', 'yticklabel', {'10^{-8}' '10^{-7}' '10^{-6}' '10^{-5}' ''})
    xlabel('1000 / Temperature (K)')
    ylabel('HF conductivity (S m^{-1})')
    lg                      = legend(phf, 'single relaxation', 'fast relaxation', 'slow relaxation', 'location', 'northeast');
    set(lg, 'fontsize', 20)
    box off
    text(3.55, 1e-4, 'a', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold')
    ax                      = gca;
    axes('position', get(gca, 'position'), 'color', 'none', 'fontsize', 20, 'xaxislocation', 'top', 'yaxislocation', 'right', 'xlim', get(gca, 'xlim'), 'yscale', 'log')
    xlabel('Temperature (\circC)')
    ylabel('Resistivity ({\Omega}{\cdot}m)', 'rotation', 270, 'position', [5.65 1e6 1.00005])
    ylim(fliplr(1 ./ get(ax, 'ylim')))
    set(gca, 'xtick', fliplr(temp_inv_vec), 'xticklabel', fliplr(temp_vec), 'ydir', 'reverse', 'ytick', logspace(4, 9, 6))
    subplot('position', [0.11 0.05 0.8 0.4])
    axis([(1e3 / 273.15) (1e3 / (273.15 - 90)) 1e0 1e5])
    hold on
    prfo                    = NaN(1, 2);
    for ii = 1:size(plot_data, 1)
        plot((1e3 ./ (plot_data{ii, 4} + 273.15)), exp(polyval(activ_poly{plot_data{ii, 5}}{plot_data{ii, 7}, plot_data{ii, 6}}, (1e3 ./ (plot_data{ii, 4} + 273.15)))), '--', 'linewidth', 2, 'color', plot_data{ii, 3})
        ind_curr            = find((temp_inv{plot_data{ii, 5}} <= (1e3 ./ (plot_data{ii, 4}(1) + 273.15))) & (temp_inv{plot_data{ii, 5}} >= (1e3 ./ (plot_data{ii, 4}(2) + 273.15))));
        if ~isempty(freq_relax_min{plot_data{ii, 5}})
            errorbar(temp_inv{plot_data{ii, 5}}(ind_curr), freq_relax{plot_data{ii, 5}}(plot_data{ii, 6}, ind_curr), (freq_relax{plot_data{ii, 5}}(plot_data{ii, 6}, ind_curr)' - freq_relax_min{plot_data{ii, 5}}{plot_data{ii, 6}}(ind_curr)), ...
                     (freq_relax_max{plot_data{ii, 5}}{plot_data{ii, 6}}(ind_curr) - freq_relax{plot_data{ii, 5}}(plot_data{ii, 6}, ind_curr)'), 'color', plot_data{ii, 3}, 'linestyle', 'none');
        end
        prfo(ii)            = plot(temp_inv{plot_data{ii, 5}}(ind_curr), freq_relax{plot_data{ii, 5}}(plot_data{ii, 6}, ind_curr), 'k', 'marker', plot_data{ii, 2}, 'markerfacecolor', plot_data{ii, 3}, 'markersize', 10, 'linestyle', 'none', 'linewidth', 1);
    end
    for ii = 1:num_core_fig3
        if (all(core_fig3{ii}.num_relax == 1) || (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1))
            prf(1)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), 'k', 'marker', marker_relax{1}, 'markerfacecolor', color_relax{1}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), (core_fig3{ii}.freq_relax_cat(1, :) - core_fig3{ii}.freq_relax_std_cat{1}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{1}(:, 2)' - core_fig3{ii}.freq_relax_cat(1, :)), 'ko', 'linestyle', 'none')
        else
            prf(2)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), 'k', 'marker', marker_relax{2}, 'markerfacecolor', color_relax{2}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
            prf(3)          = semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(2, :), 'k', 'marker', marker_relax{3}, 'markerfacecolor', color_relax{3}, 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(1, :), (core_fig3{ii}.freq_relax_cat(1, :) - core_fig3{ii}.freq_relax_std_cat{1}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{1}(:, 2)' - core_fig3{ii}.freq_relax_cat(1, :)), color_relax{2}, ...
                     'linestyle', 'none')
            errorbar((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(2, :), (core_fig3{ii}.freq_relax_cat(2, :) - core_fig3{ii}.freq_relax_std_cat{2}(:, 1)'), (core_fig3{ii}.freq_relax_std_cat{2}(:, 2)' - core_fig3{ii}.freq_relax_cat(2, :)), color_relax{3}, ...
                     'linestyle', 'none')
        end
        if (core_fig3{ii}.freq_relax_cat(2, interp1(core_fig3{ii}.temp_mean, 1:core_fig3{ii}.num_temp, -40, 'nearest', 'extrap')) < 1e1)
            for jj = 2:core_fig3{ii}.num_relax
                semilogy((1e3 ./ (273.15 + core_fig3{ii}.temp_mean)), core_fig3{ii}.freq_relax_cat(jj, :), 'k', 'marker', 'p', 'markerfacecolor', [0.75 0.75 0.75], 'markersize', 10, 'linewidth', 1);%, 'linestyle', 'none');
            end
        end
    end
    set(gca, 'fontsize', 20, 'yscale', 'log', 'ydir', 'reverse', 'yticklabel', {'' '10^1' '10^2' '10^3' '10^4' '10^5'})
    xlabel('1000 / Temperature (K)')
    ylabel('Relaxation frequency (Hz)')
    lg                      = legend(prfo([1:2 4]), plot_data([1:2 4], 1), 'location', 'northwest');
    set(lg, 'fontsize', 20)    
    box off
    text(3.55, 1e0, 'b', 'color', 'k', 'fontsize', 22, 'fontweight', 'bold')
    ax                      = gca;
    axes('position', get(gca, 'position'), 'color', 'none', 'fontsize', 20, 'xaxislocation', 'top', 'yaxislocation', 'right', 'xlim', get(gca, 'xlim'), 'yscale', 'log')
    xlabel('Temperature (\circC)')
    ylabel('Relaxation time (s)', 'rotation', 270, 'position', [5.65 2e-4 1.00005])
    ylim(fliplr(1 ./ ((2 * pi) .* get(ax, 'ylim'))))
    set(gca, 'xtick', fliplr(temp_inv_vec), 'xticklabel', fliplr(temp_vec))
    
%%
end
